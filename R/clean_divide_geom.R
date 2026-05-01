#' Clean and normalize divide geometries (POLYGON-only, valid, simplified)
#'
#' @description
#' Standardizes and repairs per-flowpath divide polygons so each
#' `member_COMID` maps to a single, valid `POLYGON` in EPSG:5070.
#' Handles `GEOMETRYCOLLECTION` (extracts polygons), `MULTIPOLYGON`
#' (keeps the largest part per `member_COMID`), fixes invalid geometry,
#' optionally dissolves/splits extras by `ID`, re-attaches attributes,
#' and (optionally) simplifies geometry while preserving shapes.
#'
#' If multiple polygon parts exist for a given divide, the function
#' selects a primary part (optionally informed by intersecting
#' `flowlines`) and intelligently re-associates smaller orphan parts.
#'
#' @param divides `sf` object of divide geometries in any CRS with at least:
#'   - `member_COMID` (identifier used to track completeness),
#'   - `ID` (grouping key used for dissolve/union steps),
#'   and a `POLYGON`/`MULTIPOLYGON`/`GEOMETRYCOLLECTION` geometry column.
#' @param ... (currently unused) reserved for future options.
#' @param keep Numeric in (0, 1] or `NULL`. Fraction of vertices to retain
#'   during simplification (passed to `rmapshaper::ms_simplify(keep=...)`).
#'   Use `NULL` to skip simplification. Default: `0.2`.
#' @param flowlines Optional `sf` LINESTRING layer (same features domain) with
#'   column `flowpath_id`. When supplied, intersect tests help choose a
#'   “prime” polygon among multiple parts for a given `ID`.
#'
#' @return An `sf` object in EPSG:5070 with valid `POLYGON` geometries:
#'   - If multiple parts per `ID` were present, the result includes the
#'     cleaned geometry with original attributes re-joined and an
#'     `areasqkm` column (via `hydrofab::add_areasqkm()`).
#'   - If only single-part polygons were present, returns the input
#'     attributes with cleaned geometry (and drops internal helper fields
#'     like `n`, `tmpID` if they exist).
#'
#' @details
#' **Workflow:**
#' 1. Transform to EPSG:5070; tag geometry type.
#' 2. Extract polygons from `GEOMETRYCOLLECTION`; pick largest part from
#'    `MULTIPOLYGON` per `member_COMID`.
#' 3. Snap to grid (tiny jitter) and validate; repair invalid features
#'    with `sf::st_make_valid()`.
#' 4. If any `member_COMID` has >1 polygon, dissolve by `ID`, explode,
#'    optionally mark a primary part using intersects with `flowlines`,
#'    and reconnect smaller parts to the most plausible main polygon.
#' 5. Optionally simplify with `ms_simplify(keep_shapes = TRUE)`.
#' 6. Guarantee no `member_COMID` loss; return cleaned geometry and
#'    attributes (plus `areasqkm`).
#'
#' **Assumptions:**
#' - `divides` contains `member_COMID` and `ID`.
#' - When provided, `flowlines` contains `flowpath_id` and is spatially
#'   compatible (will be transformed internally).
#'
#' @section CRS:
#' Output is returned in EPSG:5070.
#'
#' @examples
#' \dontrun{
#' # Minimal usage
#' clean <- clean_divide_geometry(divides)
#'
#' # Skip simplification
#' clean_nosimp <- clean_divide_geometry(divides, keep = NULL)
#'
#' # Use flowlines to guide "prime" part selection
#' clean_guided <- clean_divide_geometry(divides, flowlines = fl, keep = 0.3)
#' }
#'
#' @seealso
#' [rmapshaper::ms_simplify()], [rmapshaper::ms_dissolve()],
#' [rmapshaper::ms_explode()], [lwgeom::st_snap_to_grid()],
#' [sf::st_make_valid()], [sf::st_is_valid()]
#'
#' @importFrom sf st_transform st_geometry_type st_as_sf st_collection_extract
#' @importFrom sf st_cast st_area st_make_valid st_is_valid st_geometry
#' @importFrom sf st_zm st_crs st_intersects st_length st_drop_geometry
#' @importFrom sf st_intersection
#' @importFrom lwgeom st_snap_to_grid
#' @importFrom rmapshaper ms_dissolve ms_explode ms_simplify
#' @importFrom dplyr as_tibble mutate filter select slice_max bind_rows
#' @importFrom dplyr add_count row_number arrange desc left_join group_by
#' @importFrom dplyr summarise ungroup anti_join any_of
#' @importFrom tidyr replace_na
#' @importFrom hydrofab add_areasqkm union_polygons
#' @export

clean_divide_geometry <- function(divides, ..., keep = 0.2, flowlines = NULL) {
  # original_rows <- nrow(divides)
  original_comids <- unique(divides$member_COMID)

  divides <-
    divides |>
    sf::st_transform(5070) |>
    dplyr::as_tibble() |>
    dplyr::mutate(geometry_type = sf::st_geometry_type(geometry))

  # Handle GEOMETRYCOLLECTION
  collections <-
    dplyr::filter(divides, geometry_type == "GEOMETRYCOLLECTION")

  if (nrow(collections) > 0) {
    collections <-
      collections |>
      sf::st_as_sf(crs = 5070) |>
      sf::st_collection_extract("POLYGON", warn = FALSE) |>
      dplyr::as_tibble()
  }

  # Handle MULTIPOLYON
  multis <-
    dplyr::filter(divides, geometry_type == "MULTIPOLYGON")

  if (nrow(multis) > 0) {
    multis <-
      multis |>
      sf::st_as_sf(crs = 5070) |>
      sf::st_cast("POLYGON", warn = FALSE) |>
      dplyr::as_tibble() |>
      dplyr::mutate(area = sf::st_area(geometry)) |>
      dplyr::slice_max(
        order_by = area,
        n = 1,
        by = member_COMID,
        with_ties = FALSE
      ) |>
      dplyr::select(-area)
  }

  polygons <-
    dplyr::filter(divides, geometry_type == "POLYGON") |>
    dplyr::bind_rows(collections) |>
    dplyr::bind_rows(multis) |>
    sf::st_as_sf(crs = 5070) |>
    lwgeom::st_snap_to_grid(size = 0.0009) |>
    dplyr::select(-geometry_type)

  rm(collections, multis)

  is_invalid <-
    sf::st_geometry(polygons) |>
    # Don't remove this call!
    # For some reason, sf::st_is_valid is throwing
    # a segfault when assigning Rcpp attributes.
    # This occurs in CPL_geos_is_valid as of 6/20/25.
    sf::st_make_valid() |>
    sf::st_is_valid() |>
    which() |>
    setdiff(seq_len(nrow(polygons)))

  if (length(is_invalid) > 0) {
    sf::st_geometry(polygons)[is_invalid] <-
      sf::st_geometry(polygons)[is_invalid] |>
      sf::st_make_valid() |>
      sf::st_cast("POLYGON")
  }

  polygons <-
    polygons |>
    dplyr::add_count(member_COMID) |>
    dplyr::mutate(
      areasqkm = hydrofab::add_areasqkm(geometry),
      tmpID = dplyr::row_number()
    )

  if (any(sf::st_geometry_type(polygons) != "POLYGON")) {
    # separate polygons with more than 1 feature counts
    # dissolve, and explode if necessary
    extra_parts <-
      dplyr::filter(polygons, n > 1)

    # if (!is.null(flowlines)) ...

    extra_parts <-
      extra_parts |>
      rmapshaper::ms_dissolve(
        field = "ID",
        copy_fields = names(extra_parts),
        sys = TRUE
      ) |>
      rmapshaper::ms_explode(sys = TRUE) |>
      dplyr::mutate(tmpID = dplyr::row_number())

    if (!is.null(flowlines)) {
      flowlines <-
        dplyr::filter(flowlines, flowpath_id %in% unique(extra_parts$ID)) |>
        sf::st_zm() |>
        sf::st_transform(sf::st_crs(extra_parts))

      imap <- sf::st_intersects(extra_parts, flowlines)
      l <- lengths(imap)

      df <-
        data.frame(
          tmpID = rep(extra_parts$tmpID, times = l),
          uid = rep(extra_parts$ID, times = l),
          touch_id = flowlines$flowpath_id[unlist(imap)]
        ) |>
        dplyr::group_by(tmpID) |>
        dplyr::summarise(prime = any(uid == touch_id))

      extra_parts <-
        dplyr::left_join(extra_parts, df, by = "tmpID") |>
        dplyr::mutate(prime = tidyr::replace_na(prime, FALSE))

      rm(flowlines, imap, l, df)
    } else {
      extra_parts$prime <- FALSE
    }

    extra_parts <-
      extra_parts |>
      dplyr::mutate(areasqkm = hydrofab::add_areasqkm(geometry)) |>
      dplyr::arrange(dplyr::desc(prime), dplyr::desc(areasqkm)) |>
      dplyr::mutate(newID = dplyr::row_number()) |>
      dplyr::as_tibble()

    main_parts <-
      dplyr::slice(extra_parts, 1, .by = ID)

    small_parts <-
      dplyr::anti_join(extra_parts, main_parts, by = "newID")

    # if(!sum(nrow(main_parts)) + nrow(filter(polygons, n == 1)) == MASTER_COUNT){
    #   stop()
    # }

    if (
      !setequal(
        unique(c(
          dplyr::filter(polygons, n == 1)$member_COMID,
          main_parts$member_COMID
        )),
        original_comids
      )
    ) {
      stop("lost COMIDs")
    }

    main_parts <-
      dplyr::bind_rows(
        main_parts,
        dplyr::filter(polygons, n == 1)
      )

    if (!is.null(small_parts) && nrow(small_parts) > 0) {
      small_parts <-
        small_parts |>
        sf::st_as_sf(crs = 5070) |>
        rmapshaper::ms_dissolve(
          field = "ID",
          copy_fields = names(small_parts),
          sys = TRUE
        ) |>
        rmapshaper::ms_explode(sys = TRUE) |>
        dplyr::mutate(
          areasqkm = hydrofab::add_areasqkm(geometry),
          newID = dplyr::row_number()
        ) |>
        dplyr::select(newID, geometry)

      out <-
        main_parts |>
        sf::st_as_sf(crs = 5070) |>
        sf::st_make_valid() |>
        sf::st_intersection(x = small_parts, warn = FALSE) |>
        sf::st_collection_extract("LINESTRING", warn = FALSE) |>
        suppressWarnings()

      ints <-
        out |>
        dplyr::mutate(l = sf::st_length(geometry)) |>
        dplyr::slice_max(l, with_ties = FALSE, by = newID)

      tj <-
        ints |>
        sf::st_drop_geometry() |>
        dplyr::select(ID, newID) |>
        dplyr::left_join(small_parts, by = "newID") |>
        dplyr::bind_rows(main_parts) |>
        dplyr::select(-areasqkm, -tmpID, -newID) |>
        dplyr::group_by(ID) |>
        dplyr::mutate(n = dplyr::n()) |>
        dplyr::ungroup() |>
        sf::st_as_sf(crs = 5070) |>
        hfutils::rename_geometry("geometry")

      in_cat <-
        dplyr::filter(tj, n > 1) |>
        hydrofab::union_polygons("ID") |>
        dplyr::bind_rows(
          dplyr::filter(tj, n == 1) |>
            dplyr::select("ID")
        ) |>
        dplyr::mutate(tmpID = dplyr::row_number())
    } else {
      in_cat <- main_parts
    }

    in_cat_invalid <-
      sf::st_geometry(in_cat) |>
      sf::st_make_valid() |>
      sf::st_is_valid() |>
      which() |>
      setdiff(seq_len(nrow(in_cat)))

    if (length(in_cat_invalid) > 0) {
      sf::st_geometry(in_cat)[in_cat_invalid] <-
        sf::st_geometry(in_cat)[in_cat_invalid] |>
        sf::st_make_valid() |>
        sf::st_cast("POLYGON")
    }

    if (!is.null(keep)) {
      in_cat <-
        rmapshaper::ms_simplify(
          in_cat,
          keep = keep,
          sys = TRUE,
          keep_shapes = TRUE
        )
    }

    x <-
      divides |>
      sf::st_drop_geometry() |>
      dplyr::mutate(ID = as.numeric(ID))

    x2 <-
      in_cat |>
      dplyr::mutate(areasqkm = hydrofab::add_areasqkm(geometry)) |>
      sf::st_transform(5070) |>
      dplyr::mutate(ID = as.numeric(ID)) |>
      dplyr::select(ID, areasqkm) |>
      dplyr::left_join(x, by = "ID")

    return(x2)
  } else {
    if (!is.null(keep)) {
      polygons <-
        rmapshaper::ms_simplify(
          polygons,
          keep = keep,
          sys = TRUE,
          keep_shapes = TRUE
        )
    }

    return(dplyr::select(polygons, -dplyr::any_of(c("n", "tmpID"))))
  }
}
