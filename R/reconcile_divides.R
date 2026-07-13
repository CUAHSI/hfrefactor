#' Reconcile refactored flowpaths with catchment divide polygons
#'
#' @description
#' Produces a reconciled set of catchment *divides* (MULTIPOLYGON, EPSG:5070)
#' consistent with a refactored flowpath network. Flowpaths that were split
#' during refactoring (e.g., have dotted `COMID` such as `"12345.1"`) trigger
#' geometric subdivision of their parent divide using flow-direction (FDR) and
#' flow-accumulation (FAC) rasters. Unsplit divides are passed through.
#' When multiple refactored flowpaths map to one reconciled ID, the function
#' unions the corresponding pieces. A log of per-FEATUREID outcomes is written
#' under `/tmp/reconcile_divides/`.
#'
#' @param divides `sf` polygon layer of original catchment divides; must contain
#'   `FEATUREID` (numeric or numeric-like character). Geometry can be any
#'   polygonal type; will be transformed/cast as needed.
#' @param refactored_fp `sf` line layer of **refactored** flowpaths with a
#'   `COMID` column (character). Dotted `COMID` values (e.g., `"1001.2"`)
#'   indicate splits of a base feature id (`1001`).
#' @param reconciled_fp `sf` (any geometry) containing the mapping between
#'   `ID` (reconciled divide identifier) and `member_COMID` (comma-separated
#'   list of member flowpath COMIDs). Only the columns `ID` and `member_COMID`
#'   are used.
#' @param fdr Flow direction raster **or** a path readable by
#'   [terra::rast()]. CRS may differ; it is reprojected/cropped internally.
#' @param fac Flow accumulation raster **or** a path readable by
#'   [terra::rast()]. CRS may differ; it is reprojected/cropped internally.
#' @param dem_vrt Optional path to a DEM VRT/TIFF used when computing FDR/FAC
#'   on the fly (WhiteboxTools). In the current configuration (`use_dem = FALSE`)
#'   this is not used; kept for parity and future toggling.
#' @param min_area_m Numeric. Minimum polygon area (m²) to retain during split
#'   operations (passed through to the splitting helper).
#' @param snap_distance_m Numeric. Snapping distance (meters) used by the
#'   splitting helper when aligning boundaries to flowlines/rasters.
#' @param simplify_tolerance_m Numeric. Douglas–Peucker tolerance (meters)
#'   applied by the splitting helper for output smoothing.
#' @param keep Numeric in (0,1] or `NULL`. Intended fraction of vertices to
#'   retain for a final simplification pass. **Currently unused** in this
#'   function; reserved for future refinement to mirror other cleaners.
#'
#' @return An `sf` MULTIPOLYGON layer in EPSG:5070 with columns:
#' \itemize{
#'   \item `ID` — reconciled divide identifier (from `reconciled_fp`)
#'   \item `member_COMID` — comma-separated COMIDs associated to `ID`
#'   \item `geometry` — reconciled divide geometry (MULTIPOLYGON)
#' }
#' The result contains:
#' \enumerate{
#'   \item divides that required splitting (derived via FDR/FAC),
#'   \item divides that did not require splitting (passed through),
#'   \item divides that were unioned where a single `ID` has multiple members.
#' }
#'
#' @details
#' **High-level workflow**
#' \enumerate{
#'   \item Normalize `reconciled_fp` to a distinct mapping of `ID` → `member_COMID`.
#'   \item Keep only refactored flowpaths whose base `COMID` exists in `divides`.
#'   \item Identify base FEATUREIDs that were split (COMIDs like `"x.y"`).
#'   \item For each such FEATUREID, crop/mask FDR/FAC (or compute from `dem_vrt`
#'         using WhiteboxTools if toggled), then invoke a splitting helper
#'         (`split_divide()`, falling back to `hydrofab:::split_catchment_divide()`).
#'   \item Cast to `MULTIPOLYGON` and de-duplicate per FEATUREID.
#'   \item Collect divides that needed no split and still match refactored flowpaths.
#'   \item Where `reconciled_fp` maps multiple `member_COMID` to one `ID`,
#'         union split pieces via `hydrofab:::get_split_cats()`.
#'   \item Join back to `reconciled_fp` so output carries both `ID` and `member_COMID`.
#' }
#'
#' **Logging**
#' A TSV log is appended at a timestamped path under
#' `/tmp/reconcile_divides/`. Each row captures start time, status,
#' FEATUREID, DEM source hint (NED vs GEDTM30), and whether a retry path
#' was used. This is helpful for post-hoc QA of problematic basins.
#'
#' **Assumptions & Notes**
#' * `FEATUREID` in `divides` corresponds to the **base** of dotted `COMID`
#'   values in `refactored_fp` (e.g., `FEATUREID == 1001` matches `"1001.1"`).
#' * Rasters may be reprojected/cropped/masked internally to a buffered area
#'   around the target divide.
#' * This function references non-exported helpers in **hydrofab** via
#'   `hydrofab:::` as a fallback path; these are expected to be available
#'   when running inside the package context.
#'
#' @examples
#' \dontrun{
#' out <- reconcile_divides(
#'   divides = divides_sf,
#'   refactored_fp = refactored_flowlines_sf,
#'   reconciled_fp = id_to_members_sf,  # columns: ID, member_COMID
#'   fdr = "/path/to/fdr.tif",
#'   fac = "/path/to/fac.tif",
#'   dem_vrt = "/path/to/dem.vrt",      # not used unless DEM path is enabled
#'   min_area_m = 2e4,
#'   snap_distance_m = 30,
#'   simplify_tolerance_m = 5,
#'   keep = NULL
#' )
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{split_divide}} — primary splitting helper used here.
#'   \item Non-exported fallbacks: \code{hydrofab:::split_catchment_divide},
#'         \code{hydrofab:::get_split_cats}.
#'   \item Raster ops: [terra::rast()], [terra::crop()], [terra::mask()].
#' }
#'
#' @importFrom sf st_drop_geometry st_as_sf st_set_crs st_transform st_cast
#' @importFrom sf st_is_empty st_geometry st_dimension
#' @importFrom terra rast crop vect buffer mask project crs
#' @importFrom whitebox wbt_breach_depressions_least_cost
#' @importFrom whitebox wbt_fill_depressions_wang_and_liu
#' @importFrom whitebox wbt_d8_flow_accumulation wbt_d8_pointer
#' @importFrom hfutils rename_geometry
#' @importFrom dplyr select distinct filter group_by summarise ungroup
#' @importFrom dplyr as_tibble bind_rows mutate slice_head left_join
#' @importFrom tidyr separate_rows
#' @importFrom logger log_info
#' @importFrom glue glue_data
#' @export

reconcile_divides <- function(
    divides,
    refactored_fp,
    reconciled_fp,
    fdr,
    fac,
    dem_vrt,
    min_area_m2,
    snap_distance_m,
    simplify_tolerance_m,
    keep,
    n_cores = 1L
) {
  reconciled <-
    reconciled_fp |>
    sf::st_drop_geometry() |>
    dplyr::select(ID, member_COMID)

  comid <- refactored_fp$COMID
  featureid <- divides$FEATUREID
  comid_with_catchment <- comid[trunc(as.numeric(comid)) %in% featureid]

  reconciled <-
    reconciled |>
    dplyr::distinct() |>
    tidyr::separate_rows(member_COMID, sep = ",") |>
    dplyr::filter(trunc(as.numeric(member_COMID)) %in% comid_with_catchment) |>
    dplyr::group_by(ID) |>
    dplyr::summarise(member_COMID = paste(member_COMID, collapse = ",")) |>
    dplyr::ungroup()

  refactored_fp <-
    dplyr::as_tibble(refactored_fp) |>
    dplyr::filter(trunc(as.numeric(COMID)) %in% featureid) |>
    sf::st_as_sf(crs = 5070)

  to_split_bool <-
    grepl(".", refactored_fp$COMID, fixed = TRUE)

  to_split_ids <-
    refactored_fp$COMID[which(to_split_bool)]

  to_split_featureids <-
    unique(trunc(as.numeric(to_split_ids)))

  logger::log_info(paste("catchments to split:", length(to_split_featureids)))

  # --- begin split_cat -------------------------------------------------------
  split_cat <- function(feature_id) {
    # message(glue::glue("{Sys.time()} split_cat: {feature_id}"))

    split_set <- to_split_ids[
      as.character(feature_id) |>
        startsWith(x = to_split_ids) |>
        which()
    ]

    to_split_fline <-
      dplyr::filter(refactored_fp, COMID %in% split_set)

    to_split_cat <-
      dplyr::filter(divides, FEATUREID == !!feature_id) |>
      sf::st_set_crs(5070)

    # Check if divide falls in the FDR
    cropped <-
      suppressWarnings(
        terra::crop(
          terra::rast(fdr),
          terra::vect(to_split_cat)
        )
      )


    # Use NED but fallback to GEDTM
    # ---
    #> use_dem <- is.null(tryCatch({
    #>   # Check if we can crop/mask to the divide
    #>   masked <- suppressWarnings(
    #>     terra::rast(fdr) |>
    #>     terra::crop(terra::vect(to_split_cat), snap = "out") |>
    #>     terra::mask(terra::vect(to_split_cat))
    #>   )
    #>
    #>   # Then, check that we actually have values and not just NaNs
    #>   if (all(terra::values(is.nan(masked) | is.na(masked)))) {
    #>     # If only NaNs, return NULL so we use the DEM
    #>     # Note: this usually happens when we're in the bounds of the
    #>     #       NED-derived FDR, but in an area where there's no values
    #>     NULL
    #>   } else {
    #>     # Otherwise, we have values, so use NED
    #>     logical(0)
    #>   }
    #> },
    #>   error = function(e) {
    #>     if (grepl("[crop] too few values for writing", conditionMessage(e), fixed = TRUE)) {
    #>       NULL
    #>     } else {
    #>       conditionCall(e)
    #>     }
    #>   }
    #> ))

    # Always use GEDTM
    # ---
    #> use_dem <- TRUE

    # Always use NED
    # ---
    use_dem <- FALSE

    buffered_area <-
      terra::vect(to_split_cat) |>
      terra::buffer(250)

    # If we use the given DEM, we are computing
    # FAC/FDR ourselves, using whitebox.
    if (use_dem) {
      workdir <- file.path(tempdir(), as.character(to_split_cat$FEATUREID))
      dir.create(workdir)
      on.exit(unlink(workdir, TRUE, TRUE), TRUE)
      path_dem <- file.path(workdir, "dem.tif")
      path_breached <- file.path(workdir, "breached.tif")
      path_filled <- file.path(workdir, "filled.tif")
      path_fac <- file.path(workdir, "fac.tif")
      path_fdr <- file.path(workdir, "fdr.tif")

      dem <-
        terra::rast(dem_vrt)

      buffered_area <-
        buffered_area |>
        terra::project(terra::crs(dem))

      # Write out DEM
      terra::crop(dem, buffered_area, snap = "out") |>
        terra::mask(buffered_area) |>
        terra::writeRaster(path_dem)

      whitebox::wbt_breach_depressions_least_cost(
        dem = path_dem,
        output = path_breached,
        dist = 5,
        fill = TRUE,
        verbose_mode = FALSE
      )

      whitebox::wbt_fill_depressions_wang_and_liu(
        dem = path_breached,
        output = path_filled,
        verbose_mode = FALSE
      )

      whitebox::wbt_d8_flow_accumulation(
        input = path_filled,
        output = path_fac,
        verbose_mode = FALSE
      )

      whitebox::wbt_d8_pointer(
        dem = path_filled,
        output = path_fdr,
        verbose_mode = FALSE
      )

      fdr <- terra::rast(path_fdr)
      fac <- terra::rast(path_fac)
    } else {
      fdr <- terra::rast(fdr)
      fac <- terra::rast(fac)
    }

    buffered_area <-
      buffered_area |>
      terra::project(terra::crs(fdr))

    # Ensure our masked areas are in the projection of the divide
    fdr <-
      terra::crop(fdr, buffered_area, snap = "out") |>
      terra::mask(buffered_area) |>
      terra::project(terra::crs(to_split_cat))

    fac <-
      terra::crop(fac, buffered_area, snap = "out") |>
      terra::mask(buffered_area) |>
      terra::project(terra::crs(to_split_cat))



    # split_cats <- hydrofab:::split_catchment_divides(...)
    # -----------------------------------------------------
    first_error <- NULL
    split_cats <- try(
      split_divide(
        feature_id,
        to_split_cat,
        to_split_fline,
        fdr,
        fac,
        min_area_m2,
        snap_distance_m,
        simplify_tolerance_m
      ),
      silent = TRUE
    )

    if (inherits(split_cats, "try-error")) {
      first_error <- split_cats
      split_cats <- try(
        hydrofab:::split_catchment_divide(
          catchment = to_split_cat,
          fline = to_split_fline,
          fdr = fdr,
          fac = fac,
          lr = FALSE,
          min_area_m = min_area_m2,
          snap_distance_m = snap_distance_m,
          simplify_tolerance_m = simplify_tolerance_m,
          vector_crs = "EPSG:5070"
        ),
        silent = TRUE
      )

      if (inherits(split_cats, "try-error")) {
        second_error <- split_cats
        stop(
          "splitting failed with both methods:",
          "\n(1)", first_error,
          "\n(2)", second_error
        )
      }
    }

    # -----------------------------------------------------

    sf::st_sf(
      FEATUREID = to_split_fline$COMID,
      retried = !is.null(first_error),
      gedtm30 = use_dem,
      geometry = split_cats |>
        sf::st_transform("EPSG:5070") |>
        sf::st_cast("MULTIPOLYGON")
    )
  }
  # --- end split_cat ---------------------------------------------------------

  log_path <- paste0(
    strftime(Sys.time(), "/tmp/reconcile_divides/%Y%m%dT%H%M%S"), ".log"
  )

  if (!dir.exists(dirname(log_path))) {
    dir.create(dirname(log_path), recursive = TRUE)
  }

  if (!file.exists(log_path)) {
    cat("time\tstatus\tfeatureid\tmessage\n", file = log_path)
  }

  log_fmt <- "{arg_time}\t{arg_status}\t{arg_fid}\t{arg_msg}\n"

  n_cores <- max(1L, as.integer(n_cores))

  split_cats <-
    parallel::mclapply(
      to_split_featureids,
      function(fid) {
        status <- "SUCCESS"
        start_time <- Sys.time()
        msg <- ""

        result <- try(split_cat(fid), silent = TRUE)
        if (inherits(result, "try-error")) {
          status <- "ERROR"
          msg <- as.character(result)[1]
          result <-
            dplyr::filter(divides, FEATUREID == !!fid) |>
            dplyr::select(FEATUREID, geometry) |>
            dplyr::mutate(FEATUREID = as.character(FEATUREID)) |>
            sf::st_set_crs(5070) |>
            sf::st_cast("MULTIPOLYGON")
        } else {
          if (result$gedtm30[1]) {
            msg <- paste(msg, "dem=GEDTM30")
          } else {
            msg <- paste(msg, "dem=Input")
          }

          result$gedtm30 <- NULL

          if (result$retried[1]) {
            msg <- paste(msg, "[retried]")
          }

          result$retried <- NULL
        }

        end_time <- Sys.time()
        diff_time <- format(end_time - start_time)
        attr(result, ".log_entry") <- list(
          arg_time   = start_time,
          arg_status = status,
          arg_fid    = fid,
          arg_msg    = trimws(paste0("(", diff_time, ") ", msg), "right")
        )

        result
      },
      mc.cores = n_cores
    )

  # Workers can't share a file connection, so write log entries after collection
  log_file <- file(log_path, "a")
  for (item in split_cats) {
    entry <- attr(item, ".log_entry")
    if (!is.null(entry)) {
      log_msg <- glue::glue_data(entry, log_fmt)
      cat(log_msg, "\n", file = log_file)
      message(log_msg)
    }
  }
  close(log_file)

  if (length(split_cats) == 0) {
    split_cats <-
      sf::st_sf(FEATUREID = NA, geometry = list(sf::st_multipolygon()))
  } else {
    split_cats <-
      split_cats |>
      Filter(f = function(x) !is.null(x)) |>
      dplyr::bind_rows() |>
      sf::st_as_sf() |>
      dplyr::filter(!sf::st_is_empty(geometry))
  }

  split_cats <-
    sf::st_cast(split_cats, "MULTIPOLYGON") |>
    dplyr::slice_head(n = 1, by = FEATUREID)

  # unsplit <- hydrofab:::get_cat_unsplit(
  #   divides, refactored_fp, to_split_featureids)
  # ----------------------------------------------

  # These are the divides that did not require splitting,
  # and that have a corresponding refactored flowpath.
  unsplit <-
    divides |>
    dplyr::select(FEATUREID) |>
    dplyr::filter(
      !FEATUREID %in% to_split_featureids,
      !FEATUREID %in% trunc(as.numeric(split_cats$FEATUREID))
    ) |>
    dplyr::mutate(FEATUREID = as.character(FEATUREID)) |>
    dplyr::filter(FEATUREID %in% refactored_fp$COMID) |>
    sf::st_transform(5070)

  # faster than sf::st_cast

  sf::st_geometry(unsplit) <-
    tryCatch(
      sf::st_geometry(unsplit) |>
        lapply(list) |>
        lapply(sf::st_multipolygon) |>
        sf::st_as_sfc(crs = sf::st_crs(unsplit)),
      error = function(err) {
        sf::st_geometry(unsplit) |>
          sf::st_cast("MULTIPOLYGON")
      })

  if (nrow(unsplit) > 0) {
    split_cats <-
      dplyr::bind_rows(unsplit, split_cats) |>
      sf::st_as_sf()
  }

  # List of reconciled divide IDs where there are more than 1 member COMIDs
  combinations <- reconciled$member_COMID[grepl(",", reconciled$member_COMID)]

  logger::log_info("checking for unioned cats")
  unioned_cats <- lapply(
    combinations,
    hydrofab:::get_split_cats,
    split_cats = split_cats,
    cache = NULL
  )

  if (length(unioned_cats) > 0) {
    unioned_cats <-
      unioned_cats |>
      dplyr::bind_rows() |>
      sf::st_as_sf() |>
      hfutils::rename_geometry("geometry")

    # At this stage, split_cats contains:
    # - The divides requiring a split
    # - The divides not requiring a split (unsplit)
    # - The divides not reqiring a split, but required unioning (unioned_cats)
    split_cats <-
      split_cats |>
      dplyr::filter(
        !FEATUREID %in% unioned_cats$FEATUREID,
        !FEATUREID %in%
          unique(unlist(strsplit(unioned_cats$FEATUREID, ",", fixed = TRUE)))
      ) |>
      dplyr::bind_rows(unioned_cats) |>
      sf::st_as_sf()
  }

  if (nrow(split_cats) > 0) {
    out <-
      dplyr::select(split_cats, member_COMID = FEATUREID) |>
      dplyr::left_join(reconciled, by = "member_COMID") |>
      sf::st_as_sf()

    missing <- is.na(sf::st_dimension(sf::st_geometry(out)))
    if (any(missing)) {
      out_mp <-
        dplyr::filter(out, !missing) |>
        sf::st_cast("MULTIPOLYGON")

      out <-
        dplyr::select(divides, member_COMID = FEATUREID) |>
        dplyr::filter(
          member_COMID %in% unique(floor(as.numeric(out$member_COMID[missing])))
        ) |>
        dplyr::mutate(member_COMID = paste0(member_COMID, ".1")) |>
        dplyr::mutate(ID = out$ID[match(member_COMID, out$member_COMID)]) |>
        dplyr::select(ID, member_COMID) |>
        hfutils::rename_geometry("geometry") |>
        dplyr::bind_rows(out_mp)
    }

    out
  }
}
# =============================================================================
#' Split a catchment divide into upstream/downstream pieces using FDR/FAC
#'
#' @description
#' Given a single catchment polygon (`to_split_cat`) and one or more refactored
#' flowline segments traversing it (`to_split_fline`), this function iteratively
#' carves the divide into **upstream** pieces for each intermediate outlet and a
#' final **downstream remainder**. It uses a pre-extracted flow-direction (FDR)
#' and flow-accumulation (FAC) raster window aligned to the catchment to
#' delineate contributing-area masks, then converts those masks to polygons and
#' performs geometric difference operations to avoid overlaps.
#'
#' @param feature_id Integer-like. Base FEATUREID of the catchment being split;
#'   used only for logging/error messages.
#' @param to_split_cat `sf` POLYGON/MULTIPOLYGON (EPSG:5070 preferred). The
#'   catchment/divide geometry to be split.
#' @param to_split_fline `sf` LINESTRING(s) crossing `to_split_cat`. The line
#'   endpoints (computed via `lwgeom::st_endpoint()`) are used as ordered
#'   outlets for successive splits.
#' @param fdr A flow-direction raster (or a path readable by `terra::rast()`),
#'   already cropped/projected to a neighborhood that covers `to_split_cat`.
#' @param fac A flow-accumulation raster (or a path readable by `terra::rast()`),
#'   co-registered with `fdr`. Must have identical dimensions after prep.
#' @param min_area_m2 Numeric. Minimum polygon area in **m²** to retain when
#'   extracting the downstream remainder (filters slivers).
#' @param snap_distance_m Numeric. Distance (meters) used to snap the
#'   upstream piece to the outlet point to ensure topological contact.
#' @param simplify_tolerance_m Numeric. Douglas–Peucker tolerance (meters)
#'   used to simplify the upstream polygon before differencing.
#'
#' @return An `sfc` (EPSG:5070) of polygonal geometries in order:
#'   all **upstream** pieces for each intermediate outlet, followed by the
#'   **downstream** remainder as the last element. Each element is a valid
#'   polygonal geometry; overall the vector may contain `POLYGON` or
#'   `MULTIPOLYGON` members depending on topology.
#'
#' @details
#' **Algorithm sketch**
#' \enumerate{
#'   \item Compute outlet coordinates from `to_split_fline` endpoints and verify
#'         each falls within `to_split_cat`.
#'   \item Prepare aligned FDR/FAC matrices via
#'         `hydrofab:::prep_cat_fdr_fac(to_split_cat, fdr, fac)`.
#'   \item For each outlet (except the terminal one), locate its raster cell
#'         (`hydrofab:::get_row_col()`), collect all upstream cells
#'         (`hydrofab:::collect_upstream()`), and rasterize a binary mask.
#'   \item Convert the mask to polygons, simplify (`simplify_tolerance_m`),
#'         snap to the outlet (`snap_distance_m`), and make valid.
#'   \item Difference the current `catchment` by the upstream piece to obtain
#'         a downstream remainder; drop tiny remnants `< min_area_m2`, and keep
#'         the largest polygonic remainder.
#'   \item Append the upstream piece to the output list, update `catchment` to
#'         the remainder, and continue to the next outlet.
#'   \item After all outlets are processed, append the final remainder.
#' }
#'
#' If at any step an outlet does not lie within the catchment, the function
#' logs an error and stops. A consistency check ensures FDR/FAC matrices share
#' identical dimensions.
#'
#' @examples
#' \dontrun{
#' pieces <- split_divide(
#'   feature_id = 1001,
#'   to_split_cat = cat_sf,      # single catchment polygon
#'   to_split_fline = flines_sf, # refactored segment(s) inside the catchment
#'   fdr = fdr_rast,             # terra SpatRaster or path
#'   fac = fac_rast,             # terra SpatRaster or path
#'   min_area_m2 = 2e4,
#'   snap_distance_m = 30,
#'   simplify_tolerance_m = 5
#' )
#' }
#'
#' @seealso
#' Non-exported helpers used internally:
#' \code{hydrofab:::prep_cat_fdr_fac},
#' \code{hydrofab:::get_row_col},
#' \code{hydrofab:::collect_upstream}.
#'
#' @importFrom lwgeom st_endpoint
#' @importFrom sf st_coordinates st_point st_sfc st_within st_as_sf st_transform
#' @importFrom sf st_simplify st_snap st_make_valid st_geometry st_difference
#' @importFrom sf st_cast st_combine st_is_empty
#' @importFrom terra setValues as.polygons
#' @importFrom dplyr as_tibble mutate filter slice_max
#' @importFrom units set_units
#' @importFrom logger log_error
#' @export

split_divide <- function(
    feature_id,
    to_split_cat,
    to_split_fline,
    fdr,
    fac,
    min_area_m2,
    snap_distance_m,
    simplify_tolerance_m
    ) {
  outlets <-
    lwgeom::st_endpoint(to_split_fline) |>
    sf::st_coordinates() |>
    dplyr::as_tibble() |>
    dplyr::mutate(L1 = dplyr::row_number())

  fdr_matrix <-
    hydrofab:::prep_cat_fdr_fac(to_split_cat, fdr, fac)

  fdr_cat <- fdr_matrix$fdr
  fac_cat <- fdr_matrix$fac
  fac_matrix <- fdr_matrix$fac_matrix
  fdr_matrix <- fdr_matrix$fdr_matrix

  if (any(dim(fac_matrix) != dim(fdr_matrix))) {
    stop("flow direction and flow accumulation must be the same size")
  }

  return_cats <- list()
  catchment <- to_split_cat
  smaller_than_one_pixel <- units::set_units(min_area_m2, "m^2")
  snap_distance <- units::set_units(snap_distance_m, "m")
  for (i in seq_len(nrow(outlets) - 1)) {
    in_out <-
      (as.matrix(outlets[i, 1:2]) |>
         sf::st_point(dim = "XY") |>
         sf::st_sfc(crs = "EPSG:5070") |>
         sf::st_within(to_split_cat, sparse = FALSE))[1, 1]

    if (length(in_out) == 0 || !in_out) {
      logger::log_error(sprintf(
        "not in_out <feature_id=%d, event_id=%s>",
        feature_id, to_split_fline$event_identifier[
          to_split_fline$event_identifier != ""
        ]
      )[1])
      stop("not in_out")
    }

    row_col <- hydrofab:::get_row_col(
      fdr_cat,
      start = as.matrix(outlets[i, 1:2]),
      fac_matrix
    )

    upstream_cells <-
      hydrofab:::collect_upstream(row_col, fdr_matrix)

    mask <-
      matrix(
        0,
        nrow = nrow(fdr_matrix),
        ncol = ncol(fdr_matrix)
      )

    mask[upstream_cells] <- 1

    out <- terra::setValues(fdr_cat, mask)
    names(out) <- "cats"
    out <-
      terra::as.polygons(out) |>
      sf::st_as_sf() |>
      dplyr::filter(cats == 1) |>
      sf::st_transform(5070)

    snap_outlet <-
      outlets[i, c("X", "Y")] |>
      sf::st_as_sf(coords = 1:2, crs = 5070)

    us_catchment <-
      sf::st_geometry(out) |>
      sf::st_simplify(dTolerance = simplify_tolerance_m) |>
      sf::st_snap(snap_outlet, tolerance = snap_distance) |>
      # Don't snap to catchment divide
      # sf::st_snap(sf::st_geometry(catchment), tolerance = snap_distance) |>
      sf::st_make_valid()

    ds_catchment <-
      catchment |>
      sf::st_geometry() |>
      sf::st_make_valid() |>
      sf::st_difference(us_catchment) |>
      sf::st_cast("POLYGON") |>
      sf::st_as_sf(crs = 5070) |>
      dplyr::mutate(area = sf::st_area(x)) |>
      dplyr::filter(area > smaller_than_one_pixel) |>
      # Added this:
      dplyr::slice_max(order_by = area, n = 1) |>
      sf::st_combine() |>
      sf::st_make_valid()

    # This ensures we have no overlapping upstream
    # catchment area.
    us_catchment <-
      catchment |>
      sf::st_geometry() |>
      sf::st_make_valid() |>
      sf::st_difference(ds_catchment) |>
      sf::st_transform(5070) |>
      sf::st_make_valid()

    catchment <-
      ds_catchment

    return_cats <- c(return_cats, us_catchment)
  }

  if (sf::st_is_empty(sf::st_geometry(catchment))) {
    stop("nothing left over. split too small??")
  }

  sf::st_sfc(
    c(return_cats, sf::st_geometry(catchment)),
    crs = "EPSG:5070"
  )
}
