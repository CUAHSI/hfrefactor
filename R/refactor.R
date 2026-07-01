#' Prepare inputs for flowline refactoring and divide reconciliation
#'
#' @description
#' Normalizes column names and geometries for `flowpaths` and `divides`,
#' builds event- and topology-aware exclusion sets, and derives auxiliary
#' tables (POI/downstream relationships) needed by the refactoring workflow.
#' Returns a ready-to-use list for downstream calls.
#' @param gpkg Optional path to a GeoPackage containing the layers: flowpaths, divides, events
#' @param flowpaths `sf` LINESTRINGs of flowpaths. Columns may use common
#'   NHDPlus-style names; they are normalized to:
#'   `id` (`COMID`), `toid` (`toCOMID`), `levelpathi` (`mainstemlp`),
#'   `lengthkm`, `areasqkm`, `hydroseq`, `reachcode`, `frommeas`, `tomeas`,
#'   `totdasqkm`, `terminalpa`. Geometry is renamed to `"geometry"`.
#' @param divides `sf` POLYGON/MULTIPOLYGON catchment divides. Columns may use
#'   `FEATUREID` which will be normalized to `divide_id`. Geometry is renamed
#'   to `"geometry"`.
#' @param events Optional `sf`/`data.frame` of points of interest (POIs)
#'   with at least `flowpath_id` (matching the original IDs in `flowpaths`)
#'   and `event_id`. Used to prevent collapsing immediately downstream of POIs
#'   and to assemble `outlets`.
#' @param avoid Optional integer/numeric vector of `COMID`s to exclude from
#'   collapse/aggregation (will be intersected with available flowpath IDs).
#'
#' @return A list with elements:
#' \describe{
#'   \item{rfc_flowpaths}{`sf` flowpaths normalized and re-labeled to the
#'         refactoring schema (e.g., `COMID`, `toCOMID`, `LENGTHKM`, etc.).}
#'   \item{events}{The input `events` (or `NULL`), unchanged.}
#'   \item{ex}{Integer vector of IDs to exclude from collapse, including
#'         user-provided `avoid`, highly elongated/large catchments, and
#'         flowlines immediately downstream of POIs.}
#'   \item{outlets}{Joined POI/outlet metadata for eligible flowpaths
#'         (or `NULL` if `events` is `NULL`).}
#'   \item{divides}{Input divides as `sf` with normalized names.}
#' }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Renames geometry to `"geometry"` for both layers.
#'   \item Lower-cases field names and applies lookup renames to harmonize
#'         common NHDPlus attributes.
#'   \item Derives an exclusion set (`ex`) including:
#'         (a) user-specified `avoid`,
#'         (b) flowpaths with extreme shape ratios `sqrt(areasqkm)/lengthkm > 3`
#'             and `areasqkm > 1`,
#'         (c) flowpaths immediately downstream of POIs.
#'   \item Derives `outlets` by joining `events` to flowpath metadata.
#'   \item Casts flowpaths to EPSG:5070 and removes Z/M dimensions.
#' }
#'
#' @keywords internal
#' @noRd
#' @importFrom hfutils rename_geometry
#' @importFrom dplyr rename select any_of filter mutate inner_join
#' @importFrom dplyr group_by ungroup distinct
#' @importFrom sf st_as_sf st_drop_geometry st_zm read_sf

.setup = function(gpkg, flowpaths, divides, events, avoid = NULL) {

  cat('\nRUNNING .SETUP()', file = stderr())
  
  if(!is.null(gpkg)){
    divides = sf::read_sf(gpkg, "divides")
    flowpaths = sf::read_sf(gpkg, "flowpaths")
    events = if(hfutils::layer_exists(gpkg, "events")){
      sf::read_sf(gpkg, "events")
    } else if(!is.null(events)){
      events
    } else {
      NULL
    }
  }

  divides   <- hfutils::rename_geometry(divides,   "geometry")
  flowpaths <- hfutils::rename_geometry(flowpaths, "geometry")

  names(divides)    <- tolower(names(divides))
  names(flowpaths)  <- tolower(names(flowpaths))

  # tony: changing the column names so they work with the NHD-style labels
  # as well as the new label used in the v2.3 parquet files
  #     New format (flowpath_id, flowpath_toid) → gets renamed to id, toid 
  #     NHD-style  (comid, tocomid) → gets renamed to id, toid 
  # I'm adding new columns so that we don't break the flowpath_id
  if (all(c("flowpath_id", "flowpath_toid") %in% names(flowpaths))) {
    flowpaths$id   <- flowpaths$flowpath_id
    flowpaths$toid <- flowpaths$flowpath_toid
  } else if (all(c("comid", "tocomid") %in% names(flowpaths))) {
    flowpaths$id   <- flowpaths$comid
    flowpaths$toid <- flowpaths$tocomid
  }
  if ("mainstemlp" %in% names(flowpaths)) {
    flowpaths$levelpathi <- flowpaths$mainstemlp
  }
  ## levelpathi -> mainstemlp doesn't need to be preserved, so rename is fine here
  #flowpaths <- dplyr::rename(flowpaths, any_of(c(levelpathi = "mainstemlp")))
  
  
  #fl_lookup  <- c(id = "flowpath_id", toid = "flowpath_toid", 
  #                id = "comid", toid = "tocomid", levelpathi = "mainstemlp")
  #flowpaths <- dplyr::rename(as.data.frame(flowpaths), any_of(fl_lookup))
  
  div_lookup <- c(FEATUREID = "divide_id", toid = "tocomid", levelpathi = "mainstemlp")
  divides   <- dplyr::rename(as.data.frame(divides), any_of(div_lookup)) |>
    st_as_sf()

  avoid_int <- filter(flowpaths, (sqrt(areasqkm) / lengthkm) > 3 & areasqkm > 1)
  avoid     <- c(avoid_int$id, avoid)
  avoid     <- avoid[avoid %in% flowpaths$id]

  if (!is.null(events)) {

    #match_id <- ifelse("flowpath_id" %in% names(flowpaths), "flowpath_id", "comid")

    match_id <- if ("flowpath_id" %in% names(flowpaths)) {
      "flowpath_id"
    } else if ("comid" %in% names(flowpaths)) {
      "comid"
    } else {
      "id"
    }
    
    outlets <-
      dplyr::inner_join(
        dplyr::mutate(events, flowpath_id = as.character(flowpath_id)),
        dplyr::select(sf::st_drop_geometry(flowpaths),
                      totdasqkm, dplyr::all_of(match_id),
                      dnhydroseq) |>
          dplyr::mutate(dplyr::across(dplyr::all_of(match_id), as.character)),
        by = c("flowpath_id" = match_id)
      )

    events <- dplyr::mutate(events, flowpath_id = as.character(flowpath_id))

    # Need to avoid modification to flowlines immediately downstream of POIs
    #      This can cause some hydrologically-incorrect catchment aggregation
    POI_downstream <- dplyr::filter(
      flowpaths,
      hydroseq %in% outlets$dnhydroseq,
      areasqkm > 0
    )

    ex <- unique(c(outlets$flowpath_id, avoid, POI_downstream$id))
  } else {
    events <- NULL
    outlets <- NULL
    ex <- unique(avoid)
  }

  # derive list of unique terminal paths
  TerminalPaths <- unique(flowpaths$terminalpa)

  flowpaths <- flowpaths |>
    dplyr::mutate(refactor = as.integer(terminalpa %in% unique(flowpaths$terminalpa)))  |>
    dplyr::rename(
      COMID = id, # was flowpath_id
      toCOMID = toid, # was flowpath_toid
      LENGTHKM = lengthkm,
      REACHCODE = reachcode,
      FromMeas = frommeas,
      ToMeas = tomeas,
      Hydroseq = hydroseq,
      LevelPathI = levelpathi,
      TotDASqkm = totdasqkm) |>
    st_as_sf(crs = 5070) |>
    sf::st_zm()

  return(list(rfc_flowpaths = flowpaths,
              events = events,
              ex = ex,
              outlets = outlets,
              divides = divides,
              vpuid = divides$vpuid[1]))
}


#' Refactor flowlines and reconcile catchment divides (end-to-end)
#'
#' @description
#' Wrapper that:
#' \enumerate{
#'   \item Normalizes inputs via an internal setup,
#'   \item Calls `refactor_flowpaths()` to split/collapse flowlines,
#'   \item Builds a lookup (`reconciled_ID` ↔ `member_COMID`),
#'   \item Writes refactored flowlines and reconciliation products to `outfile`,
#'   \item Reconciles and cleans divides against the refactored network using
#'         FDR/FAC (if provided),
#'   \item Emits auxiliary POI/outlet tables for QA and downstream use.
#' }
#'
#' @param flowpaths `sf` LINESTRINGs; reference flowline features.
#' @param divides `sf` POLYGON/MULTIPOLYGON; reference catchment divides.
#' @param events Optional POI table (`sf` or `data.frame`) with at least
#'   `event_id` and `flowpath_id` (pre-refactor). Used to steer splits and
#'   avoid collapsing immediately downstream of POIs.
#' @param avoid Optional integer/numeric vector of flowpath IDs (COMIDs)
#'   to exclude from collapse operations.
#' @param split_flines_meters Numeric. Maximum desired segment length in the
#'   refactored network (default `10000`).
#' @param collapse_flines_meters Numeric. Minimum inter-confluence length to
#'   retain during collapse (default `1000`).
#' @param collapse_flines_main_meters Numeric. Minimum between-confluence
#'   length for mainstem collapse (default `1000`).
#' @param min_area_m2 Numeric. Minimum polygon area in m² to keep when
#'   splitting divides (default `800`).
#' @param snap_distance_m Numeric. Snapping distance (meters) for aligning
#'   split polygons to outlets (default `100`).
#' @param simplify_tolerance_m Numeric. Douglas–Peucker tolerance (meters)
#'   used during divide splitting (default `40`).
#' @param fac Path to a flow accumulation raster. If supplied together with
#'   `fdr`, divide reconciliation is performed.
#' @param fdr Path to a flow direction raster. If supplied together with
#'   `fac`, divide reconciliation is performed.
#' @param keep Proportion of points to retain during **final** geometry
#'   simplification of divides (passed to `ms_simplify`). Use `NULL` to
#'   skip final simplification (default `NULL`).
#' @param outfile Path to a GeoPackage to write outputs. Layers written:
#'   `"refactored_flowpaths"`, `"reconciled_flowpaths"`, `"lookup_table"`,
#'   `"outlets"`, `"pois"`, and `"refactored_divides"`.
#'
#' @return (Invisibly) the path to `outfile` after writing the layers; also
#'   produces temporary files for intermediate refactor/reconcile steps.
#'
#' @section Outputs:
#' \describe{
#'   \item{refactored_flowpaths}{Refactored flowline network in EPSG:5070.}
#'   \item{reconciled_flowpaths}{Per-flowline reconciliation geometries (EPSG:5070).}
#'   \item{lookup_table}{Mapping of `reconciled_ID` ↔ `member_COMID` and base COMIDs.}
#'   \item{outlets}{Event-tied outlets on the refactored network.}
#'   \item{pois}{Original POI/outlet joins from setup.}
#'   \item{refactored_divides}{Unioned, reconciled divides labeled by `ID`
#'         with concatenated `member_COMID`; a `split_error` flag indicates
#'         fallback geometries where splitting failed.}
#' }
#'
#' @details
#' This wrapper writes outputs as it progresses. It constructs a reconciliation
#' table from the intermediate `"reconciled_flowpaths"` product to ensure
#' divides are split/unioned consistently with the refactored network.
#' If unioning via `hydrofab::union_polygons()` fails due to memory constraints,
#' a safe `sf::st_union()` fallback is used.
#'
#' @examples
#' \dontrun{
#' refactor(
#'   flowpaths = fl,
#'   divides   = dv,
#'   events    = pois,
#'   fac = "/path/fac.tif",
#'   fdr = "/path/fdr.tif",
#'   outfile = "refactor_out.gpkg"
#' )
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{reconcile_divides}} — divide reconciliation engine.
#'   \item \code{\link{clean_divide_geometry}} — post-reconciliation cleaner.
#'   \item \code{refactor_flowpaths} — network refactoring (package internal/external).
#' }
#'
#' @importFrom hfutils rename_geometry
#' @importFrom sf st_as_sf st_transform st_drop_geometry read_sf write_sf
#' @importFrom sf st_make_valid st_union
#' @importFrom dplyr select mutate rename filter inner_join left_join distinct
#' @importFrom dplyr group_by ungroup bind_rows slice_max as_tibble summarise n
#' @importFrom tidyr unnest separate_rows fill
#' @importFrom logger log_info
#' @importFrom tools file_path_sans_ext
#' @export

refactor <- function(gpkg      = NULL,
                     flowpaths = NULL,
                     divides   = NULL,
                     events    = NULL,
                     avoid     = NULL,
                     split_flines_meters = 10000,
                     collapse_flines_meters = 1000,
                     collapse_flines_main_meters = 1000,
                     min_area_m2 = 800,
                     snap_distance_m = 100,
                     simplify_tolerance_m = 40,
                     term_add = 1e9L,
                     modification_file = NULL,
                     fac  = NULL,
                     fdr  = NULL,
                     keep = NULL,
                     outfile = NULL) {

  if(is.null(outfile)){
    stop("Please provide an output geopackage path via the 'outfile' argument.")
  }

  tf <- tempfile(pattern = "refactored", fileext = ".gpkg")
  on.exit(unlink(tf), add = TRUE)

  tr <- tempfile(pattern = "reconciled", fileext = ".gpkg")
  on.exit(unlink(tf), add = TRUE)

  outfile_internal <- tempfile(pattern = "internal", fileext = ".gpkg")
  on.exit(unlink(outfile_internal), add = TRUE)

  ll <- .setup(gpkg, flowpaths, divides, events, avoid = avoid)

  suppressWarnings({
    refactor_flowpaths(
      flowpaths  = ll$rfc_flowpaths,
      max_length = split_flines_meters,
      collapse_meters = collapse_flines_meters,
      collapse_main_meters = collapse_flines_main_meters,
      out_refactored = tf,
      out_reconciled = tr,
      events = ll$events,
      exclude_ids = ll$ex#,
      #append = FALSE
    )
  })

  rec <- hfutils::rename_geometry(st_transform(read_sf(tr), 5070), "geometry")

  logger::log_info("creating lookup table")
  lookup_table <-
    sf::st_drop_geometry(rec) |>
    dplyr::select(ID, member_COMID) |>
    dplyr::mutate(member_COMID = strsplit(member_COMID, ",")) |>
    tidyr::unnest(cols = member_COMID) |>
    dplyr::mutate(NHDPlusV2_COMID = trunc(as.numeric(member_COMID))) |>
    dplyr::rename(reconciled_ID = ID) |>
    dplyr::mutate(vpuid = ll$vpuid) |>
    dplyr::distinct()

  logger::log_info("loading refactored flowpaths")
  refactored <-
    paste0(
      "SELECT
        COMID AS member_COMID,
        Hydroseq,
        event_identifier,
        NULLIF(event_REACHCODE, '') AS event_REACHCODE
      FROM ",
      basename(tools::file_path_sans_ext(tf))
    ) |>
    sf::read_sf(tf, query = _) |>
    dplyr::inner_join(
      dplyr::select(
        sf::st_drop_geometry(ll$rfc_flowpaths),
        orig_COMID = COMID,
        Hydroseq
      ),
      by = "Hydroseq"
    )

  # Refactored flowlines containing a split event
  event_outlets <-
    refactored |>
    sf::st_drop_geometry() |>
    # Keep rows where its COMID group has 1 row with an event identifier
    dplyr::group_by(orig_COMID) |>
    dplyr::filter(any(!is.na(event_identifier))) |>
    dplyr::mutate(event_identifier = as.integer(event_identifier)) |>
    # Ensure each row in the group has the event identifier and reachcode
    tidyr::fill(event_identifier, event_REACHCODE, .direction = "downup") |>
    dplyr::ungroup()

   if(nrow(event_outlets) > 0 & !is.null(ll$events)){
     event_outlets = event_outlets |>
        dplyr::inner_join(ll$events, by = c("event_identifier" = "identifier")) |>
        dplyr::select(flowpath_id, poi_id, member_COMID, event_identifier, geom)
   } else {
     event_outlets = NULL
   }

  # rm(events)

  # Refactored flowlines with POI that did not split
  if(!is.null(ll$outlets)){
    nonevent_outlets <-
      refactored |>
      sf::st_drop_geometry() |>
      dplyr::anti_join(event_outlets, by = "member_COMID") |>
      dplyr::select(member_COMID, orig_COMID) |>
      dplyr::mutate(orig_COMID = as.character(orig_COMID)) |>
      dplyr::inner_join(ll$outlets, by = c("orig_COMID" = "flowpath_id")) |>
      dplyr::slice_max(order_by = as.numeric(member_COMID), n = 1, by = orig_COMID, with_ties = FALSE) |>
      dplyr::select(flowpath_id = orig_COMID, poi_id, member_COMID, geom)
  } else {
    nonevent_outlets = NULL
  }


  final_outlets <- dplyr::bind_rows(event_outlets, nonevent_outlets)

  if(nrow(final_outlets) == 0){
    final_outlets <- NULL
  } else {
    final_outlets <- final_outlets |>
      sf::st_as_sf(crs = 5070) |>
      dplyr::left_join(
        dplyr::select(lookup_table, member_COMID, reconciled_ID),
        by = "member_COMID"
      )

    check_dups_poi <-
      final_outlets |>
      dplyr::group_by(reconciled_ID) |>
      dplyr::filter(dplyr::n() > 1) |>
      dplyr::ungroup()

    if (nrow(check_dups_poi) > 1) {
      warning("double-check for double POIs (check_dups_poi, refactor)")
    }
  }

  logger::log_info("writing refactored flowpaths")
  sf::read_sf(tf) |>
    sf::st_transform(5070) |>
    dplyr::mutate(vpuid = ll$vpuid[1]) |>
    sf::write_sf(outfile_internal, "refactored_flowpaths")

  rec |>
    sf::st_transform(5070) |>
    dplyr::rename(flowpath_id = ID, flowpath_toid = toID) |>
    dplyr::mutate(vpuid = ll$vpuid[1]) |>
    sf::write_sf(outfile_internal, "reconciled_flowpaths")

  logger::log_info("writing poi/outlets tables")
  sf::write_sf(lookup_table, outfile_internal, "lookup_table")
  if(!is.null(final_outlets)){
    sf::write_sf(final_outlets, outfile_internal, "outlets")
  }

  if(!is.null(ll$outlets)){
    sf::write_sf(ll$outlets, outfile_internal, "pois")
  }

  rm(lookup_table, final_outlets)

  # facfdr
  logger::log_info("starting divide reconcilation")

  rec2 <-
    rec |>
    sf::st_drop_geometry() |>
    tidyr::separate_rows(member_COMID, sep = ",") |>
    dplyr::filter(trunc(as.numeric(member_COMID)) %in% ll$divides$FEATUREID) |>
    dplyr::group_by(ID) |>
    dplyr::summarise(member_COMID = paste(member_COMID, collapse = ",")) |>
    dplyr::ungroup() |>
    dplyr::left_join(
      dplyr::select(rec, ID, geometry),
      by = "ID"
    ) |>
    sf::st_as_sf(crs = 5070)

  ref <-
    sf::read_sf(tf) |>
    dplyr::filter(trunc(as.numeric(COMID)) %in% ll$divides$FEATUREID)

  logger::log_info("reconciling divides")

  divides <-
   reconcile_divides(
      ll$divides,
      refactored_fp = ref,
      reconciled_fp = rec2,
      fdr = fdr,
      fac = fac,
      min_area_m = min_area_m2,
      snap_distance_m = snap_distance_m,
      simplify_tolerance_m = simplify_tolerance_m,
      keep = keep
    ) |>
    clean_divide_geometry(div,keep = NULL)

  divides <-
    divides |>
    dplyr::filter(is.na(ID)) |>
    dplyr::select(-ID) |>
    dplyr::left_join(
      sf::st_drop_geometry(rec) |>
        dplyr::select(ID, member_COMID) |>
        tidyr::separate_rows(member_COMID, sep = ","),
      by = "member_COMID"
    ) |>
    dplyr::bind_rows(dplyr::filter(divides, !is.na(ID)))

  # Handle divides that errored
  divides_errors <-
    dplyr::filter(divides, is.na(ID)) |>
    dplyr::mutate(split_error = TRUE) |>
    dplyr::summarise(
      member_COMID = paste(member_COMID, collapse = ","),
      geometry = dplyr::first(geometry),
      .by = member_COMID
    )

  rm(rec, ref, rec2)

  logger::log_info("writing refactored divides")
  # Join overlapping divides into a seamless divide
  tryCatch(
    hydrofab::union_polygons(divides, ID = "ID"),
    error = function(e) {
      # If the above fails (with a std::bad_alloc from terra's C++ code)
      # then try to do it with just sf + dplyr.
      divides |>
        sf::st_make_valid() |>
        hfutils::rename_geometry("geometry") |>
        dplyr::summarise(
          geometry = sf::st_union(geometry),
          .by = ID
        )
    }
  ) |>
    # Bring back member_COMID from before the union
    dplyr::left_join(
      sf::st_drop_geometry(divides) |>
        dplyr::select(ID, member_COMID),
      by = "ID"
    ) |>
    # Concat the member_COMIDs all together
    dplyr::as_tibble() |>
    dplyr::summarise(
      member_COMID = paste(member_COMID, collapse = ","),
      geometry = dplyr::first(geometry),
      .by = ID
    ) |>
    dplyr::mutate(split_error = FALSE) |>
    dplyr::bind_rows(divides_errors) |>
    sf::st_as_sf(crs = 5070) |>
    dplyr::mutate(vpuid = ll$vpuid[1]) |>
    # Finalize/output
    sf::write_sf(outfile_internal, "refactored_divides")

  logger::log_info("finished refactor")

  .clean(outfile_internal,
         gpkg = gpkg,
         term_add = term_add,
         modification_file = modification_file,
         out_gpkg = outfile)

  return(invisible(outfile))

}


#' Build reconciled HydroFabric outputs from a reference GeoPackage
#'
#' INTERNAL: Reads reconciliation/lookups from a reference GeoPackage, applies
#' optional cross-VPU modifications, computes a topological hydrosequence, fixes
#' geometries, and writes the final layers (flowpaths, divides, network,
#' hydrolocations, pois, events) to `out_gpkg`.
#'
#' @param outfile Path to the **reference** GeoPackage to read from. Must contain
#'   at least the tables/layers: `lookup_table`, `reconciled_flowpaths`,
#'   `refactored_divides`, `pois` (and any others referenced by the SQL).
#' @param term_add Integer suffix used to create synthetic terminal `toid`s when
#'   computing hydrosequence. Defaults to `1e9L`.
#' @param modification_file Optional path to a CSV of cross-VPU modifications.
#'   When provided, records are interpreted as `(from COMID) -> (to COMID)` and
#'   used to rewrite upstream `new_flowpath_toid` targets across VPUs before
#'   hydrosequence is finalized.
#'
#' @return Invisibly returns nothing; called for its side-effects:
#'   \itemize{
#'     \item Writes `flowpaths`, `divides`, `network`, `hydrolocations`, `pois`,
#'       and `events` layers to `out_gpkg`.
#'     \item Emits progress messages via \pkg{cli}.
#'   }
#'
#' @section Required objects (captured from calling environment):
#' \describe{
#'   \item{`out_gpkg`}{Destination GeoPackage path to write final layers.}
#'   \item{`ref_gpkg`}{GeoPackage providing reference geometries (e.g. `divides`).}
#'   \item{`rfc_gpkg`}{GeoPackage providing RFC `outlets` for `events`.}
#' }
#'
#' @details
#' Workflow (high level):
#' \enumerate{
#'   \item Build an identity map and new sequential `new_flowpath_id` per VPU.
#'   \item Optionally apply cross-VPU rewiring from `modification_file`.
#'   \item Assign synthetic terminal `toid`s (`term_add`) and compute a
#'         hydrosequence via \code{.topo_sort()} on the \code{(id,toid)} topology.
#'   \item Read, merge, and clean spatial layers:
#'         \code{reconciled_flowpaths} → flowpaths, \code{refactored_divides} →
#'         divides (handle collections, empty geometries, biggest polygon wins).
#'   \item Build `network` (joins hydrolocations + divide metrics), then derive
#'         final `flowpath`/`divide` tables and write all outputs to `out_gpkg`.
#' }
#'
#' Notes:
#' \itemize{
#'   \item Assumes SRID EPSG:5070 for spatial outputs.
#'   \item Uses `use_stream = TRUE` to stream SQL from `outfile` with \pkg{sf}.
#'   \item Expects helper \code{.topo_sort()} to be available in the package.
#' }
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom sf read_sf write_sf st_as_sf st_drop_geometry st_geometry_type
#'   st_is_empty st_collection_extract st_cast st_area st_difference st_combine
#'   st_geometry
#' @importFrom data.table as.data.table rbindlist
#' @importFrom dplyr select distinct right_join arrange mutate filter rename
#'   everything relocate group_by ungroup slice_max row_number left_join
#' @importFrom tidyr separate_longer_delim
#' @importFrom dtplyr lazy_dt
#' @importFrom cli cli_progress_step cli_progress_done cli_alert_info
#' @importFrom rlang .data
#' @importFrom utils head

.clean = function(outfile, gpkg = NULL, term_add = 1e9L, modification_file = NULL, out_gpkg) {

  lookup = read_sf(outfile, "lookup_table")

  identity_map <-
    sf::read_sf(outfile, query = "
    WITH
      crosswalk AS (
        SELECT
          a.reconciled_ID AS flowpath_id,
          b.flowpath_toid,
          a.member_COMID,
          a.NHDPlusV2_COMID,
          a.vpuid
        FROM
          lookup_table a
        INNER JOIN
          (SELECT flowpath_id, flowpath_toid, vpuid FROM reconciled_flowpaths) b
          ON a.reconciled_ID = b.flowpath_id AND a.vpuid = b.vpuid
      ),
      new_identifiers AS (
        SELECT DISTINCT
          ROW_NUMBER() OVER (ORDER BY vpuid, flowpath_id) AS new_flowpath_id,
          flowpath_id,
          vpuid
        FROM
          (SELECT flowpath_id, vpuid FROM reconciled_flowpaths)
      )
    SELECT
      b.new_flowpath_id,
      a.flowpath_id,
      a.flowpath_toid,
      a.member_COMID,
      a.NHDPlusV2_COMID,
      a.vpuid
    FROM
      crosswalk a
    LEFT JOIN
      new_identifiers b ON a.vpuid = b.vpuid AND a.flowpath_id = b.flowpath_id
  ", use_stream = TRUE) |>
    data.table::as.data.table()

  cli::cli_progress_step("building modified identity map")
  identity_map <-
    identity_map |>
    dplyr::select(
      flowpath_toid = flowpath_id,
      new_flowpath_toid = new_flowpath_id,
      vpuid
    ) |>
    dplyr::distinct() |>
    dplyr::right_join(identity_map, by = c("vpuid", "flowpath_toid")) |>
    dplyr::select(
      flowpath_id,
      flowpath_toid,
      new_flowpath_id,
      new_flowpath_toid,
      vpuid,
      member_COMID,
      NHDPlusV2_COMID
    ) |>
    dplyr::arrange(vpuid, flowpath_id)

  # Crosswalk between new flowpath IDs and NHDPlusV2 COMIDs
  cli::cli_progress_step("creating crosswalk between new flowpath IDs and NHDPlusV2 COMIDs")
  new_identifiers_crosswalk <-
    identity_map |>
    dplyr::select(new_flowpath_id, member_COMID = NHDPlusV2_COMID) |>
    dplyr::distinct()

  if(is.null(modification_file)){
    cli::cli_alert_info("No modification file provide, skipping modifications")
    mods = data.frame(to = NA, from = NA)
  } else {
    cli::cli_progress_step("applying modifications to crosswalk")
    mods <-
      read.csv(modification_file) |>
      dplyr::filter(VPUID != toVPUID) |>
      dplyr::rename(from = COMID, to = toCOMID) |>
      dplyr::select(from, to) |>
      dplyr::mutate(connection = dplyr::row_number()) |>
      dplyr::left_join(new_identifiers_crosswalk, by = c("from" = "member_COMID")) |>
      dplyr::rename(new_from = new_flowpath_id) |>
      dplyr::left_join(new_identifiers_crosswalk, by = c("to" = "member_COMID")) |>
      dplyr::rename(new_to = new_flowpath_id) |>
      filter(!is.na(new_to))

    if(nrow(mods) == 0){
      cli::cli_alert_info("No modifications to apply")
      mods = data.frame(to = NA, from = NA)
    } else {
      identity_map$new_flowpath_toid[
        match(mods$new_from, identity_map$new_flowpath_id)
      ] <- mods$new_to
    }
  }

  cli::cli_progress_step("finalizing crosswalk terminals")
  to_map <- which(is.na(identity_map$new_flowpath_toid))

  identity_map$new_flowpath_toid[to_map] <- seq(1L, length(to_map)) + term_add

  cli::cli_progress_step("computing hydrosequence")
  topology <-
    identity_map |>
    dtplyr::lazy_dt() |>
    dplyr::select(-NHDPlusV2_COMID) |>
    dplyr::group_by(new_flowpath_id) |>
    dplyr::mutate(
      member_COMID = paste0(unique(member_COMID), collapse=",")
    ) |>
    dplyr::ungroup() |>
    dplyr::distinct() |>
    data.table::as.data.table()

  topology <-
    topology |>
    dplyr::select(id = new_flowpath_id, toid = new_flowpath_toid) |>
    .topo_sort() |>
    merge(x = topology, by.x = c("new_flowpath_id", "new_flowpath_toid"), by.y = c("id", "toid"), all.x = TRUE) |>
    dplyr::rename(
      reconciled_id = flowpath_id,
      reconciled_toid = flowpath_toid,
      flowpath_id = new_flowpath_id,
      flowpath_toid = new_flowpath_toid
    )

  cli::cli_progress_step("processing flowpaths")
  flowpaths <-
    sf::read_sf(outfile, query = "
    SELECT
      flowpath_id AS reconciled_id,
      vpuid,
      LevelPathID AS levelpathid,
      member_COMID,
      geom AS geometry
    FROM
      reconciled_flowpaths
  ", use_stream = TRUE) |>
    data.table::as.data.table() |>
    merge(
      topology[, .(reconciled_id, vpuid, flowpath_id, flowpath_toid, hydroseq)],
      by = c("vpuid", "reconciled_id")
    ) |>
    dplyr::select(
      flowpath_id,
      flowpath_toid,
      reconciled_id,
      vpuid,
      levelpathid,
      hydroseq,
      member_COMID,
      geometry
    ) |>
    dplyr::mutate(
      lengthkm = hydrofab::add_lengthkm(geometry),
      ibt = flowpath_id %in% mods$from
    ) |>
    sf::st_as_sf(crs = 5070)

  cli::cli_progress_step("processing divides")
  divides <-
    sf::read_sf(outfile, query = "
    SELECT
      ID AS divide_id,
      vpuid,
      geom AS geometry
    FROM
      refactored_divides
    WHERE
      ID IS NOT NULL
  ", use_stream = TRUE) |>
    data.table::as.data.table() |>
    merge(topology[, .(divide_id = reconciled_id, vpuid, flowpath_id, member_COMID)], by = c("vpuid", "divide_id")) |>
    dplyr::select(reconciled_id = divide_id, divide_id = flowpath_id, vpuid, member_COMID, geometry)

  divides_mask <-
    which(sf::st_geometry_type(divides$geometry) == "GEOMETRYCOLLECTION")

  if (length(divides_mask) > 0) {
    divides <-
      divides[divides_mask, ] |>
      sf::st_as_sf(crs = 5070) |>
      sf::st_collection_extract("POLYGON", warn = FALSE) |>
      hydrofab::union_polygons("divide_id") |>
      data.table::as.data.table() |>
      dplyr::mutate(geometry = sf::st_cast(geometry, "MULTIPOLYGON")) |>
      merge(divides[divides_mask, .(divide_id, vpuid, member_COMID)], by = "divide_id") |>
      list(divides[!divides_mask, .(divide_id, geometry, vpuid, member_COMID)]) |>
      data.table::rbindlist(ignore.attr = TRUE) |>
      sf::st_as_sf(crs = 5070)
  }

  divides_mask <-
    which(sf::st_is_empty(divides$geometry))

  if (length(divides_mask) > 0) {
    comids <- trunc(as.numeric(divides[divides_mask, ]$member_COMID))

    ref_divides <-
      sf::read_sf(
        ref_gpkg,
        query = paste0(
          "SELECT divide_id AS NHDPlusV2_COMID, geom AS geometry FROM divides WHERE divide_id IN (",
          paste0(comids, collapse = ","),
          ")"
        )
      )

    divides2 <-
      divides |>
      dplyr::filter(trunc(as.numeric(member_COMID)) %in% comids) |>
      suppressWarnings() |>
      dplyr::mutate(
        NHDPlusV2_COMID = trunc(as.numeric(member_COMID)),
        empty = sf::st_is_empty(geometry)
      )

    for (comid in comids) {
      divides2$geometry[divides2$NHDPlusV2_COMID == comid & divides2$empty] <-
        sf::st_difference(
          sf::st_cast(sf::st_geometry(ref_divides)[ref_divides$NHDPlusV2_COMID == comid], "MULTIPOLYGON"),
          sf::st_combine(sf::st_cast(divides2[divides2$NHDPlusV2_COMID == comid & !divides2$empty, geometry], "MULTIPOLYGON"))
        ) |>
        sf::st_as_sf() |>
        sf::st_cast("POLYGON") |>
        dplyr::mutate(area = sf::st_area(x)) |>
        dplyr::slice_max(order_by = area, n = 1, with_ties = FALSE, na_rm = TRUE) |>
        sf::st_geometry() |>
        sf::st_cast("MULTIPOLYGON")
    }

    for (mask in divides_mask) {
      divides$geometry[mask] <- divides2[divides2$empty & match(divides2$member_COMID, divides$member_COMID[mask]), ]$geometry
    }

    rm(divides2)
  }

  rm(divides_mask)

  divides <-
    divides |>
    data.table::as.data.table() |>
    dplyr::distinct() |>
    sf::st_as_sf(crs = 5070)

  divides_mask <- which(vapply(sf::st_geometry(divides), length, integer(1)) > 1)
  if (length(divides_mask) > 0) {
    divides$geometry[divides_mask] <-
      divides[divides_mask, ] |>
      sf::st_cast("POLYGON", warn = FALSE) |>
      dplyr::mutate(area = sf::st_area(geometry)) |>
      dplyr::slice_max(order_by = area, n = 1, by = divide_id, with_ties = FALSE, na_rm = TRUE) |>
      dplyr::select(-area) |>
      sf::st_as_sf(crs = 5070) |>
      sf::st_cast("MULTIPOLYGON") |>
      sf::st_geometry()
  }
  rm(divides_mask)

  divides_clean <-
    hydrofab::clean_geometry(
      catchments = divides,
      flowlines = flowpaths,
      fl_ID = "flowpath_id",
      ID = "divide_id",
      keep = NULL,
      crs = 5070,
      grid = .0009,
      gb = 8,
      force = FALSE,
      sys = NULL
    ) |>
    dplyr::mutate(areasqkm = hydrofab::add_areasqkm(geometry))

  cli::cli_progress_step("processing lookup table")
  lookup <-
    topology |>
    dplyr::select(vpuid, flowpath_id, reconciled_id, member_COMID) |>
    tidyr::separate_longer_delim(member_COMID, delim = ",") |>
    data.table::as.data.table() |>
    dplyr::mutate(NHDPlusV2_COMID = trunc(as.numeric(member_COMID)))


  if(!is.null(gpkg) & 'hydrolocations' %in% st_layers(gpkg)$name){

      cli::cli_progress_step("processing hydrolocations")
      hydrolocations <-
        sf::read_sf(gpkg, query = "
        SELECT
          flowpath_id AS NHDPlusV2_COMID,
          poi_id,
          hl_id,
          hl_link,
          hl_type,
          hl_reference,
          hl_class,
          hl_source,
          hl_uri,
          pathmeas,
          reachmeas,
          ST_X(geom) AS X,
          ST_Y(geom) AS Y
        FROM
          hydrolocations
      ", use_stream = TRUE) |>
        data.table::as.data.table() |>
        merge(lookup, by = "NHDPlusV2_COMID")

      cli::cli_progress_step("processing network")
      divide_network <-
        sf::st_drop_geometry(divides_clean) |>
        data.table::as.data.table() |>
        dplyr::select(divide_id, areasqkm) |>
        merge(
          data.table::as.data.table(sf::st_drop_geometry(flowpaths)),
          by.x = "divide_id",
          by.y = "flowpath_id",
          all = TRUE
        ) |>
        dplyr::mutate(flowpath_id = divide_id)

      network <-
        hydrolocations |>
        dplyr::select(
          flowpath_id,
          poi_id,
          hl_id,
          hl_link,
          hl_type,
          hl_reference,
          hl_class,
          hl_source,
          hl_uri,
          pathmeas,
          reachmeas
        ) |>
        merge(divide_network, by = "flowpath_id", all = TRUE) |>
        dplyr::mutate(topo = "fl-fl") |>
        tidyr::separate_longer_delim(member_COMID, ",") |>
        data.table::as.data.table() |>
        dplyr::mutate(hf_id = trunc(as.numeric(member_COMID))) |>
        dplyr::select(-member_COMID) |>
        dplyr::relocate(flowpath_id, flowpath_toid, divide_id) |>
        dplyr::arrange(hydroseq)
  } else {
    cli::cli_progress_step("skipping hydrolocations/network (no hydrolocations provided or in gpkg)")
    
    # Create a minimal network from divides for the flowpath processing
    # if the hydrolocations/network step is being skipped.
    # This ensures that the 'finalizing flowpaths' step doesn't fail.
    network <-
      sf::st_drop_geometry(divides_clean) |>
      data.table::as.data.table() |>
      dplyr::select(divide_id, areasqkm) |>
      merge(
        data.table::as.data.table(sf::st_drop_geometry(flowpaths)),
        by.x = "divide_id",
        by.y = "flowpath_id",
        all = TRUE
      ) |>
      dplyr::mutate(flowpath_id = divide_id) |>
      dplyr::relocate(flowpath_id, flowpath_toid, divide_id) |>
      dplyr::arrange(hydroseq)
  }

  cli::cli_progress_step("finalizing flowpaths")
  flowpath <-
    network |>
    dplyr::select(divide_id, flowpath_id, areasqkm) |>
    dplyr::distinct() |>
    dplyr::filter(
      dplyr::everything() |>
        dplyr::pick() |>
        complete.cases()
    ) |>
    merge(x = data.table::as.data.table(flowpaths), y = _, by = "flowpath_id", all.x = TRUE) |>
    dplyr::arrange(hydroseq)

  suppressWarnings({
    flowpath[, totdasqkm := cumsum(areasqkm), by = "levelpathid"]
    flowpath[, has_divide := !is.na(divide_id)]
    colnames(flowpath) <- tolower(colnames(flowpath))
  })


  cli::cli_progress_step("finalizing divides")
  divide <-
    network |>
    dplyr::select(divide_id, flowpath_id) |>
    dplyr::distinct() |>
    dplyr::filter(
      dplyr::everything() |>
        dplyr::pick() |>
        complete.cases()
    ) |>
    merge(x = data.table::as.data.table(divides_clean), y = _, by = "divide_id", all.x = TRUE)

  divide[, has_flowpath := !is.na(flowpath_id)]
  divide$member_COMID <- NULL
  colnames(divide) <- tolower(colnames(divide))

  mask <- sf::st_geometry_type(divide$geometry) != "POLYGON"

  if (any(mask)) {
    tmp <-
      divide[which(mask), .(divide_id, geometry = sf::st_as_sfc(lapply(geometry, identity)))] |>
      sf::st_as_sf(crs = 5070) |>
      sf::st_cast("POLYGON", warn = FALSE) |>
      rmapshaper::ms_dissolve(field = "divide_id", snap = TRUE) |>
      rmapshaper::ms_explode() |>
      dplyr::slice_max(order_by = sf::st_area(geometry), n = 1, by = divide_id, with_ties = FALSE)

    divide$geometry[which(mask)] <- sf::st_geometry(tmp)

    rm(tmp)
  }

  cli::cli_progress_step("writing flowpaths")
  flowpath |>
    sf::st_as_sf(crs = 5070) |>
    sf::write_sf(out_gpkg, "flowpaths")

  cli::cli_progress_step("writing divides")
  divide |>
    sf::st_as_sf(crs = 5070) |>
    sf::write_sf(out_gpkg, "divides")

  cli::cli_progress_step("writing network")
  network |>
    sf::write_sf(out_gpkg, "network")

  if(!is.null(gpkg) & 'hydrolocations' %in% st_layers(gpkg)$name){
    cli::cli_progress_step("writing hydrolocations")
    hydrolocations |>
      dplyr::select(-member_COMID) |>
      dplyr::rename(hf_id = NHDPlusV2_COMID) |>
      sf::st_as_sf(coords = c("X", "Y"), crs = 5070, remove = FALSE) |>
      sf::write_sf(out_gpkg, "hydrolocations")
  } else {
    cli::cli_progress_step("skipping writing hydrolocations, none provided")
  }

  if(!is.null(gpkg) & 'pois' %in% st_layers(gpkg)$name){
    cli::cli_progress_step("writing pois")
    read_sf(gpkg, "pois") |>
      select(-vpuid) |>
      rename(NHDPlusV2_COMID = flowpath_id) |>
      inner_join(lookup, by = "NHDPlusV2_COMID") |>
      sf::write_sf(out_gpkg, "pois")
  } else {
    cli::cli_progress_step("skipping writing pois, none provided")
  }

#   cli::cli_progress_step("writing events")
#   sf::read_sf(outfile, query = "
#   SELECT DISTINCT
#     event_identifier,
#     flowpath_id AS NHDPlusV2_COMID,
#     poi_id,
#     vpuid,
#     geom AS geometry
#   FROM
#     outlets
# ", use_stream = TRUE, promote_to_multi = FALSE) |>
#     data.table::as.data.table() |>
#     merge(lookup, by = c("vpuid", "NHDPlusV2_COMID")) |>
#     sf::write_sf(out_gpkg, "events")

  cli::cli_progress_done()

}


#' Topologically sort a hydrological network
#'
#' This internal helper performs a depth-first traversal of a transposed
#' network to generate a consistent hydrosequence ordering. The function
#' is designed for hydrological fabric data, where traversing the
#' transposed network corresponds to moving upstream.
#'
#' @param topology A `data.table` with at least two integer columns:
#'   \describe{
#'     \item{id}{Unique node identifier.}
#'     \item{toid}{Downstream node identifier.}
#'     \item{term_add}{(optional) integer threshold used to add terminal
#'     connections to a synthetic root node.}
#'   }
#'
#' @return A `data.table` with columns:
#'   \describe{
#'     \item{id}{Node identifier.}
#'     \item{toid}{Downstream node identifier.}
#'     \item{hydroseq}{Hydrological sequence index, giving a topological
#'     ordering of the network.}
#'   }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Builds a transposed edge list where edges point upstream.
#'   \item Adds a synthetic root (id = 0) connecting all terminal nodes.
#'   \item Runs a depth-first search (DFS) from the synthetic root to
#'   assign an initial order.
#'   \item Resolves ties and reassigns sequential hydrosequence indices.
#' }
#'
#' This ensures a deterministic upstream-to-downstream ordering suitable
#' for hydrofabric workflows.
#'
#' @keywords internal

.topo_sort <- function(topology, term_add = 1e9L) {
  # Create a transpose network, where traversing the network
  # is equivalent to traversing the hydrological network upstream.
  edgelist <-
    data.table::rbindlist(list(
      topology[, .(id = toid, toid = id)],
      topology[toid >= term_add, .(id = 0L, toid)]
    ))

  # Perform DFS from each terminal upstream to get a
  # distinct topological sort for the hydrosequence.
  sorted <- data.table::data.table(
    node = as.integer(
      names(
        igraph::dfs(
          igraph::graph_from_data_frame(edgelist),
          root = "0",
          mode = "out"
        )$order
      )[-1]
    )
  )[, .(node, hydroseq = .I)]

  # Merge the initial hydrosequence to the eddgelist, removing the fake "0"
  # connections for the terminals, and handling ties in the hydrosequence.
  merge(edgelist, sorted, by.x = "id", by.y = "node", all.x = TRUE)[
    !is.na(hydroseq)
  ][
    order(hydroseq, id, toid)
  ][
    , .(id, toid, hydroseq = .I)
  ][
    , .(id = toid, toid = id, hydroseq)
  ]
}









