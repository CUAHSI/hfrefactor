refactor_flowpaths <- function(
    flowpaths,
    max_length,
    exclude_ids,
    events,
    collapse_meters,
    collapse_main_meters,
    out_refactored,
    out_reconciled
) {

  flowpaths <- flowpaths |>
    sf::st_cast("LINESTRING", warn = FALSE) |>
    sf::st_transform(5070) |>
    hydrofab::split_flowlines(
      max_length = 10000,
      events = events,
      avoid = exclude_ids
    )

  exclude_ids <- c(
    exclude_ids,
    dplyr::filter(flowpaths, !is.na(event_REACH_meas)) |>
      dplyr::select(COMID, toCOMID) |>
      unlist()
  )

  # three_pass
  collapsed_flowpaths <-
    collapse_flowlines(
      flines = sf::st_set_geometry(flowpaths, NULL),
      thresh = (0.1 * collapse_meters / 1000),
      add_category = TRUE,
      mainstem_thresh = (0.1 * collapse_main_meters / 1000),
      exclude_cats = exclude_ids) |>
    collapse_flowlines(
      thresh = 0.5 * collapse_meters / 1000,
      add_category = TRUE,
      mainstem_thresh = 0.5 * collapse_main_meters / 1000,
      exclude_cats = exclude_ids) |>
    collapse_flowlines(
      thresh = collapse_meters / 1000,
      add_category = TRUE,
      mainstem_thresh = collapse_main_meters / 1000,
      exclude_cats = exclude_ids
    )

    dplyr::select(flowpaths, COMID) |>
      mutate(COMID = as.character(COMID)) |>
      dplyr::inner_join(collapsed_flowpaths, by = "COMID") |>
      sf::st_as_sf() |>
      sf::st_transform(5070) |>
      sf::st_write(out_refactored, layer_options = "OVERWRITE=YES", append = FALSE)

  collapsed <-
    reconcile_collapsed_flowlines(
      flines = collapsed_flowpaths,
      geom = dplyr::select(flowpaths, COMID),
      id = "COMID"
    )

  collapsed[["member_COMID"]] <-
    collapsed[["member_COMID"]] |>
    lapply(paste, collapse = ",") |>
    unlist()

  collapsed |>
    sf::st_transform(5070) |>
    sf::st_write(out_reconciled, layer_options = "OVERWRITE=YES", append = FALSE)
}
