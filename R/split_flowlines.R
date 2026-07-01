#' Split flowlines (sequential, flowpath_id-based)
#' @description Split flowlines by events and/or a max segment length.
#' Returns semantic IDs with .part suffix ordered upstream (.1) -> downstream (.N).
#'
#' @param flines sf LINESTRINGs with columns:
#'   - flowpath_id (coerced to character)
#'   - one downstream column: flowpath_toid | toCOMID | toid  (kept as-is)
#'   - LENGTHKM (will be recomputed from geometry)
#'   CRS should be projected in meters for length-based splitting.
#' @param max_length numeric (meters). If NULL, only event splits apply.
#' @param events data.frame of event matches with **per-event measures**:
#'   required: flowpath_id, REACHCODE, FromMeas, ToMeas, REACH_meas
#'   optional: identifier (passed through to output)
#' @param avoid vector of flowpath_id to exclude from *length* splitting (events still apply).
#' @return sf with split flowlines; includes event_* fields and recomputed LENGTHKM.
#' @export
split_flowlines <- function(flines,
                            max_length = NULL,
                            events = NULL,
                            avoid = NA) {

  stopifnot(inherits(flines, "sf"))
  flines <- .coerce_and_check_fp(flines)  # ensures flowpath_id, finds downstream col
  dn_col <- attr(flines, "downstream_col")

  # Geometry sanity
  if (is.na(sf::st_crs(flines))) {
    warning("Input has no CRS; lengths assumed to be meters.")
  } else if (sf::st_is_longlat(flines)) {
    warning("CRS appears geographic (lon/lat). Length-based splitting expects meters.")
  }

  # Compute split windows (fractions in [0,1])
  split_pts <- .split_points_fp(flines,
                                max_length = max_length,
                                events = events,
                                avoid = avoid)

  # Nothing to split: normalize and return
  if (nrow(split_pts) == 0) {
    out <- dplyr::mutate(flines,
                         flowpath_id = as.character(flowpath_id),
                         !!dn_col := as.character(.data[[dn_col]]),
                         event_REACHCODE = NA_character_,
                         event_REACH_meas = NA_real_,
                         event_identifier = NA_character_)
    return(out)
  }

  # Build geometries sequentially
  built <- .build_split_lines_fp(flines, split_pts)

  # Bring back network attrs and compute parts/LENGTHKM
  base_attrs <- flines |>
    sf::st_set_geometry(NULL) |>
    dplyr::select(flowpath_id, dplyr::all_of(dn_col),
                  dplyr::any_of(c("mainstemlp","hydroseq","totdasqkm")))

  split <- dplyr::left_join(built, base_attrs, by = "flowpath_id") |>
    dplyr::group_by(flowpath_id) |>
    dplyr::arrange(end, .by_group = TRUE) |>
    dplyr::mutate(part = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::group_by(flowpath_id) |>
    dplyr::mutate(
      !!dn_col := ifelse(part == max(part),
                         as.character(.data[[dn_col]]),
                         paste0(flowpath_id, ".", part + 1))
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(flowpath_id = paste0(flowpath_id, ".", part)) |>
    dplyr::select(-part)

  split <- dplyr::mutate(split, lengthkm = .length_km(sf::st_geometry(split)))

  # Originals that weren’t split
  base_ids <- unique(.base_id(split$flowpath_id))
  not_split <- flines |>
    dplyr::filter(!flowpath_id %in% base_ids) |>
    dplyr::mutate(
      event_REACHCODE = NA_character_,
      event_REACH_meas = NA_real_,
      event_identifier = NA_character_,
      flowpath_id = as.character(flowpath_id),
      !!dn_col := as.character(.data[[dn_col]])
    )

  # Merge and redirect downstream pointers that target a newly split base id -> ".1"
  out <- dplyr::bind_rows(not_split, split) |>
    dplyr::mutate(
      !!dn_col := ifelse(.data[[dn_col]] %in% base_ids,
                         paste0(.data[[dn_col]], ".1"),
                         .data[[dn_col]])
    )

  out
}

# ---- helpers ---------------------------------------------------------------

.coerce_and_check_fp <- function(x) {
  dn_candidates <- c("flowpath_toid", "toCOMID", "toid")
  dn_col <- dn_candidates[dn_candidates %in% names(x)][1]
  if (is.na(dn_col)) stop("Downstream ID column not found. Provide one of: ", paste(dn_candidates, collapse = ", "))

  id_candidates = c("COMID", "ID", "flowpath_id")
  id_col <- id_candidates[id_candidates %in% names(x)][1]
  if (is.na(id_col)) stop("ID column not found. Provide one of: ", paste(id_candidates, collapse = ", "))

  ln_candidates = c("LENGTHKM", "lengthkm")
  ln_col <- ln_candidates[ln_candidates %in% names(x)][1]
  if (is.na(ln_col)) stop("Length column not found. Provide one of: ", paste(ln_candidates, collapse = ", "))

  req <- c(id_col, dn_col, ln_col)
  miss <- setdiff(req, names(x))
  if (length(miss)) stop("Missing required columns in `flines`: ", paste(miss, collapse = ", "))

  x$flowpath_id <- as.character(x[[id_col]])
  x$flowpath_toid   <- as.character(x[[dn_col]])
  attr(x, "downstream_col") <- dn_col
  x
}

# Build unified split points from events + length rules ---------------------------------

.split_points_fp <- function(flines, max_length, events, avoid) {
  ev <- .split_by_event_fp(flines, events)
  lb <- .split_by_length_fp(flines, max_length, ev, avoid)

  if (nrow(ev) > 0 && nrow(lb) > 0) {
    # When events exist, lb rows carry event_split_fID where applicable; adjust event-local windows
    lb <- lb |>
      dplyr::left_join(
        dplyr::select(ev, e_start = start, e_end = end, event_split_fID = event_split_fID),
        by = "event_split_fID"
      ) |>
      dplyr::mutate(
        start = ifelse(!is.na(e_start), start * (e_end - e_start) + e_start, start),
        end   = ifelse(!is.na(e_start), end   * (e_end - e_start) + e_start, end)
      ) |>
      dplyr::mutate(
        event_REACHCODE = ifelse(
          !is.na(event_REACHCODE) && event_REACHCODE != "" &&
            round(end, 5) == round(e_end, 5),
          event_REACHCODE, NA_character_
        ),
        event_REACH_meas = ifelse(is.na(event_REACHCODE), NA_real_, event_REACH_meas),
        event_identifier = ifelse(is.na(event_identifier), NA_character_, event_identifier)
      ) |>
      dplyr::select(-e_start, -e_end)
  } else if (nrow(ev) == 0 && nrow(lb) > 0 && !"event_identifier" %in% names(lb)) {
    lb$event_identifier <- seq_len(nrow(lb))
  }
  lb
}

# Event windows (now sourced entirely from `events`, not `flines`) ----------------------

.split_by_event_fp <- function(flines, events) {
  if (is.null(events) || NROW(events) == 0) {
    return(data.frame(
      flowpath_id = character(0), start = numeric(0), end = numeric(0),
      event_REACHCODE = character(0), event_REACH_meas = numeric(0),
      event_identifier = character(0), event_split_fID = integer(0)
    ))
  }

  # required on events
  req_ev <- c("flowpath_id", "REACHCODE", "FromMeas", "ToMeas", "REACH_meas")
  miss <- setdiff(req_ev, names(events))
  if (length(miss)) stop("Missing required columns in `events`: ", paste(miss, collapse = ", "))

  if (!"identifier" %in% names(events)) {
    events$identifier <- as.character(seq_len(nrow(events)))
  } else {
    events$identifier <- as.character(events$identifier)
  }

  # only events for flowpaths present in flines
  keep_ids <- unique(flines$flowpath_id)
  ev_all <- dplyr::filter(events, flowpath_id %in% keep_ids)

  if (nrow(ev_all) == 0) {
    return(data.frame(
      flowpath_id = character(0), start = numeric(0), end = numeric(0),
      event_REACHCODE = character(0), event_REACH_meas = numeric(0),
      event_identifier = character(0), event_split_fID = integer(0)
    ))
  }

  # Build per-flowpath windows from 0..1 using event measures and the event's own From/To
  out_list <- lapply(split(ev_all, ev_all$flowpath_id), function(df) {
    # upstream (larger measure) to downstream (smaller)
    df <- dplyr::arrange(df, dplyr::desc(REACH_meas))
    # normalize: map REACH_meas from [FromMeas, ToMeas] -> [0, 100]
    # then invert so 0 is outlet, 100 inlet, to match substring behavior
    norm <- 100 * (df$REACH_meas - df$FromMeas) / pmax(1e-9, (df$ToMeas - df$FromMeas))

    data.frame(
      flowpath_id      = df$flowpath_id[1],
      start            = c(0, 100 - norm) / 100,
      end              = c(100 - norm, 100) / 100,
      event_REACHCODE  = c(df$REACHCODE, ""),
      event_REACH_meas = c(df$REACH_meas, NA_real_),
      event_identifier = c(df$identifier, "")
    )
  })

  out <- do.call(rbind, out_list)
  out$event_split_fID <- seq_len(nrow(out))
  out
}

# Length-based splitting (respects event windows when present) -------------------------

.split_by_length_fp <- function(flines, max_length, event_pts, avoid) {
  if (is.null(max_length)) {
    if (nrow(event_pts) == 0) return(data.frame())
    return(dplyr::mutate(event_pts, fID = event_split_fID))
  }

  tmp <- dplyr::select(flines, flowpath_id)
  len_tab <- dplyr::mutate(tmp, geom_len = .length_m(sf::st_geometry(tmp))) |>
    sf::st_set_geometry(NULL)
  rm(tmp)

  if (nrow(event_pts) > 0) {
    event_pts <- dplyr::left_join(event_pts,
                                  dplyr::select(len_tab, flowpath_id, geom_len),
                                  by = "flowpath_id") |>
      dplyr::mutate(geom_len = geom_len * (end - start))

    base_rows <- len_tab |>
      dplyr::filter(!flowpath_id %in% unique(event_pts$flowpath_id)) |>
      dplyr::mutate(
        start = 0, end = 1,
        event_REACHCODE = NA_character_,
        event_REACH_meas = NA_real_,
        event_identifier = NA_character_,
        event_split_fID = NA_integer_
      )

    lines <- dplyr::bind_rows(
      dplyr::select(event_pts, flowpath_id, start, end,
                    event_REACHCODE, event_REACH_meas,
                    event_identifier, event_split_fID, geom_len),
      base_rows
    )
  } else {
    lines <- len_tab |>
      dplyr::mutate(
        start = 0, end = 1,
        event_REACHCODE = NA_character_,
        event_REACH_meas = NA_real_,
        event_identifier = NA_character_,
        event_split_fID = NA_integer_
      )
  }

  if (max_length < 50) warning("Very small max_length detected (<50 m). Are your units meters?)")

  avoid <- as.character(avoid)
  need <- lines |>
    dplyr::filter((geom_len >= max_length & !flowpath_id %in% avoid) |
                    !is.na(event_split_fID)) |>
    dplyr::mutate(
      pieces = pmax(1L, ceiling(geom_len / max_length)),
      fID = seq_len(dplyr::n())
    ) |>
    dplyr::select(-geom_len)

  if (nrow(need) == 0) {
    if (nrow(event_pts) == 0) return(data.frame())
    return(dplyr::mutate(
      dplyr::select(event_pts, flowpath_id, start, end,
                    event_REACHCODE, event_REACH_meas, event_split_fID),
      fID = event_split_fID
    ))
  }

  need_rep <- need[rep(seq_len(nrow(need)), need$pieces), , drop = FALSE] |>
    dplyr::select(-pieces) |>
    dplyr::group_by(fID) |>
    dplyr::mutate(
      piece = dplyr::row_number(),
      nP    = dplyr::n(),
      start = (piece - 1) / nP,
      end   =  piece / nP
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-piece, -nP) |>
    dplyr::group_by(flowpath_id) |>
    dplyr::arrange(end, .by_group = TRUE) |>
    dplyr::mutate(
      split_fID = ifelse(dplyr::row_number() == 1L,
                         flowpath_id,
                         paste0(flowpath_id, ".", dplyr::row_number() - 1L))
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(fID, end)

  need_rep
}

# Geometry build ----------------------------------------------------------------------

.build_split_lines_fp <- function(flines, split_pts) {
  gtab <- flines |>
    dplyr::select(flowpath_id) |>
    sf::st_as_sf()

  idx <- match(split_pts$flowpath_id, gtab$flowpath_id)
  if (anyNA(idx)) stop("Internal mismatch: split points reference unknown flowpath_id(s).")

  geoms <- sf::st_geometry(gtab)[idx]

  new_geoms <- mapply(
    FUN = function(g, s, e) lwgeom::st_linesubstring(g, from = s, to = e),
    g   = geoms,
    s   = split_pts$start,
    e   = split_pts$end,
    SIMPLIFY = FALSE
  )

  sf::st_sf(
    split_pts[, c("flowpath_id", "start", "end",
                  "event_REACHCODE", "event_REACH_meas", "event_identifier",
                  "split_fID", "fID")],
    geometry = sf::st_sfc(new_geoms, crs = sf::st_crs(flines))
  )
}

# small utils ------------------------------------------------------------------------

.length_m  <- function(geom) as.numeric(sf::st_length(geom))
.length_km <- function(geom) (.length_m(geom) / 1000)
.base_id   <- function(x) sub("\\..*$", "", x)
