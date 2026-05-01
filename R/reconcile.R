#' @title Reconcile Collapsed Flowlines
#' @description Reconciles output of collapse_flowlines giving a unique ID to
#' each new flowpath and providing a mapping to NHDPlus COMIDs.
#' @param flines data.frame with COMID, toCOMID, LENGTHKM, LevelPathI, Hydroseq,
#' and TotDASqKM columns
#' @param geom sf data.frame for flines
#' @param id character id collumn name.
#' @return reconciled flowpaths with new ID, toID, LevelPathID, and Hydroseq identifiers.
#' Note that all the identifiers are new integer IDs. LevelPathID and Hydroseq are consistent
#' with the LevelPathID and Hydroseq from the input NHDPlus flowlines.
#' @importFrom dplyr group_by ungroup filter left_join select rename
#' mutate distinct summarise arrange desc
#' @seealso The \code{\link{refactor_nhdplus}} function implements a complete
#' workflow using `reconcile_collapsed_flowlines()`.
#' @export


reconcile_collapsed_flowlines <- function(flines, geom = NULL, id = "COMID") {

  new_flines <-
    mutate(flines,
           becomes =
             ifelse((is.na(joined_fromCOMID) | joined_fromCOMID == -9999),
                    ifelse((is.na(joined_toCOMID) | joined_toCOMID == -9999),
                           COMID, joined_toCOMID
                    ),
                    joined_fromCOMID
             )
    )

  # In the case that something is first joined to then the thing it joins with joins from
  # we have to do a little reassignment.
  joined_from <- new_flines[new_flines$joined_fromCOMID %in% new_flines$becomes, ]
  joined_to <- new_flines[new_flines$joined_toCOMID %in% new_flines$becomes, ]

  joined_tofrom <- joined_to[joined_to$becomes %in% joined_from$COMID, ]

  update_tofrom <- left_join(select(joined_tofrom, COMID, becomes),
                             select(joined_from, COMID, new_becomes = becomes),
                             by = c("becomes" = "COMID")
  )

  if (nrow(update_tofrom) > 0) {
    new_flines <- left_join(new_flines, select(update_tofrom, COMID, new_becomes), by = "COMID") %>%
      mutate(becomes = ifelse(!is.na(new_becomes), new_becomes, becomes)) %>%
      select(-new_becomes)
  }

  new_flines <- new_flines %>%
    group_by(becomes) %>%
    mutate(
      LENGTHKM = max(LENGTHKM),
      Hydroseq = min(Hydroseq),
      LevelPathI = min(LevelPathI)
    ) %>%
    ungroup() %>%
    tidyr::separate(COMID, c("orig_COMID", "part"),
                    sep = "\\.", remove = FALSE, fill = "right"
    ) %>%
    mutate(
      new_Hydroseq = ifelse(is.na(part),
                            as.character(Hydroseq),
                            paste(as.character(Hydroseq),
                                  part,
                                  sep = "."
                            )
      ),
      part = ifelse(is.na(part), "0", part)
    ) %>%
    ungroup() %>%
    select(-joined_fromCOMID, -joined_toCOMID)

  new_flines <-
    left_join(new_flines,
              data.frame(
                becomes = unique(new_flines$becomes),
                ID = seq_len(length(unique(new_flines$becomes))),
                stringsAsFactors = FALSE
              ),
              by = "becomes"
    )

  tocomid_updater <- filter(
    select(new_flines, becomes, toCOMID),
    !is.na(toCOMID)
  )

  new_flines <- distinct(left_join(select(new_flines, -toCOMID),
                                   tocomid_updater,
                                   by = "becomes"
  ))

  new_flines <- left_join(new_flines,
                          select(new_flines, becomes, toID = ID),
                          by = c("toCOMID" = "becomes")
  ) %>%
    arrange(Hydroseq, desc(part))

  new_flines <- left_join(new_flines,
                          data.frame(
                            ID = unique(new_flines$ID),
                            ID_Hydroseq = seq_len(length(unique(new_flines$ID)))
                          ),
                          by = "ID"
  )

  new_lp <- group_by(new_flines, LevelPathI) %>%
    filter(Hydroseq == min(Hydroseq)) %>% # Get the outlet by hydrosequence.
    ungroup() %>%
    group_by(Hydroseq) %>%
    # Get the outlet if the original was split.
    filter(as.numeric(part) == max(as.numeric(part))) %>%
    ungroup() %>%
    select(ID_LevelPathID = ID_Hydroseq, LevelPathI)

  if (!"event_identifier" %in% names(new_flines)) {
    new_flines$event_identifier <- rep(NA, nrow(new_flines))
  }

  new_flines <- left_join(distinct(new_flines), distinct(new_lp), by = "LevelPathI") %>%
    select(ID, toID, LENGTHKM,
           member_COMID = COMID,
           LevelPathID = ID_LevelPathID, Hydroseq = ID_Hydroseq,
           event_identifier, orig_levelpathID = LevelPathI
    )

  if (!is.null(geom)) {
    geom_column <- attr(geom, "sf_column")

    if (is.null(geom_column)) stop("geom must contain an sf geometry column")

    new_flines <- select(geom, member_COMID = id, geom_column) |>
      mutate(member_COMID = as.character(member_COMID)) |>
      right_join(new_flines,by = "member_COMID")

    new_flines <- new_flines %>%
      st_drop_geometry() %>%
      group_by(ID) %>%
      summarise(
        toID = toID[1],
        LENGTHKM = LENGTHKM[1],
        LevelPathID = LevelPathID[1],
        Hydroseq = Hydroseq[1],
        event_identifier = event_identifier[1],
        orig_levelpathID = orig_levelpathID[1],
        member_COMID = list(unique(member_COMID))
      ) %>%
      ungroup() %>%
      left_join(
        hfutils::union_linestrings(
          select(new_flines[!sf::st_is_empty(new_flines), ], ID), "ID"
        ),
        by = "ID"
      ) %>%
      sf::st_as_sf()
  }

  return(new_flines)

}
