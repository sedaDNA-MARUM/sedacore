#' Classify a core's drilling program from its ID or label
#'
#' Regex router used to populate `core_registry$source_program` and to
#' decide whether an IODP campaign search applies to a given core (see
#' Phase 5, PANGAEA IODP campaign search).
#' Recognizes MARUM/AWI-style `GeoB` core names, IODP hole IDs (`U1385A`,
#' optionally prefixed with an expedition number, e.g. `"361-U1479A"`),
#' legacy ODP leg-site-hole IDs (e.g. `"108-659A"`), and DSDP IDs. Anything
#' else is classified as `"other"`.
#'
#' @param core_id A character vector of core IDs or labels.
#' @return A character vector the same length as `core_id`, one of
#'   `"GeoB"`, `"IODP"`, `"ODP"`, `"DSDP"`, `"other"`.
#' @export
classify_core_source <- function(core_id) {
  core_id <- as.character(core_id)

  dplyr::case_when(
    is.na(core_id) ~ NA_character_,
    stringr::str_detect(core_id, stringr::regex("^GeoB", ignore_case = TRUE)) ~ "GeoB",
    stringr::str_detect(core_id, stringr::regex("IODP", ignore_case = TRUE)) ~ "IODP",
    stringr::str_detect(core_id, "^[0-9]{2,3}-U[0-9]{4}[A-Z]?$|^U[0-9]{4}[A-Z]?$") ~ "IODP",
    stringr::str_detect(core_id, stringr::regex("DSDP", ignore_case = TRUE)) ~ "DSDP",
    stringr::str_detect(core_id, stringr::regex("ODP", ignore_case = TRUE)) ~ "ODP",
    stringr::str_detect(core_id, "^[0-9]{2,3}-[0-9]{3,4}[A-Z]?$") ~ "ODP",
    TRUE ~ "other"
  )
}

#' Search PANGAEA's IODP/ODP/DSDP campaign index
#'
#' PANGAEA hosts an IODP data portal (<https://iodp.pangaea.de/>) and
#' indexes ODP/DSDP/IODP expedition data through the *same* search
#' infrastructure as regular datasets, filterable with a `campaign:`
#' query prefix -- confirmed empirically, e.g.
#' `pg_search(query = "campaign:361")` returns real IODP Expedition 361
#' (Site U1478) datasets. This means IODP support needs no separate
#' parser: results feed straight into the same [pg_fetch_dataset()] /
#' [harmonize_columns()] pipeline used for GeoB cores.
#'
#' @param expedition An expedition/campaign code (e.g. `"361"`), or
#'   `NULL` to search by `bbox` alone.
#' @param bbox Optional bounding box (see [bbox_from_radius()]) to
#'   restrict the search geographically.
#' @param count Max results to request.
#' @return A tibble: `score`, `doi`, `size`, `size_measure`, `citation`.
#' @export
search_iodp_campaigns <- function(expedition = NULL, bbox = NULL, count = 20) {
  if (is.null(expedition) && is.null(bbox)) {
    rlang::abort("Provide at least one of `expedition` or `bbox`.")
  }
  query <- if (!is.null(expedition)) paste0("campaign:", expedition) else "*"

  res <- tryCatch(
    pangaear::pg_search(query = query, bbox = bbox, count = count),
    error = function(e) NULL
  )
  if (is.null(res) || nrow(res) == 0) {
    return(tibble::tibble(
      score = numeric(), doi = character(), size = numeric(),
      size_measure = character(), citation = character()
    ))
  }
  tibble::as_tibble(res)[, c("score", "doi", "size", "size_measure", "citation")]
}
