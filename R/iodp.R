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

#' Import IODP contextual data from a manual export
#'
#' Documented fallback for IODP/ODP/DSDP data not (yet) mirrored to
#' PANGAEA -- e.g. a still-in-moratorium expedition. IODP's LIMS Online
#' Reports Explorer (<https://web.iodp.tamu.edu/LORE/>) has no documented
#' public API for automated retrieval, so this is an occasional-use manual
#' path, not automation: export the report you need from LORE, reshape it
#' by hand into the template at
#' `system.file("extdata", "iodp_manual_template.csv", package =
#' "sedacore")` (`depth_cm, parameter_group, parameter_name, species,
#' value, unit`), and import it here into the standard [new_param_long()]
#' schema so it can be combined with PANGAEA-sourced data downstream.
#'
#' @param path Path to a CSV shaped like `iodp_manual_template.csv`.
#' @param core_id The core this import belongs to.
#' @param doi Optional provenance DOI/identifier for the manual source.
#' @return A `param_long` tibble.
#' @export
iodp_manual_import <- function(path, core_id, doi = NA_character_) {
  raw <- readr::read_csv(path, show_col_types = FALSE)

  required <- c("depth_cm", "parameter_group", "parameter_name", "value", "unit")
  missing_cols <- setdiff(required, names(raw))
  if (length(missing_cols) > 0) {
    rlang::abort(sprintf(
      "`path` is missing required column(s): %s", paste(missing_cols, collapse = ", ")
    ))
  }
  if (!"species" %in% names(raw)) {
    raw$species <- NA_character_
  }

  out <- tibble::tibble(
    core_id = core_id,
    pangaea_doi = doi,
    depth_cm = as.numeric(raw$depth_cm),
    parameter_group = as.character(raw$parameter_group),
    parameter_name = as.character(raw$parameter_name),
    species = as.character(raw$species),
    value = as.numeric(raw$value),
    unit = as.character(raw$unit),
    raw_column_name = NA_character_,
    is_interpolated = FALSE,
    interp_gap_cm = NA_real_,
    retrieved_at = Sys.time()
  )
  new_param_long(out)
}
