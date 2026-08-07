#' Import contextual data from a manual export
#'
#' Documented fallback for data that isn't (yet) reachable through any of
#' the package's automated sources (PANGAEA, Zenodo's IODP community) --
#' e.g. a still-in-moratorium IODP expedition, a table copied out of a
#' published paper's supplementary material, or a one-off export from
#' some other database. IODP's LIMS Online Reports Explorer
#' (<https://web.iodp.tamu.edu/LORE/>) in particular has no documented
#' public API for automated retrieval, so this is the recommended path
#' for it specifically too: export the report you need, reshape it by
#' hand into the template at `system.file("extdata",
#' "external_import_template.csv", package = "sedacore")`
#' (`depth_cm, parameter_group, parameter_name, species, value, unit`),
#' and import it here into the standard [new_param_long()] schema so it
#' can be combined with PANGAEA/Zenodo-sourced data downstream, flowing
#' through the exact same harmonized/depth-matched/checklist/export
#' pipeline with no special-casing.
#'
#' @param path Path to a CSV shaped like `external_import_template.csv`.
#' @param core_id The core this import belongs to.
#' @param source Optional free-text label for where this came from (e.g.
#'   `"LORE export"`, `"Smith et al. 2020, Table 2"`) -- kept as a
#'   passthrough column for provenance, distinct from `doi` since a
#'   manual source often has no DOI at all.
#' @param doi Optional provenance DOI/identifier for the manual source,
#'   if it has one.
#' @return A `param_long` tibble, with an extra `source` column.
#' @export
import_external_data <- function(path, core_id, source = NA_character_, doi = NA_character_) {
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
    retrieved_at = Sys.time(),
    source = source
  )
  new_param_long(out)
}
