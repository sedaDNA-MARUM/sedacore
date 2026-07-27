#' Load the dataset-override table
#'
#' @param path Path to `dataset_overrides.R`.
#' @return A named list of `function(raw_data) -> raw_data`, keyed by DOI.
#' @noRd
get_dataset_overrides <- function(
    path = system.file("extdata", "dataset_overrides.R", package = "sedacore")) {
  env <- new.env(parent = baseenv())
  sys.source(path, envir = env)
  env$dataset_overrides
}

#' Apply a dataset-specific override, if one exists, before harmonization
#'
#' @param doi The dataset DOI.
#' @param raw_data The raw parameter table (as returned in
#'   `pg_fetch_dataset()`'s `data` list-column).
#' @param overrides The override table, from `get_dataset_overrides()`.
#' @return `raw_data`, transformed by the override function if `doi` has
#'   one registered, otherwise unchanged.
#' @export
apply_dataset_override <- function(doi, raw_data, overrides = get_dataset_overrides()) {
  fn <- overrides[[doi]]
  if (is.null(fn)) {
    return(raw_data)
  }
  fn(raw_data)
}
