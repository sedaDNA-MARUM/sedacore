#' Search PANGAEA for datasets matching each core
#'
#' Runs a [pangaear::pg_search()] text query per core (by `core_id`) and
#' returns all candidates for review -- PANGAEA search relevance varies
#' enough (see the sibling `sedaDNA_atlantic_LGM` project, where a single
#' core has isotope, TOC, susceptibility, and water-density data split
#' across separate datasets) that auto-picking the top hit per core is not
#' reliable. Feed the result to [pg_fetch_dataset()] once the right DOI(s)
#' are identified (interactively, or via the Shiny app's review step).
#'
#' @param x A `core_registry`, as produced by [collapse_to_core_registry()].
#' @param count Max results to request per core from PANGAEA.
#' @return A tibble with one row per (core, search hit): `core_id`, `score`,
#'   `doi`, `size`, `size_measure`, `citation`.
#' @export
pg_search_cores <- function(x, count = 20) {
  validate_core_registry(x)

  purrr::map_dfr(x$core_id, function(core_id) {
    res <- tryCatch(
      pangaear::pg_search(query = core_id, count = count),
      error = function(e) NULL
    )
    if (is.null(res) || nrow(res) == 0) {
      return(tibble::tibble(
        core_id = core_id, score = NA_real_, doi = NA_character_,
        size = NA_real_, size_measure = NA_character_, citation = NA_character_
      ))
    }
    tibble::as_tibble(res)[, c("score", "doi", "size", "size_measure", "citation")] |>
      dplyr::mutate(core_id = core_id, .before = 1)
  })
}
