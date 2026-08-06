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
#' @param required_parameters Optional character vector of
#'   `parameter_name`s (see [harmonize_columns()]/`column_synonyms.csv`);
#'   if given, each result is additionally fetched and harmonized to check
#'   whether it actually has at least one of these parameters (OR
#'   semantics), and results with none of them are dropped -- e.g. "only
#'   show search results that have TOC data". `NULL` (default) skips this
#'   entirely, keeping the default search exactly as fast/cheap as a plain
#'   text search (`pg_search()` alone doesn't know about parameter
#'   content, so checking it means fetching every result's actual data --
#'   real, additive cost, unlike [find_neighbours()]'s equivalent filter,
#'   which reuses fetches it has to make anyway). When set, an
#'   `available_parameters` column always records what was found (`NA`
#'   when not checked, i.e. `required_parameters` is `NULL`).
#' @param progress Optional `function(i, n, doi)`, called after each of
#'   the `n` results is checked when `required_parameters` is set (see
#'   [find_neighbours()]'s equivalent argument) -- otherwise never called,
#'   since nothing is fetched.
#' @return A tibble with one row per (core, search hit): `core_id`, `score`,
#'   `doi`, `size`, `size_measure`, `citation`, `available_parameters`.
#' @export
pg_search_cores <- function(x, count = 20, required_parameters = NULL, progress = NULL) {
  validate_core_registry(x)

  purrr::map_dfr(x$core_id, function(core_id) {
    res <- tryCatch(
      pangaear::pg_search(query = core_id, count = count),
      error = function(e) NULL
    )
    if (is.null(res) || nrow(res) == 0) {
      return(tibble::tibble(
        core_id = core_id, score = NA_real_, doi = NA_character_,
        size = NA_real_, size_measure = NA_character_, citation = NA_character_,
        available_parameters = NA_character_
      ))
    }
    out <- tibble::as_tibble(res)[, c("score", "doi", "size", "size_measure", "citation")] |>
      dplyr::mutate(core_id = core_id, .before = 1)

    if (is.null(required_parameters)) {
      out$available_parameters <- NA_character_
      return(out)
    }

    # Only reached when the caller actually asked for this -- pg_search()
    # itself has no idea what parameters a result contains, so finding
    # out means fetching each one's real data, which plain pg_search_cores()
    # never had to do before.
    n <- nrow(out)
    out$available_parameters <- purrr::imap_chr(out$doi, function(doi, i) {
      result <- if (is.na(doi)) {
        NA_character_
      } else {
        fetched <- tryCatch(pg_fetch_dataset(doi, overwrite = FALSE), error = function(e) NULL)
        if (is.null(fetched) || nrow(fetched) == 0) {
          NA_character_
        } else {
          params <- unlist(lapply(seq_len(nrow(fetched)), function(j) {
            h <- tryCatch(
              harmonize_columns(fetched$data[[j]], core_id = core_id, doi = fetched$doi[j], unmatched = "silent"),
              error = function(e) NULL
            )
            if (is.null(h) || nrow(h) == 0) character() else unique(h$parameter_name)
          }))
          if (length(params) == 0) NA_character_ else paste(sort(unique(params)), collapse = ", ")
        }
      }
      if (!is.null(progress)) {
        progress(i, n, doi)
      }
      result
    })

    has_required <- vapply(out$available_parameters, function(ap) {
      if (is.na(ap)) {
        return(FALSE)
      }
      any(trimws(strsplit(ap, ",")[[1]]) %in% required_parameters)
    }, logical(1))
    out[has_required, ]
  })
}
