# Shared plumbing for the Shiny app's modules.
#
# The app is a set of modules (R/mod-*.R) over one shared `state`
# (a reactiveValues) plus two reactives passed in by the modules that own
# them (`active_core`, `measurement_filter`). Keeping all of this in the
# package (rather than in inst/app/) is what makes each module testable
# with shiny::testServer() from the ordinary test suite.

#' Create the shared app state
#'
#' @return A `shiny::reactiveValues` with the registries, the accumulated
#'   `param_long`, the derived `depth_matched_long`, transient search
#'   results, and the log.
#' @keywords internal
#' @export
new_app_state <- function() {
  shiny::reactiveValues(
    core_registry = NULL,
    sample_registry = NULL,
    param_long = empty_param_long(),
    depth_matched_long = empty_depth_matched_long(),
    search_results = empty_search_results(),
    zenodo_search_results = empty_zenodo_search_results(),
    neighbour_results = empty_neighbour_registry(),
    log = character()
  )
}

#' Prepend a timestamped line to the app log
#'
#' @param state A [new_app_state()].
#' @param msg The message.
#' @keywords internal
#' @export
log_message <- function(state, msg) {
  state$log <- c(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg), state$log)
  invisible(state)
}

#' Every parameter name any synonym table can produce
#'
#' Drawn from *every* synonym table (PANGAEA and Zenodo/IODP), so
#' source-specific names (e.g. Zenodo's `total_carbon`) are selectable in
#' the measurement filter from the start.
#'
#' @return A sorted character vector.
#' @keywords internal
#' @export
known_parameter_names <- function() {
  sort(unique(c(load_synonym_map()$parameter_name, load_zenodo_synonym_map()$parameter_name)))
}

#' Fetch-and-harmonize services bound to the app state
#'
#' Returns the two retrieval routes (PANGAEA, Zenodo) as functions that
#' fetch a record, harmonize every table in it with the right source
#' conventions, and append the result to `state$param_long`, logging as
#' they go. `measurement_filter` is read *when a fetch runs*, so changing
#' the filter between fetches takes effect without rebuilding anything.
#'
#' @param state A [new_app_state()].
#' @param measurement_filter A reactive (or function) returning the
#'   currently selected `parameter_name`s.
#' @param fetch_pangaea,fetch_zenodo The network calls, injectable so the
#'   modules can be tested offline.
#' @return A list with `pangaea(doi, core_id)` and `zenodo(doi, core_id)`.
#' @keywords internal
#' @export
make_fetchers <- function(state, measurement_filter,
                          fetch_pangaea = pg_fetch_dataset,
                          fetch_zenodo = zenodo_fetch_dataset) {
  store <- function(raw_data, core_id, doi, depth_col, synonym_map, synonym_csv_name) {
    res <- harmonize_for_store(
      raw_data, core_id, doi, depth_col, synonym_map,
      measurement_filter = measurement_filter(), synonym_csv_name = synonym_csv_name
    )
    if (!is.null(res$data)) {
      state$param_long <- dplyr::bind_rows(state$param_long, res$data)
    }
    for (msg in res$messages) log_message(state, msg)
  }

  list(
    pangaea = function(doi, core_id) {
      log_message(state, sprintf("Fetching %s for %s...", doi, core_id))
      fetched <- tryCatch(fetch_pangaea(doi), error = function(e) {
        log_message(state, paste("ERROR fetching:", conditionMessage(e)))
        NULL
      })
      if (is.null(fetched)) {
        return(invisible())
      }
      for (i in seq_len(nrow(fetched))) {
        store(fetched$data[[i]], core_id, fetched$doi[i],
          depth_col = "Depth sed [m]", synonym_map = load_synonym_map(),
          synonym_csv_name = "column_synonyms.csv"
        )
      }
    },

    # Zenodo/IODP records can bundle several files per DOI (e.g.
    # data_by_hole.zip and data_by_splice.zip on one record); the fetcher
    # already expands those to one row per file, same shape as PANGAEA's,
    # so this loop matches the one above -- pointed at Zenodo's own depth
    # column (CSF-A scale, not PANGAEA's "Depth sed") and synonym table.
    zenodo = function(doi, core_id) {
      log_message(state, sprintf("Fetching Zenodo record %s for %s...", doi, core_id))
      fetched <- tryCatch(fetch_zenodo(doi), error = function(e) {
        log_message(state, paste("ERROR fetching from Zenodo:", conditionMessage(e)))
        NULL
      })
      if (is.null(fetched) || nrow(fetched) == 0) {
        log_message(state, sprintf("No tabular data file found in Zenodo record %s.", doi))
        return(invisible())
      }
      for (i in seq_len(nrow(fetched))) {
        store(fetched$data[[i]], core_id, fetched$doi[i],
          depth_col = "Top depth CSF-A (m)", synonym_map = load_zenodo_synonym_map(),
          synonym_csv_name = "zenodo_column_synonyms.csv"
        )
      }
    }
  )
}
