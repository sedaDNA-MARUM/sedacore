#' Assemble source adapters' output into one bundle's tables
#'
#' Each adapter ([build_palmod()], and the atlas adapter to come) returns
#' tables keyed by its own site names. This stacks them, folds the sites onto
#' canonical ones ([build_crosswalk()]), re-keys every table, drops links that
#' became duplicates, and marks one preferred series per measurement
#' ([mark_preferred_series()]), so that coverage counts each record once.
#'
#' @param sources A list of adapter results, each with `tables` and a
#'   `source` description (as returned by [build_palmod()] and
#'   [read_wpa_sources()]). A source needs only `sites` and whichever other
#'   tables it has.
#' @param priority Source ids, highest priority first.
#' @param max_km,overrides Passed to [build_crosswalk()]; overrides default
#'   to the shipped `site_crosswalk_overrides.csv`.
#' @return A list: `tables` (ready for [write_bundle()] with `crosswalked =
#'   TRUE`), `sources` (manifest entries) and `review` (what a person should
#'   look at; see [build_crosswalk()]).
#' @keywords internal
#' @export
assemble_bundle <- function(sources, priority, max_km = 25, overrides = load_crosswalk_overrides()) {
  ids <- vapply(sources, function(s) s$source$id, character(1))
  if (anyDuplicated(ids)) rlang::abort("Source ids must be unique.")
  unknown <- setdiff(ids, priority)
  if (length(unknown) > 0) {
    rlang::abort(sprintf("Source(s) missing from `priority`: %s", paste(unknown, collapse = ", ")))
  }

  stacked <- dplyr::bind_rows(lapply(sources, function(s) {
    dplyr::mutate(s$tables$sites, source = s$source$id)
  }))
  cw <- build_crosswalk(stacked, priority, max_km, overrides)

  rekeyed <- lapply(sources, function(s) apply_crosswalk(s$tables, cw$site_aliases, s$source$id))
  bind <- function(nm) {
    out <- dplyr::bind_rows(lapply(rekeyed, `[[`, nm))
    # no source supplied this table: a typed, empty one rather than no columns
    if (ncol(out) == 0 && !is.null(bundle_schema[[nm]])) out <- empty_from_spec(as.list(bundle_schema[[nm]]))
    out
  }

  catalog <- mark_preferred_series(bind("series_catalog"), priority)
  series <- add_preferred_to_series(bind("series"), catalog)

  links <- dplyr::distinct(bind("site_links"), .data$site_id, .data$link_type, .data$url, .keep_all = TRUE)

  sites <- cw$sites
  sites$n_aliases <- NULL
  tables <- list(
    sites = sites, site_aliases = cw$site_aliases, series = series, series_catalog = catalog,
    age_control = bind("age_control"), site_links = links
  )
  records <- bind("indexed_records")
  if (nrow(records) > 0) tables$indexed_records <- records
  for (nm in names(tables)) check_bundle_table(tables[[nm]], nm)

  # nothing may point at a site that does not exist (an index record with no
  # site yet is allowed: its site_id is NA)
  for (nm in intersect(c("series", "series_catalog", "age_control", "site_links", "indexed_records"), names(tables))) {
    orphans <- setdiff(stats::na.omit(tables[[nm]]$site_id), sites$site_id)
    if (length(orphans) > 0) {
      rlang::abort(sprintf("Table `%s` refers to unknown site(s): %s", nm, paste(utils::head(orphans, 5), collapse = ", ")))
    }
  }
  list(tables = tables, sources = lapply(sources, `[[`, "source"), review = cw$review)
}
