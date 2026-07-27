#' Extract water depth (m, positive downward) from a PANGAEA elevation string
#'
#' PANGAEA event metadata reports `ELEVATION` as a negative-metres string
#' like `"-980.0 m"` (elevation relative to sea level, so seafloor depth is
#' negative). This returns the absolute value, matching `core_registry`'s
#' `water_depth_m` convention (positive downward).
#'
#' @param x A character scalar, e.g. `"-980.0 m"`, or `NA`.
#' @return A numeric scalar, or `NA_real_` if unparseable.
#' @noRd
parse_elevation_to_water_depth <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) {
    return(NA_real_)
  }
  val <- suppressWarnings(as.numeric(stringr::str_extract(x, "-?[0-9.]+")))
  if (is.na(val)) {
    return(NA_real_)
  }
  abs(val)
}

#' @noRd
event_field <- function(events, field) {
  if (is.null(events) || !field %in% names(events)) {
    return(NA_character_)
  }
  val <- events[[field]]
  if (is.null(val) || length(val) == 0) NA_character_ else as.character(val[[1]])
}

#' @noRd
known_event_fields <- c(
  "LATITUDE", "LONGITUDE", "DATE/TIME", "ELEVATION", "RECOVERY",
  "LOCATION", "CAMPAIGN", "BASIS", "METHOD/DEVICE", "COMMENT"
)

#' Extract the event/core label from PANGAEA event metadata
#'
#' PANGAEA's `events` metadata list carries the event/core label (e.g.
#' `"GeoB1701-4"`) as its own list *name*, with a `NA` value -- confirmed
#' empirically, it's always the one entry that isn't one of the known
#' fixed fields (`LATITUDE`, `CAMPAIGN`, etc.).
#'
#' @param events The `metadata$events` list from a `pg_data()` entry.
#' @return A character scalar, or `NA_character_` if not found.
#' @noRd
event_label <- function(events) {
  if (is.null(events)) {
    return(NA_character_)
  }
  candidates <- setdiff(names(events), known_event_fields)
  if (length(candidates) == 0) NA_character_ else candidates[[1]]
}

#' @noRd
one_pangaea_entry_to_row <- function(entry) {
  events <- entry$metadata$events
  tibble::tibble(
    doi           = entry$doi %||% NA_character_,
    citation      = entry$citation %||% NA_character_,
    event_label   = event_label(events),
    latitude      = suppressWarnings(as.numeric(event_field(events, "LATITUDE"))),
    longitude     = suppressWarnings(as.numeric(event_field(events, "LONGITUDE"))),
    water_depth_m = parse_elevation_to_water_depth(event_field(events, "ELEVATION")),
    campaign      = event_field(events, "CAMPAIGN"),
    date_time     = suppressWarnings(as.POSIXct(
      event_field(events, "DATE/TIME"),
      format = "%Y-%m-%dT%H:%M:%S", tz = "UTC"
    )),
    data          = list(entry$data)
  )
}

#' Fetch a PANGAEA dataset by DOI
#'
#' Thin wrapper around [pangaear::pg_data()] with retry and a result shape
#' that's uniform whether `doi` points at a single dataset or a
#' parent/collection DOI spanning multiple events (PANGAEA returns the
#' latter as a list of one-event-each child datasets, e.g. all cores from
#' one campaign) -- confirmed empirically, since single-event and
#' collection DOIs are not otherwise distinguishable up front. Also pulls
#' latitude/longitude/water depth/campaign out of each child's event
#' metadata, since (also confirmed empirically) these are *not* embedded as
#' columns in `$data` for single-event PANGAEA datasets -- they live in
#' `$metadata$events` and must be parsed out.
#'
#' @param doi A PANGAEA DOI (e.g. `"10.1594/PANGAEA.54904"`).
#' @param retries Number of retry attempts on failure.
#' @param overwrite Passed to [pangaear::pg_data()]. `TRUE` (default)
#'   always re-downloads, guaranteeing fresh data -- important for the
#'   harmonization pipeline, where stale cached data would silently go
#'   unnoticed. Neighbour search ([find_neighbours()]) passes `FALSE`
#'   instead, since it may fetch the same candidate DOI's coordinates
#'   repeatedly across searches and only needs the coordinates, not
#'   guaranteed-fresh measurements.
#' @return A tibble with one row per child dataset: `doi`, `citation`,
#'   `event_label`, `latitude`, `longitude`, `water_depth_m`, `campaign`,
#'   `date_time` (collection date, from event metadata; `NA` if absent),
#'   and `data` (a list-column holding each child's raw parameter table).
#' @export
pg_fetch_dataset <- function(doi, retries = 2, overwrite = TRUE) {
  check_utf8_locale()

  attempt <- 0
  repeat {
    attempt <- attempt + 1
    result <- tryCatch(
      pangaear::pg_data(doi = doi, overwrite = overwrite, mssgs = FALSE),
      error = function(e) e
    )
    if (!inherits(result, "error")) break
    if (attempt > retries) {
      rlang::abort(sprintf(
        "Failed to fetch PANGAEA dataset %s after %d attempt(s): %s",
        doi, attempt, conditionMessage(result)
      ))
    }
    Sys.sleep(1)
  }

  out <- purrr::map_dfr(result, one_pangaea_entry_to_row)

  abort_if_mojibake(out$citation, doi)
  abort_if_mojibake(unlist(lapply(out$data, names), use.names = FALSE), doi)

  out
}
