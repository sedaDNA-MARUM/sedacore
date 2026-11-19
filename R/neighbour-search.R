#' Normalize a core label for identity comparison
#'
#' PANGAEA event labels for the same core vary in formatting across
#' datasets -- e.g. `"GeoB18549-2"` in some, `"GeoB18549-2 (528-2)"`
#' (with a trailing cruise station number) in others for the exact same
#' core. Confirmed empirically: an exact-string comparison silently fails
#' to recognize these as the same core, which matters for excluding cores
#' the user already has from [find_neighbours()]'s results. Strips
#' anything from the first space or opening parenthesis onward, then
#' trims/lowercases.
#'
#' @param x A character vector of core IDs/labels.
#' @return `x`, normalized.
#' @noRd
normalize_core_label <- function(x) {
  tolower(trimws(stringr::str_extract(x, "^[^\\s(]+")))
}

#' Compute a bounding box circumscribing a search radius
#'
#' PANGAEA's search only accepts a rectangular `bbox`, not a circle, so
#' this over-fetches a square around `(latitude, longitude)` for
#' [find_neighbours()] to then cut down precisely with [circle_filter()].
#' Uses the same Mollweide equal-area projection already trusted for
#' distance calculations in the sibling `sedaDNA_atlantic_LGM` project's
#' `00a--coast_proximity_score.R`, for consistency.
#'
#' @param latitude,longitude Centre point, WGS84 degrees.
#' @param radius_km Search radius in kilometres.
#' @return A numeric vector `c(minlon, minlat, maxlon, maxlat)`, suitable
#'   for [pangaear::pg_search()]'s `bbox` argument.
#' @export
bbox_from_radius <- function(latitude, longitude, radius_km) {
  pt <- sf::st_sfc(sf::st_point(c(longitude, latitude)), crs = 4326)
  pt_moll <- sf::st_transform(pt, "+proj=moll")
  buf_moll <- sf::st_buffer(pt_moll, dist = radius_km * 1000)
  buf_wgs84 <- sf::st_transform(buf_moll, 4326)
  bb <- sf::st_bbox(buf_wgs84)
  c(unname(bb[["xmin"]]), unname(bb[["ymin"]]), unname(bb[["xmax"]]), unname(bb[["ymax"]]))
}

#' Cut a candidate table down to a precise radius
#'
#' `bbox_from_radius()`'s rectangle over-includes the corners of the
#' circumscribing square; this computes true great-circle distance
#' (Mollweide projection, matching `bbox_from_radius()`) and filters to
#' `distance_km <= radius_km`.
#'
#' @param x A data frame with `latitude`/`longitude` columns (WGS84
#'   degrees).
#' @param center_lat,center_lon Centre point, WGS84 degrees.
#' @param radius_km Search radius in kilometres.
#' @return `x`, filtered to rows within `radius_km`, with a `distance_km`
#'   column added. Rows with missing coordinates are dropped.
#' @export
circle_filter <- function(x, center_lat, center_lon, radius_km) {
  has_coords <- !is.na(x$latitude) & !is.na(x$longitude)
  x <- x[has_coords, ]
  if (nrow(x) == 0) {
    x$distance_km <- numeric()
    return(x)
  }

  pts <- sf::st_as_sf(x, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  pts_moll <- sf::st_transform(pts, "+proj=moll")
  center_moll <- sf::st_transform(
    sf::st_sfc(sf::st_point(c(center_lon, center_lat)), crs = 4326),
    "+proj=moll"
  )

  x$distance_km <- as.numeric(sf::st_distance(pts_moll, center_moll)) / 1000
  x[x$distance_km <= radius_km, ]
}

#' Find PANGAEA-indexed cores within a radius of a core
#'
#' Two-step search, since PANGAEA's search API filters by bounding box but
#' doesn't return per-hit coordinates: [bbox_from_radius()] casts a wide
#' net via [pangaear::pg_search()], then each candidate is fetched (via
#' [pg_fetch_dataset()], with `overwrite = FALSE` to make repeated
#' searches cheaper) to get its real coordinates for a precise
#' [circle_filter()] cut. This is the most request-heavy operation in the
#' package -- expect it to be slow for large radii/dense regions.
#'
#' @param x A `core_registry`, as produced by [collapse_to_core_registry()].
#'   Candidates matching *any* `core_id` already in `x` are excluded from
#'   the results, not just the one being searched around -- e.g.
#'   searching around `"GeoB18549-1"` also excludes `"GeoB18549-2"` if
#'   both are already in the user's own input file, since the point of
#'   this search is to discover cores not already on hand.
#' @param core_id The `core_id` in `x` to search around.
#' @param radius_km Search radius in kilometres.
#' @param count Max bbox candidates to request from PANGAEA before the
#'   precise circle cut.
#' @param progress Optional `function(i, n, doi)` called after each of the
#'   `n` bbox candidates is fetched (`i` is the 1-based candidate index,
#'   `doi` the one just fetched) -- lets a caller (e.g. a Shiny app via
#'   [shiny::withProgress()]) show live progress during what is otherwise
#'   a long sequence of individually-fast requests. Errors from `progress`
#'   itself are not caught.
#' @param source_programs Optional character vector restricting results to
#'   these drilling programs (any of `"GeoB"`, `"IODP"`, `"ODP"`, `"DSDP"`,
#'   `"other"`, as classified by [classify_core_source()]); `NULL`
#'   (default) keeps all programs. Applied as a post-filter on the
#'   candidates within `radius_km` -- every candidate in the bbox still
#'   has to be fetched to learn its coordinates in the first place (see
#'   Details), so this does not reduce the number of PANGAEA requests,
#'   only which of the results are returned.
#' @param query PANGAEA search query used for the bbox candidate search.
#'   Defaults to the quoted phrase `"Depth sed"`, present in essentially
#'   every individual sediment-core dataset -- confirmed empirically to
#'   matter a lot: a wildcard `"*"` query relies on PANGAEA's relevance
#'   ranking within the bbox, which is not distance-aware, and in
#'   well-studied regions (thousands of indexed datasets, most of them
#'   moorings/sensors/compilations rather than cores) the genuinely
#'   closest cores can rank far outside any practical `count`. Restricting
#'   to `"Depth sed"` cuts the candidate pool from thousands to hundreds
#'   in such regions and reliably surfaces real nearby cores instead.
#' @param required_parameters Optional character vector of
#'   `parameter_name`s (see [harmonize_columns()]/`column_synonyms.csv`);
#'   if given, only candidates with at least one of these parameters are
#'   returned (OR semantics, matching `source_programs`) -- e.g. "only
#'   show neighbours that actually have TOC data". `NULL` (default) keeps
#'   all candidates. This costs no extra PANGAEA requests: every
#'   candidate's data is already fetched to read its coordinates (see
#'   Details), so checking parameter availability just reuses that same
#'   in-memory data via [harmonize_columns()]. Regardless of whether this
#'   filter is set, every candidate's actually-found parameters are always
#'   recorded in the result's `available_parameters` column.
#' @return A `neighbour_registry` tibble, sorted with `"GeoB"` cores
#'   first, then `"IODP"`, then `"ODP"`/`"DSDP"`/`"other"`, and within each
#'   of those groups by `date_time` descending (most recently collected
#'   first) -- cores from MARUM's own program are the most immediately
#'   actionable for follow-up retrieval, and more recent data is more
#'   likely to be well-documented/complete.
#' @export
find_neighbours <- function(x, core_id, radius_km, count = 200, progress = NULL,
                             source_programs = NULL, query = "\"Depth sed\"",
                             required_parameters = NULL) {
  validate_core_registry(x)
  row <- x[x$core_id == core_id, ]
  if (nrow(row) != 1) {
    rlang::abort(sprintf("`core_id` \"%s\" must match exactly one row in `x`.", core_id))
  }

  empty <- empty_neighbour_registry()

  bbox <- bbox_from_radius(row$latitude, row$longitude, radius_km)
  candidates <- tryCatch(
    pangaear::pg_search(query = query, bbox = bbox, count = count),
    error = function(e) NULL
  )
  if (is.null(candidates) || nrow(candidates) == 0) {
    return(empty)
  }

  n <- nrow(candidates)
  coords <- purrr::imap_dfr(candidates$doi, function(doi, i) {
    fetched <- tryCatch(pg_fetch_dataset(doi, overwrite = FALSE), error = function(e) NULL)
    if (!is.null(progress)) {
      progress(i, n, doi)
    }
    if (is.null(fetched) || nrow(fetched) == 0) {
      return(tibble::tibble())
    }
    # Every candidate's raw data is already in memory at this point (it
    # had to be fetched to read off coordinates) -- reuse it here to
    # record which parameters this candidate actually has, at no extra
    # network cost. A candidate failing to harmonize (e.g. no "Depth
    # sed" column despite matching the bbox query) just means "no
    # parameters known", not a hard error for the whole search.
    fetched$available_parameters <- vapply(seq_len(nrow(fetched)), function(j) {
      h <- tryCatch(
        harmonize_columns(fetched$data[[j]], core_id = "candidate", doi = fetched$doi[j], unmatched = "silent"),
        error = function(e) NULL
      )
      if (is.null(h) || nrow(h) == 0) {
        return(NA_character_)
      }
      paste(sort(unique(h$parameter_name)), collapse = ", ")
    }, character(1))
    fetched[, c("doi", "event_label", "latitude", "longitude", "date_time", "citation", "available_parameters")]
  })
  if (nrow(coords) == 0) {
    return(empty)
  }

  filtered <- circle_filter(coords, row$latitude, row$longitude, radius_km)
  known_core_ids <- normalize_core_label(x$core_id)
  is_known <- !is.na(filtered$event_label) &
    normalize_core_label(filtered$event_label) %in% known_core_ids
  filtered <- filtered[!is_known, ]
  if (nrow(filtered) == 0) {
    return(empty)
  }

  out <- tibble::tibble(
    core_id = core_id,
    neighbour_doi = filtered$doi,
    neighbour_core_id = filtered$event_label,
    latitude = filtered$latitude,
    longitude = filtered$longitude,
    distance_km = filtered$distance_km,
    source_program = classify_core_source(filtered$event_label),
    date_time = filtered$date_time,
    citation = filtered$citation,
    available_parameters = filtered$available_parameters,
    included_in_export = TRUE
  )

  if (!is.null(source_programs)) {
    out <- out[out$source_program %in% source_programs, ]
  }
  if (!is.null(required_parameters)) {
    has_required <- vapply(out$available_parameters, function(ap) {
      if (is.na(ap)) {
        return(FALSE)
      }
      any(trimws(strsplit(ap, ",")[[1]]) %in% required_parameters)
    }, logical(1))
    out <- out[has_required, ]
  }
  if (nrow(out) == 0) {
    return(empty)
  }

  new_neighbour_registry(order_neighbours(out))
}

#' Sort a neighbour candidate table by program preference, then recency
#'
#' `"GeoB"` cores first, then `"IODP"`, then `"ODP"`/`"DSDP"`/`"other"`;
#' within each of those groups, most recently collected (`date_time`)
#' first. Rows with a missing `date_time` sort last within their program
#' group. Extracted as its own function so the ordering logic is directly
#' testable without a live PANGAEA search.
#'
#' @param x A data frame with `source_program` and `date_time` columns
#'   (e.g. the `out` table being assembled inside [find_neighbours()]).
#' @return `x`, reordered.
#' @noRd
order_neighbours <- function(x) {
  program_rank <- c(GeoB = 1, IODP = 2, ODP = 3, DSDP = 4, other = 5)
  x[order(program_rank[x$source_program], -as.numeric(x$date_time)), ]
}
