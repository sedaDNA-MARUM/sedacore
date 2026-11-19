# Logic that used to live inline in inst/app/app.R, where it could not be
# tested without driving a browser. Everything here is a pure function of
# its arguments (no Shiny, no reactive state, no network), so the app is
# reduced to wiring: read inputs, call these, write the results back.
#
# Exported (but hidden from the index via @keywords internal) because the
# app runs as a separate script that only sees the package's exports.

# -- empty tables, derived from the schema rather than copied from it --------

#' @noRd
empty_from_spec <- function(spec) {
  cols <- lapply(spec, function(type) {
    switch(type,
      character = character(),
      numeric = numeric(),
      integer = integer(),
      logical = logical(),
      POSIXct = as.POSIXct(character()),
      rlang::abort(sprintf("No empty vector defined for schema type \"%s\".", type))
    )
  })
  tibble::as_tibble(cols)
}

#' Zero-row tables used to initialise app state
#'
#' Built from [schema_spec][new_param_long] so they cannot drift from the
#' schema (the app used to carry hand-copied duplicates).
#'
#' @return A zero-row tibble of the named shape.
#' @keywords internal
#' @export
empty_param_long <- function() {
  new_param_long(empty_from_spec(schema_spec$param_long))
}

#' @rdname empty_param_long
#' @export
empty_neighbour_registry <- function() {
  new_neighbour_registry(empty_from_spec(schema_spec$neighbour_registry))
}

#' @rdname empty_param_long
#' @details `empty_depth_matched_long()` is [empty_param_long()] plus
#'   `target_core_id`: which active core's sample depths a row was matched
#'   against. Usually the row's own `core_id`, but different when the
#'   measurement was *borrowed* from a neighbour -- kept as its own column
#'   so that is always visible rather than implied.
#' @export
empty_depth_matched_long <- function() {
  x <- empty_param_long()
  x$target_core_id <- character()
  x
}

#' @rdname empty_param_long
#' @export
empty_search_results <- function() {
  tibble::tibble(
    core_id = character(), score = numeric(), doi = character(),
    size = numeric(), size_measure = character(), citation = character(),
    available_parameters = character()
  )
}

#' @rdname empty_param_long
#' @export
empty_zenodo_search_results <- function() {
  tibble::tibble(
    doi = character(), title = character(), description = character(), n_files = integer()
  )
}

# -- measurement filter -------------------------------------------------------

#' Restrict harmonized rows to the selected measurements
#'
#' Applies the "restrict to selected measurements" choice, except that
#' `age_model` rows are **always kept**. An age model is not a measurement
#' one opts in or out of: it is context needed (a) to show an `age_ka`
#' column next to whatever is selected, and (b) by
#' [report_depth_coverage()]'s age-based neighbour matching, which has
#' nothing to work with if a core's own age model was filtered away
#' upstream of it.
#'
#' @param h A `param_long` tibble.
#' @param selected Character vector of `parameter_name`s to keep; empty or
#'   `NULL` keeps everything.
#' @return `h`, restricted.
#' @keywords internal
#' @export
keep_selected_measurements <- function(h, selected) {
  if (length(selected) == 0) {
    return(h)
  }
  h[h$parameter_name %in% selected | h$parameter_name == "age_model", ]
}

#' Harmonize one raw table and apply the measurement filter
#'
#' The orchestration the app used to do inline: [harmonize_columns()], the
#' optional measurement filter (applied *after* harmonizing so nothing is
#' lost during parsing itself, only when deciding what to keep), and the
#' log lines describing what happened.
#'
#' @param raw_data One raw data table (a `data` list-column element).
#' @param core_id,doi Passed to [harmonize_columns()].
#' @param depth_col,synonym_map Source-specific depth column and synonym
#'   table (PANGAEA vs IODP/Zenodo conventions differ).
#' @param measurement_filter Character vector of `parameter_name`s, or empty.
#' @param synonym_csv_name Only used in the "not in <csv>" log message.
#' @return A list: `data` (the rows to append, or `NULL` when there is
#'   nothing to add) and `messages` (log lines, in the order to emit them).
#' @keywords internal
#' @export
harmonize_for_store <- function(raw_data, core_id, doi, depth_col, synonym_map,
                                measurement_filter = NULL, synonym_csv_name = "column_synonyms.csv") {
  h <- tryCatch(
    harmonize_columns(
      raw_data,
      core_id = core_id, doi = doi,
      depth_col = depth_col, synonym_map = synonym_map, unmatched = "warn"
    ),
    error = function(e) e
  )
  if (inherits(h, "error")) {
    return(list(data = NULL, messages = paste("ERROR harmonizing:", conditionMessage(h))))
  }
  if (nrow(h) == 0) {
    return(list(data = NULL, messages = character()))
  }

  fallback <- attr(h, "fallback_columns")
  dropped <- attr(h, "dropped_columns")
  n_before_filter <- nrow(h)
  kept <- keep_selected_measurements(h, measurement_filter)

  messages <- if (nrow(kept) > 0) {
    sprintf("Harmonized %d rows from %s.", nrow(kept), doi)
  } else {
    sprintf(
      "0 rows from %s after restricting to selected measurements (%d rows available before filtering).",
      doi, n_before_filter
    )
  }
  if (length(fallback) > 0) {
    messages <- c(messages, sprintf(
      "Not in %s, imported as-is: %s", synonym_csv_name, paste(fallback, collapse = ", ")
    ))
  }
  if (length(dropped) > 0) {
    messages <- c(messages, sprintf(
      "Excluded (non-numeric, e.g. free-text comment fields): %s", paste(dropped, collapse = ", ")
    ))
  }

  list(data = if (nrow(kept) > 0) kept else NULL, messages = messages)
}

# -- display tables -----------------------------------------------------------

#' Show `age_model` as an `age_ka` column instead of as measurement rows
#'
#' Drops the raw `age_model` rows and gives every other row the age
#' interpolated from *its own core's* age model at its depth, placed right
#' after `depth_cm`. Each core's age model is extracted once (the app used
#' to rescan the whole table per core on every render). Values already
#' present in an existing `age_ka` column are kept; only gaps are filled.
#'
#' @param x A `param_long` tibble.
#' @param gap_multiplier,snap_tolerance_cm Passed to [depth_to_age()].
#' @return `x` without `age_model` rows, plus `age_ka`.
#' @keywords internal
#' @export
attach_age_ka <- function(x, gap_multiplier = 10, snap_tolerance_cm = 4) {
  validate_param_long(x)
  age_rows <- x[x$parameter_name == "age_model" & !is.na(x$value), ]
  out <- x[x$parameter_name != "age_model", ]
  if (!"age_ka" %in% names(out)) {
    out$age_ka <- NA_real_
  }

  models <- split(age_rows, age_rows$core_id)
  for (cid in intersect(names(models), unique(out$core_id))) {
    m <- models[[cid]]
    m <- m[order(m$depth_cm), ]
    age_model <- tibble::tibble(depth_cm = m$depth_cm, age_ka = m$value)
    idx <- out$core_id == cid
    out$age_ka[idx] <- dplyr::coalesce(
      out$age_ka[idx],
      depth_to_age(age_model, out$depth_cm[idx], gap_multiplier, snap_tolerance_cm)
    )
  }
  dplyr::relocate(out, "age_ka", .after = "depth_cm")
}

#' @noRd
collapse_ids <- function(ids, max_ids) {
  u <- sort(unique(ids))
  if (length(u) <= max_ids) {
    return(paste(u, collapse = ", "))
  }
  paste0(paste(u[seq_len(max_ids)], collapse = ", "), sprintf(" ... (+%d more)", length(u) - max_ids))
}

#' Summarise what has been retrieved, per measurement type
#'
#' One row per (`parameter_group`, `parameter_name`, `unit`). The `cores`
#' column is truncated to `max_cores` ids -- with hundreds of cores loaded
#' an unbounded comma list made a single cell megabytes wide.
#'
#' @param x A `param_long` tibble.
#' @param max_cores How many core ids to list before "... (+N more)".
#' @return A tibble: `parameter_group, parameter_name, unit, n_cores, cores,
#'   n_measured, n_interpolated`.
#' @keywords internal
#' @export
summarise_measurements <- function(x, max_cores = 6) {
  x |>
    dplyr::group_by(.data$parameter_group, .data$parameter_name, .data$unit) |>
    dplyr::summarise(
      n_cores = dplyr::n_distinct(.data$core_id),
      cores = collapse_ids(.data$core_id, max_cores),
      n_measured = sum(!.data$is_interpolated),
      n_interpolated = sum(.data$is_interpolated),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$parameter_group, .data$parameter_name)
}

# -- interpolation ------------------------------------------------------------

#' Select the data an interpolation run should use
#'
#' The active core's own data plus any already-fetched neighbours found for
#' it (some cores have no data of their own, so their neighbours' series are
#' the only thing to match against), restricted by the measurement filter
#' independently of anything applied at harmonize time -- so the filter can
#' be changed freely without re-fetching.
#'
#' @param param_long A `param_long` tibble (everything harmonized so far).
#' @param neighbour_results A `neighbour_registry`.
#' @param core_id The active (target) core.
#' @param measurement_filter Character vector of `parameter_name`s, or empty.
#' @return A list: `data`, `neighbour_ids`, and `message` -- non-`NULL` (a
#'   ready-to-log sentence) when `data` is empty.
#' @keywords internal
#' @export
interp_candidate_data <- function(param_long, neighbour_results, core_id, measurement_filter = NULL) {
  neighbour_ids <- unique(neighbour_results$neighbour_core_id[neighbour_results$core_id == core_id])
  neighbour_ids <- neighbour_ids[!is.na(neighbour_ids)]
  current <- param_long[param_long$core_id %in% union(core_id, neighbour_ids), ]
  current <- keep_selected_measurements(current, measurement_filter)

  message <- NULL
  if (nrow(current) == 0) {
    message <- "No harmonized data yet for the active core"
    if (length(neighbour_ids) > 0) message <- paste0(message, " or its fetched neighbours")
    if (length(measurement_filter) > 0) message <- paste0(message, " matching the selected measurement filter")
    message <- paste0(message, "; fetch some first.")
  }
  list(data = current, neighbour_ids = neighbour_ids, message = message)
}

#' Log lines summarising a depth/age-matching report
#'
#' Turns the output of [report_depth_coverage()] into the sentences shown
#' in the app's log: an overall tally, a depth-mismatch explanation, and --
#' because cross-core integration is refused outright rather than silently
#' depth-matched -- exactly *why* each neighbour (or sample depth) was
#' skipped, with a pointer to XRF data when that is the only alternative.
#'
#' @param report Output of [report_depth_coverage()].
#' @param core_id The active (target) core.
#' @param sample_depths The target's sample depths (cm).
#' @return A character vector, in the order to log them.
#' @keywords internal
#' @export
depth_coverage_messages <- function(report, core_id, sample_depths) {
  n_measured <- sum(report$match_status == "measured")
  n_interp <- sum(report$match_status == "interpolated")
  n_no_data <- sum(report$match_status == "no data in that depth range")
  n_gap <- sum(report$match_status == "gap too large to interpolate")
  n_age <- sum(report$match_basis == "age", na.rm = TRUE)
  n_refused <- sum(report$match_basis == "refused", na.rm = TRUE)

  borrowed_ids <- unique(report$core_id[!is.na(report$core_id) & report$core_id != core_id])
  source_note <- if (n_age > 0) {
    sprintf(" (%d borrowed from neighbour core(s) by age: %s)", n_age, paste(borrowed_ids, collapse = ", "))
  } else {
    ""
  }
  messages <- sprintf(
    "Depth coverage for %s's %d sample depth(s) x %d parameter(s): %d measured, %d interpolated, %d out of retrieved range, %d gaps too large%s.",
    core_id, length(unique(sample_depths)), length(unique(report$parameter_name)),
    n_measured, n_interp, n_no_data, n_gap, source_note
  )

  if (n_no_data > 0) {
    no_data <- report[report$match_status == "no data in that depth range", ]
    no_data <- unique(no_data[, c("parameter_group", "parameter_name", "measured_depth_min", "measured_depth_max")])
    details <- sprintf(
      "%s (%s): retrieved %.0f-%.0f cm",
      no_data$parameter_name, no_data$parameter_group,
      no_data$measured_depth_min, no_data$measured_depth_max
    )
    messages <- c(messages, sprintf(
      "Depth mismatch: your samples are at %.0f-%.0f cm, but retrieved data doesn't cover that range for: %s.",
      min(sample_depths, na.rm = TRUE), max(sample_depths, na.rm = TRUE), paste(details, collapse = "; ")
    ))
  }

  if (n_refused > 0) {
    # match_basis is NA for cells no path ever attempted; a bare
    # `== "refused"` index would pull all-NA rows in alongside the real
    # ones and make the any()/if below NA (this crashed the app once).
    refused <- report[!is.na(report$match_basis) & report$match_basis == "refused", ]
    no_age_model <- unique(refused$core_id[refused$match_status == "no age model for cross-core match"])
    no_age_model <- no_age_model[!is.na(no_age_model)]
    for (nb in no_age_model) {
      has_xrf <- any(refused$has_xrf_data[refused$core_id == nb], na.rm = TRUE)
      xrf_note <- if (has_xrf) {
        " XRF data is available for this core if you want to build a depth-time correlation manually."
      } else {
        ""
      }
      messages <- c(messages, sprintf(
        "Skipped %s: no age model available for cross-core matching.%s", nb, xrf_note
      ))
    }
    no_target_age <- "target core has no age_ka for this sample depth"
    if (any(refused$match_status == no_target_age, na.rm = TRUE)) {
      messages <- c(messages, sprintf(
        "Skipped cross-core matching for %d sample depth(s): no age_ka on file for %s in sample_registry.",
        length(unique(refused$depth_cm[refused$match_status == no_target_age])), core_id
      ))
    }
  }

  messages
}

# -- small conveniences -------------------------------------------------------

#' Turn a DOI into a clickable link
#'
#' For DT columns rendered with `escape = FALSE`. NA-safe (`NA` stays `NA`,
#' not the text "NA"). `base_url` defaults to PANGAEA's own landing-page
#' resolver; pass `"https://doi.org/"` for any other DOI (e.g. Zenodo's).
#'
#' @param doi Character vector of DOIs.
#' @param base_url Resolver prefix.
#' @return A character vector of HTML anchors.
#' @keywords internal
#' @export
doi_link <- function(doi, base_url = "https://doi.pangaea.de/") {
  ifelse(is.na(doi), NA_character_,
    sprintf(
      '<a href="%s%s" target="_blank" rel="noopener noreferrer">%s</a>',
      base_url, doi, doi
    )
  )
}

#' IODP expedition number encoded in a core id
#'
#' Zenodo's IODP records are per expedition, not per core, so a search by
#' core id is meaningless -- but an id like `"362-U1480A"` carries the
#' expedition number, matching the convention [classify_core_source()]
#' already recognises.
#'
#' @param core_id Character vector.
#' @return The expedition number as a string, `NA` where `core_id` is not an
#'   expedition-prefixed IODP hole.
#' @keywords internal
#' @export
iodp_expedition_from_core_id <- function(core_id) {
  stringr::str_match(core_id, "^([0-9]{2,3})-U[0-9]{4}[A-Z]?$")[, 2]
}

#' File name for an exported table
#'
#' @param core_id The active core (sanitised for use in a file name);
#'   missing/empty becomes `"export"`.
#' @param kind `"harmonized"` or `"depth_matched"`.
#' @param date Date stamp.
#' @return A file name ending in `.tsv`.
#' @keywords internal
#' @export
export_filename <- function(core_id, kind = c("harmonized", "depth_matched"), date = Sys.Date()) {
  kind <- match.arg(kind)
  core_part <- if (!is.null(core_id) && length(core_id) == 1 && !is.na(core_id) && nzchar(core_id)) {
    gsub("[^A-Za-z0-9_-]+", "_", core_id)
  } else {
    "export"
  }
  prefix <- if (kind == "harmonized") "sedacore" else "sedacore_depth_matched"
  sprintf("%s_%s_%s.tsv", prefix, core_part, format(date))
}

#' Read a core metadata TSV into both registries
#'
#' @param path Path to a core-metadata TSV (see the bundled
#'   `inst/extdata/example_core_metadata.tsv`).
#' @return A list: `sample_registry` and `core_registry`.
#' @keywords internal
#' @export
load_core_registries <- function(path) {
  raw <- readr::read_tsv(path, show_col_types = FALSE)
  sr <- as_sample_registry(raw)
  list(sample_registry = sr, core_registry = collapse_to_core_registry(sr))
}

#' Display form of the depth/age-matched table
#'
#' Hides `age_model` rows (their information is the `age_ka` column every
#' other row already carries) and shows `matched_age_ka` as `age_ka`
#' right after `depth_cm`. Works on the empty initial table too, which has
#' no `matched_age_ka` column yet (the app used to error when that tab was
#' opened before the first interpolation).
#'
#' @param x A `depth_matched_long` table (see [empty_depth_matched_long()]).
#' @return `x` reshaped for display.
#' @keywords internal
#' @export
depth_matched_display <- function(x) {
  x <- x[x$parameter_name != "age_model", ]
  if ("matched_age_ka" %in% names(x)) {
    x <- dplyr::rename(x, age_ka = "matched_age_ka")
  } else {
    x$age_ka <- NA_real_
  }
  dplyr::relocate(x, "age_ka", .after = "depth_cm")
}

#' One-sentence summary of a target checklist
#'
#' @param checklist Output of [check_target_measurements()].
#' @param core_id The core it was computed for.
#' @return A single string.
#' @keywords internal
#' @export
checklist_summary <- function(checklist, core_id) {
  missing <- checklist$parameter_name[checklist$status == "missing"]
  in_neighbour <- checklist$parameter_name[checklist$status == "found in neighbour"]
  parts <- character()
  if (length(in_neighbour) > 0) {
    parts <- c(parts, sprintf("found only in a neighbour core: %s", paste(in_neighbour, collapse = ", ")))
  }
  if (length(missing) > 0) {
    parts <- c(parts, sprintf("missing entirely: %s", paste(missing, collapse = ", ")))
  }
  if (length(parts) == 0) {
    sprintf("All %d target measurement(s) found for %s.", nrow(checklist), core_id)
  } else {
    sprintf("For %s -- %s.", core_id, paste(parts, collapse = "; "))
  }
}
