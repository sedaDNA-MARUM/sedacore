#' Load the target-measurement checklist
#'
#' A short, fixed list of parameters of particular interest for sedaDNA
#' work (e.g. TOC/TN/CaCO3 for productivity & preservation context,
#' age_model for chronology/sedimentation rate) -- edit
#' `inst/extdata/target_measurements.csv` to add/remove entries as the
#' list is refined; nothing else needs to change. `parameter_name` values
#' must match the canonical names produced by [harmonize_columns()] (see
#' `inst/extdata/column_synonyms.csv`).
#'
#' @param path Path to `target_measurements.csv`.
#' @return A tibble: `parameter_name, notes`.
#' @export
load_target_measurements <- function(
    path = system.file("extdata", "target_measurements.csv", package = "sedacore")) {
  readr::read_csv(path, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
}

#' @noRd
summarise_target_coverage <- function(data) {
  # dplyr::summarise() evaluates its expressions once even on a 0-row
  # input (to infer output types), which triggers min()/max()'s "no
  # non-missing arguments" warning despite the result being correctly
  # empty -- a realistic case here (a core with nothing fetched yet), so
  # short-circuit instead of just living with the noise.
  if (nrow(data) == 0) {
    return(tibble::tibble(
      parameter_name = character(), n_measured = integer(),
      depth_min_cm = numeric(), depth_max_cm = numeric(), source_core_id = character()
    ))
  }
  data |>
    dplyr::group_by(.data$parameter_name) |>
    dplyr::summarise(
      n_measured = dplyr::n(),
      depth_min_cm = min(.data$depth_cm, na.rm = TRUE),
      depth_max_cm = max(.data$depth_cm, na.rm = TRUE),
      source_core_id = paste(sort(unique(.data$core_id)), collapse = ", "),
      .groups = "drop"
    )
}

#' Check which target measurements are present for a core
#'
#' Cross-references a fixed checklist (see [load_target_measurements()])
#' against a core's own harmonized data, so it's immediately visible which
#' parameters still need to be sourced (fetched from PANGAEA, or measured
#' at MARUM) for that core, rather than having to scan the full
#' harmonized table by eye. Interpolated estimates never count -- only
#' directly-measured rows.
#'
#' If `neighbour_core_ids` is given, a parameter missing from the core's
#' own data but present in one of those neighbours is reported as
#' `"found in neighbour"` rather than `"missing"` -- distinct from a
#' genuine `"found"`, since it's an approximation from a nearby core, not
#' something actually measured on this one.
#'
#' @param x A `param_long` tibble.
#' @param core_id The core to check.
#' @param targets A target-measurement checklist, from
#'   [load_target_measurements()].
#' @param neighbour_core_ids Optional character vector of neighbouring
#'   cores' IDs (e.g. `neighbour_registry$neighbour_core_id`, filtered to
#'   `core_id`'s own search) whose data in `x` should be checked as a
#'   fallback for parameters `core_id` has none of its own.
#' @return A tibble: `parameter_name, notes, status` (`"found"`, `"found
#'   in neighbour"`, or `"missing"`), `n_measured`, `depth_min_cm`,
#'   `depth_max_cm` (`NA` when missing), and `source_core_id` (the core(s)
#'   the reported data actually came from; `NA` when missing).
#' @export
check_target_measurements <- function(x, core_id, targets = load_target_measurements(),
                                       neighbour_core_ids = character()) {
  validate_param_long(x)

  own_coverage <- summarise_target_coverage(x[x$core_id == core_id & !x$is_interpolated, ])

  neighbour_core_ids <- unique(neighbour_core_ids[!is.na(neighbour_core_ids) & neighbour_core_ids != core_id])
  neighbour_coverage <- summarise_target_coverage(
    x[x$core_id %in% neighbour_core_ids & !x$is_interpolated, ]
  )

  out <- dplyr::left_join(targets, own_coverage, by = "parameter_name")
  need_neighbour <- is.na(out$n_measured) & out$parameter_name %in% neighbour_coverage$parameter_name
  if (any(need_neighbour)) {
    fill <- neighbour_coverage[match(out$parameter_name[need_neighbour], neighbour_coverage$parameter_name), ]
    out$n_measured[need_neighbour] <- fill$n_measured
    out$depth_min_cm[need_neighbour] <- fill$depth_min_cm
    out$depth_max_cm[need_neighbour] <- fill$depth_max_cm
    out$source_core_id[need_neighbour] <- fill$source_core_id
  }

  out$status <- dplyr::case_when(
    !is.na(out$source_core_id) & out$source_core_id == core_id ~ "found",
    !is.na(out$source_core_id) ~ "found in neighbour",
    TRUE ~ "missing"
  )
  out$n_measured[is.na(out$n_measured)] <- 0L

  dplyr::select(out, "parameter_name", "notes", "status", "n_measured",
    "depth_min_cm", "depth_max_cm", "source_core_id"
  )
}
