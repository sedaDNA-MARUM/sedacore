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
#' An optional `also_matches` column lists further canonical names
#' (pipe-separated, e.g. `IRD_per_g|IRD_percent`) that count towards the
#' same checklist entry. Needed where one proxy is reported in several
#' genuinely different quantities (ice-rafted debris as counts, per
#' gram, percent or flux) that must stay distinct `parameter_name`s so
#' they're never interpolated into one another, but that a user checking
#' "do I have IRD for this core?" would treat as the same thing.
#'
#' @param path Path to `target_measurements.csv`.
#' @return A tibble: `parameter_name, notes`, plus `also_matches` when
#'   present in the file.
#' @export
load_target_measurements <- function(
    path = system.file("extdata", "target_measurements.csv", package = "sedacore")) {
  readr::read_csv(path, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
}

#' @noRd
expand_target_aliases <- function(targets) {
  also <- if ("also_matches" %in% names(targets)) targets$also_matches else rep(NA_character_, nrow(targets))
  purrr::map2_dfr(targets$parameter_name, also, function(p, a) {
    extra <- if (is.na(a) || !nzchar(a)) character() else trimws(strsplit(a, "|", fixed = TRUE)[[1]])
    tibble::tibble(target = p, name = unique(c(p, extra)))
  })
}

#' Every `parameter_name` the target checklist covers
#'
#' The checklist's own `parameter_name`s plus every `also_matches`
#' alternative -- e.g. for restricting a fetch to "the target checklist"
#' without accidentally excluding an IRD dataset because it reports
#' `IRD_per_g` rather than the checklist's headline `IRD_count`.
#'
#' @param targets A target-measurement checklist, from
#'   [load_target_measurements()].
#' @return A character vector of unique `parameter_name`s.
#' @export
target_measurement_names <- function(targets = load_target_measurements()) {
  unique(expand_target_aliases(targets)$name)
}

#' @noRd
summarise_target_coverage <- function(data, aliases) {
  empty <- tibble::tibble(
    parameter_name = character(), n_measured = integer(),
    depth_min_cm = numeric(), depth_max_cm = numeric(), source_core_id = character(),
    found_as = character()
  )
  # Fold each alternative name onto the checklist entry it counts towards
  # (keeping the original in found_as), so everything downstream just
  # groups by the entry's own parameter_name.
  data <- data[data$parameter_name %in% aliases$name, ]
  # dplyr::summarise() evaluates its expressions once even on a 0-row
  # input (to infer output types), which triggers min()/max()'s "no
  # non-missing arguments" warning despite the result being correctly
  # empty -- a realistic case here (a core with nothing fetched yet, or
  # none of the checklist's parameters), so short-circuit instead of just
  # living with the noise.
  if (nrow(data) == 0) {
    return(empty)
  }
  data$found_as <- data$parameter_name
  data$parameter_name <- aliases$target[match(data$parameter_name, aliases$name)]
  data |>
    dplyr::group_by(.data$parameter_name) |>
    dplyr::summarise(
      n_measured = dplyr::n(),
      depth_min_cm = min(.data$depth_cm, na.rm = TRUE),
      depth_max_cm = max(.data$depth_cm, na.rm = TRUE),
      source_core_id = paste(sort(unique(.data$core_id)), collapse = ", "),
      found_as = paste(sort(unique(.data$found_as)), collapse = ", "),
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
#'   `depth_max_cm` (`NA` when missing), `source_core_id` (the core(s)
#'   the reported data actually came from; `NA` when missing), and
#'   `found_as` (the actual `parameter_name`(s) that satisfied the entry --
#'   differs from `parameter_name` when an `also_matches` alternative
#'   did; `NA` when missing).
#' @export
check_target_measurements <- function(x, core_id, targets = load_target_measurements(),
                                       neighbour_core_ids = character()) {
  validate_param_long(x)

  aliases <- expand_target_aliases(targets)
  own_coverage <- summarise_target_coverage(x[x$core_id == core_id & !x$is_interpolated, ], aliases)

  neighbour_core_ids <- unique(neighbour_core_ids[!is.na(neighbour_core_ids) & neighbour_core_ids != core_id])
  neighbour_coverage <- summarise_target_coverage(
    x[x$core_id %in% neighbour_core_ids & !x$is_interpolated, ], aliases
  )

  out <- dplyr::left_join(targets, own_coverage, by = "parameter_name")
  need_neighbour <- is.na(out$n_measured) & out$parameter_name %in% neighbour_coverage$parameter_name
  if (any(need_neighbour)) {
    fill <- neighbour_coverage[match(out$parameter_name[need_neighbour], neighbour_coverage$parameter_name), ]
    out$n_measured[need_neighbour] <- fill$n_measured
    out$depth_min_cm[need_neighbour] <- fill$depth_min_cm
    out$depth_max_cm[need_neighbour] <- fill$depth_max_cm
    out$source_core_id[need_neighbour] <- fill$source_core_id
    out$found_as[need_neighbour] <- fill$found_as
  }

  out$status <- dplyr::case_when(
    !is.na(out$source_core_id) & out$source_core_id == core_id ~ "found",
    !is.na(out$source_core_id) ~ "found in neighbour",
    TRUE ~ "missing"
  )
  out$n_measured[is.na(out$n_measured)] <- 0L

  dplyr::select(out, "parameter_name", "notes", "status", "n_measured",
    "depth_min_cm", "depth_max_cm", "source_core_id", "found_as"
  )
}
