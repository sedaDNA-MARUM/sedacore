#' Load the unit-conversion lookup table
#'
#' @param path Path to `unit_conversions.csv`.
#' @return A tibble: `parameter_name, unit_raw, unit_target, multiplier,
#'   offset`.
#' @export
load_unit_conversions <- function(
    path = system.file("extdata", "unit_conversions.csv", package = "sedacore")) {
  readr::read_csv(
    path,
    show_col_types = FALSE,
    col_types = readr::cols(
      parameter_name = "c", unit_raw = "c", unit_target = "c",
      multiplier = "d", offset = "d"
    )
  )
}

#' Rescale `param_long` values to a common unit per parameter
#'
#' Depth is already normalized to centimetres at harmonization time (see
#' [harmonize_columns()]); this handles the remaining case, where the same
#' `parameter_name` shows up with different raw units/scales across
#' datasets (e.g. susceptibility reported as `10**-5 SI` in one dataset and
#' `10**-6 SI` in another). Deliberately does **not** auto-convert isotope
#' reference scales (PDB vs VPDB): those should be kept as distinct
#' `parameter_name` values (e.g. `d18O_VPDB`) via `column_synonyms.csv`,
#' since silently converting between them is a real source of error in
#' this literature, not a unit-conversion problem.
#'
#' @param x A `param_long` tibble.
#' @param conversions A conversion table, from [load_unit_conversions()].
#' @param unmapped One of `"pass_through"` (default: leave `value`/`unit`
#'   unchanged), `"warn"`, or `"error"`: what to do about
#'   (`parameter_name`, `unit`) combinations not found in `conversions`.
#' @return `x`, with `value`/`unit` rescaled to each parameter's
#'   `unit_target`.
#' @export
standardize_units <- function(x, conversions = load_unit_conversions(),
                               unmapped = c("pass_through", "warn", "error")) {
  validate_param_long(x)
  unmapped <- match.arg(unmapped)

  key <- paste(x$parameter_name, x$unit, sep = "\r")
  conv_key <- paste(conversions$parameter_name, conversions$unit_raw, sep = "\r")
  idx <- match(key, conv_key)

  is_unmapped <- is.na(idx) & !is.na(x$value)
  if (any(is_unmapped)) {
    combos <- unique(paste0(x$parameter_name[is_unmapped], " [", x$unit[is_unmapped], "]"))
    msg <- sprintf(
      "No unit_conversions.csv entry for: %s. Values left unchanged.",
      paste(combos, collapse = ", ")
    )
    if (unmapped == "warn") rlang::warn(msg, class = "sedacore_unmapped_units")
    if (unmapped == "error") rlang::abort(msg, class = "sedacore_unmapped_units")
  }

  multiplier <- ifelse(is.na(idx), 1, conversions$multiplier[idx])
  offset <- ifelse(is.na(idx), 0, conversions$offset[idx])
  unit_target <- ifelse(is.na(idx), x$unit, conversions$unit_target[idx])

  x$value <- x$value * multiplier + offset
  x$unit <- unit_target
  x
}
