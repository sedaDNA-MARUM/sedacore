#' Extract a core's own age-depth model from a param_long tibble
#'
#' An age model isn't a separate table in sedacore's schema -- it's just
#' ordinary `param_long` rows with `parameter_name == "age_model"`
#' (`depth_cm`, `value` = age in ka), harmonized via the same
#' PANGAEA/Zenodo/manual-import pipeline as any other measurement (see
#' `column_synonyms.csv`'s `"Age model [ka]"`/`"Age [ka BP]"` rows). This
#' pulls just one core's own such rows out, sorted by depth, ready for
#' [depth_to_age()].
#'
#' @param x A `param_long` tibble.
#' @param core_id The core whose age model to extract.
#' @return A tibble with `depth_cm` and `age_ka`, sorted by `depth_cm`.
#'   Zero rows if this core has no `age_model` data -- callers should
#'   branch on `nrow() == 0` rather than expect an error.
#' @export
core_age_model <- function(x, core_id) {
  validate_param_long(x)
  rows <- x[x$core_id == core_id & x$parameter_name == "age_model" & !is.na(x$value), ]
  rows <- rows[order(rows$depth_cm), ]
  tibble::tibble(depth_cm = rows$depth_cm, age_ka = rows$value)
}

#' Convert depths to ages using a core's own age-depth model
#'
#' Exact match, then snap to a close-enough control point, then
#' gap-guarded linear interpolation between the two bracketing control
#' points -- the same three-tier logic as `interpolate_one_group()`
#' (`R/interpolate.R`), applied to a single age-model series instead of a
#' group of measurements. Never extrapolates beyond the age model's own
#' depth range.
#'
#' @param age_model A tibble from [core_age_model()].
#' @param depth_cm Numeric vector of depths to convert.
#' @param gap_multiplier Multiple of the age model's own median control-
#'   point spacing allowed between the two points bracketing a depth.
#' @param snap_tolerance_cm Maximum distance (cm) from a depth to an age-
#'   model control point for that point's age to be used as-is rather than
#'   requiring an exact depth match or falling through to interpolation.
#' @return A numeric vector of ages (ka), the same length as `depth_cm`,
#'   `NA` where no age could be resolved (empty age model, a depth outside
#'   its range, or bracketing control points too far apart).
#' @export
depth_to_age <- function(age_model, depth_cm, gap_multiplier = 10, snap_tolerance_cm = 4) {
  if (nrow(age_model) == 0) {
    return(rep(NA_real_, length(depth_cm)))
  }

  depths <- age_model$depth_cm
  ages <- age_model$age_ka
  spacing <- if (length(unique(depths)) >= 2) stats::median(diff(sort(unique(depths)))) else Inf
  max_gap <- gap_multiplier * spacing

  vapply(depth_cm, function(td) {
    # Age-only series (no depth) carry NA depth_cm; nothing to convert.
    if (is.na(td)) {
      return(NA_real_)
    }
    exact <- which(depths == td)
    if (length(exact) > 0) {
      return(ages[exact[1]])
    }

    gaps <- abs(depths - td)
    best <- which.min(gaps)
    if (gaps[best] <= snap_tolerance_cm) {
      return(ages[best])
    }

    if (length(depths) < 2 || td < min(depths) || td > max(depths)) {
      return(NA_real_)
    }
    lo <- max(depths[depths <= td])
    hi <- min(depths[depths >= td])
    if ((hi - lo) > max_gap) {
      return(NA_real_)
    }
    # round() to clear floating-point noise from the interpolation
    # arithmetic itself -- see the matching rationale in interpolate.R.
    round(stats::approx(depths, ages, xout = td)$y, 6)
  }, numeric(1))
}
