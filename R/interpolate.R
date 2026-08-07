#' Join sample depths to the nearest existing measurement
#'
#' Generalizes the `cdata_filtering()` nearest-depth join in
#' `00a--gather_core_contextual_data.R` (which hand-coded a ladder of
#' +/-1/2/3 cm offsets): for each `(core_id, depth_cm)` in `targets`, finds
#' the nearest measurement per `(parameter_group, parameter_name,
#' species)` in `x` within `tolerance_cm`, rather than synthesizing a new
#' value. Use [interpolate_depth_series()] instead when you want an
#' estimated value *between* measurements, not just the closest existing
#' one.
#'
#' @param x A `param_long` tibble of measured values.
#' @param targets A data frame with `core_id` and `depth_cm` columns (e.g.
#'   a `sample_registry`) giving the depths to snap to.
#' @param tolerance_cm Maximum distance (cm) to snap across.
#' @return A tibble: `x`'s columns plus `sample_depth_cm` (the requested
#'   target depth) and `gap_cm` (distance from the target to the matched
#'   measurement). One row per `(target depth, parameter_group,
#'   parameter_name, species)` combination found within tolerance.
#' @export
find_neighbours_by_depth <- function(x, targets, tolerance_cm) {
  validate_param_long(x)

  purrr::pmap_dfr(
    list(core_id = targets$core_id, target_depth_cm = targets$depth_cm),
    function(core_id, target_depth_cm) {
      candidates <- x[x$core_id == core_id & !is.na(x$value), ]
      if (nrow(candidates) == 0) {
        return(tibble::tibble())
      }
      candidates$gap_cm <- abs(candidates$depth_cm - target_depth_cm)
      candidates <- candidates[candidates$gap_cm <= tolerance_cm, ]
      if (nrow(candidates) == 0) {
        return(tibble::tibble())
      }

      candidates |>
        dplyr::group_by(.data$parameter_group, .data$parameter_name, .data$species) |>
        dplyr::slice_min(.data$gap_cm, n = 1, with_ties = FALSE) |>
        dplyr::ungroup() |>
        dplyr::mutate(sample_depth_cm = target_depth_cm)
    }
  )
}

#' @noRd
interpolate_one_group <- function(g, target_depths, method, gap_multiplier, snap_tolerance_cm) {
  # A row can carry a real depth_cm but an NA value -- e.g. PANGAEA tables
  # that share one composite/spliced depth column across many parameter
  # columns at fine resolution, where this particular parameter's own cell
  # is blank at that depth (measured only every few cm, not at every row
  # of the shared scale). Such a row is not a real measurement of this
  # group's parameter and must be dropped before any match/bracket logic,
  # or it gets picked up as a false "exact match" with value = NA -- which
  # both fabricates a bogus measured row and (worse) blocks the real
  # interpolation/snap to the actual nearby measurements from ever being
  # attempted, since the target depth looks "already resolved".
  g <- g[!is.na(g$value), ]
  if (nrow(g) == 0) {
    return(tibble::tibble())
  }

  exact_rows <- g[g$depth_cm %in% target_depths, ]
  # An exact match isn't bracketed by two different real depths -- it *is*
  # one, so both bounds collapse to the same value. Together with
  # depth_cm (the target/sample depth) and core_id (which core this row's
  # real measurement actually came from -- possibly a neighbour, not the
  # target core itself), this makes every row fully traceable back to the
  # real source data it came from.
  if (nrow(exact_rows) > 0) {
    exact_rows$bracket_depth_min_cm <- exact_rows$depth_cm
    exact_rows$bracket_depth_max_cm <- exact_rows$depth_cm
  }

  remaining_targets <- setdiff(target_depths, g$depth_cm)

  # Sample depths and a core's own measured depths are rarely recorded to
  # identical precision (e.g. a sample logged at 246 cm vs. the nearest
  # real measurement at 245.5 cm) -- snap to the closest real value within
  # snap_tolerance_cm and use it as-is (not interpolated) rather than
  # requiring exact equality or interpolating across a bracket needlessly.
  # Tried before interpolation: a very close real value is more reliable
  # than a linear estimate, and this also rescues target depths near the
  # edge of the measured range where interpolation can't reach.
  snapped_rows <- tibble::tibble()
  if (length(remaining_targets) > 0 && snap_tolerance_cm > 0) {
    snapped_rows <- purrr::map_dfr(remaining_targets, function(td) {
      gaps <- abs(g$depth_cm - td)
      best <- which.min(gaps)
      if (gaps[best] > snap_tolerance_cm) {
        return(tibble::tibble())
      }
      row <- g[best, ]
      row$bracket_depth_min_cm <- row$depth_cm
      row$bracket_depth_max_cm <- row$depth_cm
      row$depth_cm <- td
      row
    })
    if (nrow(snapped_rows) > 0) {
      remaining_targets <- setdiff(remaining_targets, snapped_rows$depth_cm)
    }
  }

  if (length(remaining_targets) == 0 || nrow(g) < 2) {
    return(dplyr::bind_rows(exact_rows, snapped_rows))
  }

  ord <- order(g$depth_cm)
  depths <- g$depth_cm[ord]
  values <- g$value[ord]
  spacing <- stats::median(diff(sort(unique(depths))))
  max_gap <- gap_multiplier * spacing

  interpolated_rows <- purrr::map_dfr(remaining_targets, function(td) {
    if (td < min(depths) || td > max(depths)) {
      return(tibble::tibble())
    }
    lo <- max(depths[depths <= td])
    hi <- min(depths[depths >= td])
    if ((hi - lo) > max_gap) {
      return(tibble::tibble())
    }
    val <- if (method == "linear") {
      stats::approx(depths, values, xout = td)$y
    } else {
      stats::spline(depths, values, xout = td)$y
    }
    # round() to clear floating-point noise from the interpolation
    # arithmetic itself (e.g. 16.740000000000002) -- see the matching
    # depth_cm rounding in harmonize_columns() for the same rationale.
    val <- round(val, 6)
    tibble::tibble(
      core_id = g$core_id[1], pangaea_doi = NA_character_, depth_cm = td,
      parameter_group = g$parameter_group[1], parameter_name = g$parameter_name[1],
      species = g$species[1], value = val, unit = g$unit[1],
      raw_column_name = NA_character_, is_interpolated = TRUE,
      interp_gap_cm = hi - lo, retrieved_at = Sys.time(),
      bracket_depth_min_cm = lo, bracket_depth_max_cm = hi
    )
  })

  dplyr::bind_rows(exact_rows, snapped_rows, interpolated_rows)
}

#' Interpolate `param_long` values at target depths
#'
#' Builds a **separate, depth-aligned view** of `x` at exactly the
#' requested `targets` depths (e.g. sedaDNA sample depths) -- it does not
#' modify or filter `x` itself, so the full raw/measured dataset in `x`
#' stays intact and nothing already retrieved is ever lost; this is purely
#' an additional, purpose-built table for depth-matching against sample
#' depths. For each `(core_id, parameter_group, parameter_name, species)`
#' group: a target depth that was already directly measured is returned
#' as-is (`is_interpolated = FALSE`); otherwise, if a real measurement
#' exists within `snap_tolerance_cm` of the target depth, that nearby
#' value is used as-is instead (also `is_interpolated = FALSE` -- sample
#' depths and a core's own measured depths are rarely recorded to
#' identical precision, so requiring exact equality would miss real,
#' usably-close data); otherwise a value is estimated between the two
#' bracketing measurements, subject to a **per-core gap-size guard** --
#' only filled if those two measurements are within `gap_multiplier`
#' times that group's own median sample spacing of each other (not a
#' single hardcoded cm threshold, unlike the +/-1/2/3 cm ladder in the
#' original `cdata_filtering()`) -- this avoids interpolating across a
#' genuine data desert in sparsely-sampled cores. Never extrapolates
#' beyond the measured depth range; a target depth outside that range, or
#' with too large a bracketing gap, is simply absent from the result.
#' Rows whose `value` is `NA` (e.g. a depth shared with another parameter
#' on a composite/spliced depth scale, where this parameter itself wasn't
#' measured there) are never treated as real measurements -- excluded
#' from exact/snap matching and interpolation alike.
#'
#' Every interpolated row is stamped `is_interpolated = TRUE` with
#' `interp_gap_cm` recording the bracketing gap, so downstream analyses
#' can trivially `dplyr::filter(!is_interpolated)` to get only the
#' directly-measured subset of this depth-matched table. Two more columns,
#' `bracket_depth_min_cm`/`bracket_depth_max_cm`, record the *real* source
#' depth(s) a row actually came from -- for an exact match these both
#' equal `depth_cm` itself; for an interpolated estimate they're the two
#' real bracketing measurements it was computed between. Combined with
#' `depth_cm` (the target/sample depth being matched to) and `core_id`
#' (which core that real measurement actually came from -- e.g. a
#' neighbour, when depth-matching borrowed data), every row is fully
#' traceable: target depth, in which core, matched from which real
#' source depth(s), in which core.
#'
#' @param x A `param_long` tibble of measured values (rows already flagged
#'   `is_interpolated = TRUE` are ignored as interpolation inputs).
#' @param targets A data frame with `core_id` and `depth_cm` columns (e.g.
#'   a `sample_registry`) giving the depths to estimate.
#' @param method `"linear"` (default) or `"spline"`. Linear is the
#'   conservative default: sediment depth series are often sparse and
#'   irregular, and splines are prone to overshoot between sparse points
#'   in a way that's misleading for physically bounded values like
#'   isotope ratios.
#' @param gap_multiplier Multiple of a group's own median sample spacing
#'   allowed between the two measurements bracketing a target depth.
#' @param snap_tolerance_cm Maximum distance (cm) from a target depth to a
#'   real measurement for that measurement to be used as-is rather than
#'   requiring an exact depth match or falling through to interpolation.
#' @return A `param_long` tibble with one row per `(target depth,
#'   parameter_group, parameter_name, species)` that could be resolved --
#'   an exact/snapped measurement or an interpolated estimate -- plus
#'   `bracket_depth_min_cm`/`bracket_depth_max_cm` (see Details). This is
#'   a new, separate table; it is not merged into `x`.
#' @export
interpolate_depth_series <- function(x, targets,
                                      method = c("linear", "spline"),
                                      gap_multiplier = 10,
                                      snap_tolerance_cm = 4) {
  validate_param_long(x)
  method <- match.arg(method)

  measured <- x[!x$is_interpolated, ]
  target_depths_by_core <- split(targets$depth_cm, targets$core_id)

  # dplyr::group_split() (not base split()/interaction()) because most
  # non-isotope parameters have species = NA -- interaction()'s default
  # NA handling drops every group entirely when one of the grouping
  # columns is all-NA, silently producing zero groups.
  groups <- measured |>
    dplyr::group_by(.data$core_id, .data$parameter_group, .data$parameter_name, .data$species) |>
    dplyr::group_split()

  out <- purrr::map_dfr(groups, function(g) {
    core_id <- g$core_id[1]
    td <- target_depths_by_core[[core_id]]
    if (is.null(td)) {
      return(tibble::tibble())
    }
    interpolate_one_group(g, td, method, gap_multiplier, snap_tolerance_cm)
  })

  if (nrow(out) == 0) {
    out <- tibble::tibble(
      core_id = character(), pangaea_doi = character(), depth_cm = numeric(),
      parameter_group = character(), parameter_name = character(),
      species = character(), value = numeric(), unit = character(),
      raw_column_name = character(), is_interpolated = logical(),
      interp_gap_cm = numeric(), retrieved_at = as.POSIXct(character()),
      bracket_depth_min_cm = numeric(), bracket_depth_max_cm = numeric()
    )
  }

  new_param_long(out)
}

#' @noRd
cross_depths_params <- function(sample_depths, params) {
  n_d <- length(sample_depths)
  n_p <- nrow(params)
  tibble::tibble(
    depth_cm = rep(sample_depths, times = n_p),
    parameter_group = rep(params$parameter_group, each = n_d),
    parameter_name = rep(params$parameter_name, each = n_d),
    species = rep(params$species, each = n_d),
    unit = rep(params$unit, each = n_d)
  )
}

#' Report depth-match coverage for every sample depth, including gaps
#'
#' [interpolate_depth_series()] only ever returns rows it could actually
#' resolve, so a sample depth with no usable data nearby is simply absent
#' from its output -- indistinguishable from "not asked for". This wraps
#' it so every `(sample depth, parameter_group, parameter_name, species)`
#' combination that exists anywhere in `x` gets a row, even when nothing
#' was found there, flagged via a `match_status` column instead of being
#' silently dropped -- e.g. samples at 40-200 cm against data retrieved
#' at 300-1000 cm should say so, not just return nothing.
#'
#' @param x A `param_long` tibble of measured values -- may span more than
#'   one `core_id` (e.g. a core's own data plus data borrowed from
#'   neighbouring cores); `sample_depths` is matched against all of it.
#' @param sample_depths Numeric vector of depths to report coverage for
#'   (e.g. one core's own sample depths).
#' @param method,gap_multiplier,snap_tolerance_cm Passed to
#'   [interpolate_depth_series()].
#' @return A `param_long` tibble, one row per `(sample depth,
#'   parameter_group, parameter_name, species)` combination present in
#'   `x`, plus:
#'   - `match_status`: `"measured"` (an exact or within-tolerance match --
#'     see `snap_tolerance_cm`), `"interpolated"`, `"no data in that depth
#'     range"` (this parameter's own retrieved depths don't overlap the
#'     sample-depth range at all), or `"gap too large to interpolate"`
#'     (they do overlap, but no bracketing pair was close enough).
#'   - `measured_depth_min`/`measured_depth_max`: this parameter's own
#'     retrieved depth range in `x` (real, non-`NA` values only), for
#'     building diagnostic messages.
#'   - `bracket_depth_min_cm`/`bracket_depth_max_cm`: the real source
#'     depth(s) (in `core_id`, which may be a neighbour core, not the one
#'     `sample_depths` belongs to) this row's value was actually matched
#'     from -- see [interpolate_depth_series()] for the exact-vs-snapped-
#'     vs-interpolated distinction. `NA` on unresolved rows.
#'   `core_id`/`pangaea_doi`/`value`/etc. are `NA` on unresolved rows.
#' @export
report_depth_coverage <- function(x, sample_depths,
                                   method = c("linear", "spline"),
                                   gap_multiplier = 10,
                                   snap_tolerance_cm = 4) {
  validate_param_long(x)
  method <- match.arg(method)

  sample_depths <- unique(sample_depths[!is.na(sample_depths)])
  if (length(sample_depths) == 0) {
    rlang::abort("`sample_depths` must contain at least one non-NA depth.")
  }
  sample_range <- range(sample_depths)

  # NA-valued rows (a depth shared with another parameter on a
  # composite/spliced depth scale, where this parameter itself wasn't
  # measured there) aren't real measurements -- excluded here too, so the
  # reported depth range and the parameter list both reflect actual
  # coverage, not the shared scale's full extent.
  measured <- x[!x$is_interpolated & !is.na(x$value), ]
  if (nrow(measured) == 0) {
    rlang::abort("`x` has no measured (non-interpolated) rows to report coverage for.")
  }

  params <- unique(measured[, c("parameter_group", "parameter_name", "species", "unit")])

  coverage <- measured |>
    dplyr::group_by(.data$parameter_group, .data$parameter_name, .data$species) |>
    dplyr::summarise(
      measured_depth_min = min(.data$depth_cm, na.rm = TRUE),
      measured_depth_max = max(.data$depth_cm, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(in_range = .data$measured_depth_min <= sample_range[2] & .data$measured_depth_max >= sample_range[1])

  core_ids <- unique(x$core_id)
  targets <- tibble::tibble(
    core_id = rep(core_ids, each = length(sample_depths)),
    depth_cm = rep(sample_depths, times = length(core_ids))
  )
  matched <- interpolate_depth_series(
    x, targets,
    method = method, gap_multiplier = gap_multiplier, snap_tolerance_cm = snap_tolerance_cm
  )

  full_grid <- cross_depths_params(sample_depths, params)

  out <- dplyr::left_join(
    full_grid, matched,
    by = c("depth_cm", "parameter_group", "parameter_name", "species", "unit")
  )
  out <- dplyr::left_join(
    out, coverage,
    by = c("parameter_group", "parameter_name", "species")
  )

  out$match_status <- dplyr::case_when(
    !is.na(out$value) & !out$is_interpolated ~ "measured",
    !is.na(out$value) & out$is_interpolated ~ "interpolated",
    !out$in_range ~ "no data in that depth range",
    TRUE ~ "gap too large to interpolate"
  )
  out$is_interpolated[is.na(out$is_interpolated)] <- FALSE

  out <- dplyr::select(
    out,
    "core_id", "pangaea_doi", "depth_cm", "parameter_group", "parameter_name",
    "species", "value", "unit", "raw_column_name", "is_interpolated",
    "interp_gap_cm", "retrieved_at", "match_status",
    "bracket_depth_min_cm", "bracket_depth_max_cm",
    "measured_depth_min", "measured_depth_max"
  )

  new_param_long(out)
}
