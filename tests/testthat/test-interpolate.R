make_series <- function(depths, values, core_id = "core1", parameter_name = "d18O",
                         species = "G. inflata") {
  new_param_long(tibble::tibble(
    core_id = core_id, pangaea_doi = "doi1", depth_cm = depths,
    parameter_group = "isotopes", parameter_name = parameter_name,
    species = species, value = values, unit = "permil PDB",
    raw_column_name = "raw", is_interpolated = FALSE,
    interp_gap_cm = NA_real_, retrieved_at = Sys.time()
  ))
}

test_that("interpolate_depth_series linearly interpolates a mid-gap depth", {
  x <- make_series(depths = c(0, 10, 20), values = c(0, 1, 2))
  targets <- tibble::tibble(core_id = "core1", depth_cm = 5)

  out <- interpolate_depth_series(x, targets, method = "linear", gap_multiplier = 10)

  expect_equal(nrow(out), 1)
  expect_equal(out$value, 0.5)
  expect_true(out$is_interpolated)
  expect_equal(out$interp_gap_cm, 10)
  # The real bracketing depths (0 and 10) this estimate was computed
  # between, distinct from depth_cm (5, the target/sample depth).
  expect_equal(out$bracket_depth_min_cm, 0)
  expect_equal(out$bracket_depth_max_cm, 10)
})

test_that("interpolate_depth_series does not extrapolate beyond the measured range", {
  x <- make_series(depths = c(0, 10, 20), values = c(0, 1, 2))
  targets <- tibble::tibble(core_id = "core1", depth_cm = c(-5, 25))

  out <- interpolate_depth_series(x, targets)
  expect_equal(nrow(out), 0)
})

test_that("interpolate_depth_series respects the gap-size guard", {
  # depths 0, 10, 20 (spacing 10), then a big desert to 200: target at 110
  # sits between 20 and 200 -- a gap of 180, way beyond 10x the median
  # spacing of 10, so it should NOT be filled.
  x <- make_series(depths = c(0, 10, 20, 200), values = c(0, 1, 2, 20))
  targets <- tibble::tibble(core_id = "core1", depth_cm = 110)

  out <- interpolate_depth_series(x, targets, gap_multiplier = 10)
  expect_equal(nrow(out), 0)
})

test_that("interpolate_depth_series returns already-measured target depths as exact matches", {
  x <- make_series(depths = c(0, 10, 20), values = c(0, 1, 2))
  targets <- tibble::tibble(core_id = "core1", depth_cm = 10)

  out <- interpolate_depth_series(x, targets)
  expect_equal(nrow(out), 1)
  expect_equal(out$value, 1)
  expect_false(out$is_interpolated)
  expect_true(is.na(out$interp_gap_cm))
  # An exact match isn't bracketed by two different depths -- both
  # bounds collapse to the real (== target) depth itself.
  expect_equal(out$bracket_depth_min_cm, 10)
  expect_equal(out$bracket_depth_max_cm, 10)
})

test_that("interpolate_depth_series keeps groups (parameter/species) separate", {
  x <- dplyr::bind_rows(
    make_series(c(0, 20), c(0, 2), parameter_name = "d18O"),
    make_series(c(0, 20), c(10, 30), parameter_name = "d13C")
  )
  targets <- tibble::tibble(core_id = "core1", depth_cm = 10)

  out <- interpolate_depth_series(x, targets)
  expect_equal(nrow(out), 2)
  expect_setequal(out$parameter_name, c("d18O", "d13C"))
  expect_equal(out$value[out$parameter_name == "d18O"], 1)
  expect_equal(out$value[out$parameter_name == "d13C"], 20)
})

test_that("interpolate_depth_series rounds interpolated values to clear floating-point noise", {
  # stats::approx()'s own arithmetic can produce noise like
  # 16.740000000000002 for otherwise-clean 1-decimal input values --
  # exposed the same way as the depth_cm rounding fix, via
  # readr::write_tsv()'s full-precision serialization.
  x <- make_series(depths = c(146, 160), values = c(16.6, 16.9), parameter_name = "SST")
  targets <- tibble::tibble(core_id = "core1", depth_cm = 150)

  # snap_tolerance_cm = 0: this test is about interpolation rounding
  # specifically, not the (closer, 4 cm-away) snap-to-nearest behavior.
  out <- interpolate_depth_series(x, targets, gap_multiplier = 10, snap_tolerance_cm = 0)
  raw_approx <- stats::approx(c(146, 160), c(16.6, 16.9), xout = 150)$y
  expect_true(raw_approx != round(raw_approx, 6)) # confirms the raw arithmetic does produce noise
  expect_equal(out$value, round(raw_approx, 6))
})

test_that("interpolate_depth_series works when species is NA (most non-isotope parameters)", {
  # Regression test: base split()/interaction() drops every group when one
  # of the grouping columns (species) is entirely NA, silently producing
  # zero rows for every non-isotope parameter (biomarkers, toc, porewater
  # chemistry, ...) -- found via live testing on a real alkenone dataset.
  x <- make_series(depths = c(0, 10, 20), values = c(0, 1, 2),
    parameter_name = "K37", species = NA_character_
  )
  targets <- tibble::tibble(core_id = "core1", depth_cm = 5)

  out <- interpolate_depth_series(x, targets, method = "linear", gap_multiplier = 10)

  expect_equal(nrow(out), 1)
  expect_equal(out$value, 0.5)
  expect_true(out$is_interpolated)
  expect_true(is.na(out$species))
})

test_that("report_depth_coverage returns one row per sample depth, flagging measured/interpolated", {
  x <- make_series(depths = c(0, 10, 20), values = c(0, 1, 2))
  out <- report_depth_coverage(x, sample_depths = c(0, 5, 10))

  expect_equal(nrow(out), 3)
  expect_equal(out$match_status[out$depth_cm == 0], "measured")
  expect_equal(out$match_status[out$depth_cm == 5], "interpolated")
  expect_equal(out$match_status[out$depth_cm == 10], "measured")
  expect_equal(out$value[out$depth_cm == 5], 0.5)
})

test_that("report_depth_coverage flags sample depths outside the retrieved depth range", {
  # The exact scenario reported by a user: samples at 40-200 cm, but the
  # only data retrieved for this parameter is at 300-1000 cm -- previously
  # this just silently returned zero rows with no explanation.
  x <- make_series(depths = c(300, 1000), values = c(1, 2))
  out <- report_depth_coverage(x, sample_depths = c(40, 200))

  expect_equal(nrow(out), 2)
  expect_true(all(out$match_status == "no data in that depth range"))
  expect_true(all(is.na(out$value)))
  expect_true(all(is.na(out$bracket_depth_min_cm)))
  expect_true(all(is.na(out$bracket_depth_max_cm)))
  expect_equal(unique(out$measured_depth_min), 300)
  expect_equal(unique(out$measured_depth_max), 1000)
})

test_that("report_depth_coverage flags a sample depth within range but across too large a gap", {
  x <- make_series(depths = c(0, 10, 20, 200), values = c(0, 1, 2, 20))
  out <- report_depth_coverage(x, sample_depths = 110, gap_multiplier = 10)

  expect_equal(nrow(out), 1)
  expect_equal(out$match_status, "gap too large to interpolate")
  expect_true(is.na(out$value))
})

test_that("report_depth_coverage matches data borrowed from a different core_id (neighbour)", {
  # e.g. GeoB18549-1 has no PANGAEA data of its own; its neighbour
  # GeoB1706-2's harmonized data is passed in instead, matched against
  # GeoB18549-1's own sample depths.
  x <- make_series(depths = c(0, 10, 20), values = c(0, 1, 2), core_id = "GeoB1706-2")
  out <- report_depth_coverage(x, sample_depths = c(0, 5))

  expect_equal(nrow(out), 2)
  expect_equal(out$core_id[out$depth_cm == 0], "GeoB1706-2")
  expect_equal(out$match_status[out$depth_cm == 0], "measured")
  expect_equal(out$match_status[out$depth_cm == 5], "interpolated")
  # The real depths this came from are the neighbour's own (0 and 10,
  # both in core_id "GeoB1706-2"), distinct from depth_cm (the target
  # core's own sample depths, 0 and 5).
  expect_equal(out$bracket_depth_min_cm[out$depth_cm == 0], 0)
  expect_equal(out$bracket_depth_max_cm[out$depth_cm == 0], 0)
  expect_equal(out$bracket_depth_min_cm[out$depth_cm == 5], 0)
  expect_equal(out$bracket_depth_max_cm[out$depth_cm == 5], 10)
})

test_that("report_depth_coverage keeps parameters separate when only one has coverage", {
  x <- dplyr::bind_rows(
    make_series(c(0, 10), c(0, 1), parameter_name = "d18O"),
    make_series(c(300, 1000), c(9, 10), parameter_name = "TOC")
  )
  out <- report_depth_coverage(x, sample_depths = 5)

  expect_equal(nrow(out), 2)
  expect_equal(out$match_status[out$parameter_name == "d18O"], "interpolated")
  expect_equal(out$match_status[out$parameter_name == "TOC"], "no data in that depth range")
})

test_that("interpolate_depth_series ignores NA-valued rows for exact matching", {
  # Regression test: a real PANGAEA dataset (MD95-2024, TOC) shares one
  # composite depth column across many parameters at fine resolution, but
  # TOC itself was only measured every ~5 cm -- so a target depth (246)
  # can land exactly on a row where *another* parameter was measured but
  # TOC's own cell is blank (NA). That row must not count as an exact
  # match, or it (a) fabricates a bogus "measured" row with value = NA and
  # (b) blocks the real bracket (245.5 -> 250.5) from ever being tried.
  x <- make_series(
    depths = c(240, 245.5, 246, 250.5, 260),
    values = c(0.36, 0.31, NA, 0.30, 0.29),
    parameter_name = "TOC", species = NA_character_
  )
  targets <- tibble::tibble(core_id = "core1", depth_cm = 246)

  out <- interpolate_depth_series(x, targets, gap_multiplier = 10, snap_tolerance_cm = 0)

  expect_equal(nrow(out), 1)
  expect_true(out$is_interpolated)
  expect_equal(out$bracket_depth_min_cm, 245.5)
  expect_equal(out$bracket_depth_max_cm, 250.5)
  expect_equal(out$value, round(stats::approx(c(245.5, 250.5), c(0.31, 0.30), xout = 246)$y, 6))
})

test_that("interpolate_depth_series snaps to a close real measurement within snap_tolerance_cm", {
  x <- make_series(depths = c(240, 245.5, 250.5, 260), values = c(0.36, 0.31, 0.30, 0.29), parameter_name = "TOC")
  targets <- tibble::tibble(core_id = "core1", depth_cm = 246)

  out <- interpolate_depth_series(x, targets, gap_multiplier = 10, snap_tolerance_cm = 4)

  expect_equal(nrow(out), 1)
  expect_false(out$is_interpolated)
  expect_equal(out$value, 0.31)
  expect_equal(out$bracket_depth_min_cm, 245.5)
  expect_equal(out$bracket_depth_max_cm, 245.5)
})

test_that("interpolate_depth_series does not snap beyond snap_tolerance_cm", {
  x <- make_series(depths = c(240, 245.5, 250.5, 260), values = c(0.36, 0.31, 0.30, 0.29), parameter_name = "TOC")
  targets <- tibble::tibble(core_id = "core1", depth_cm = 246)

  out <- interpolate_depth_series(x, targets, gap_multiplier = 10, snap_tolerance_cm = 0.4)

  expect_equal(nrow(out), 1)
  expect_true(out$is_interpolated)
  expect_equal(out$bracket_depth_min_cm, 245.5)
  expect_equal(out$bracket_depth_max_cm, 250.5)
})

test_that("report_depth_coverage marks a snapped-within-tolerance depth as measured", {
  x <- make_series(depths = c(240, 245.5, 250.5, 260), values = c(0.36, 0.31, 0.30, 0.29), parameter_name = "TOC")
  out <- report_depth_coverage(x, sample_depths = 246, snap_tolerance_cm = 4)

  expect_equal(out$match_status, "measured")
  expect_false(out$is_interpolated)
  expect_equal(out$value, 0.31)
})

test_that("find_neighbours_by_depth snaps to the nearest measurement within tolerance", {
  x <- make_series(depths = c(0, 10, 20), values = c(0, 1, 2))
  targets <- tibble::tibble(core_id = "core1", depth_cm = c(9, 50))

  out <- find_neighbours_by_depth(x, targets, tolerance_cm = 3)

  expect_equal(nrow(out), 1)
  expect_equal(out$depth_cm, 10)
  expect_equal(out$sample_depth_cm, 9)
  expect_equal(out$gap_cm, 1)
})

test_that("find_neighbours_by_depth ignores NA-valued rows", {
  # Same composite-depth-scale issue as interpolate_depth_series(): a row
  # closer to the target but with value = NA must not be snapped to over
  # a real (non-NA) value that's further away.
  x <- make_series(depths = c(0, 9, 10, 20), values = c(0, NA, 1, 2))
  targets <- tibble::tibble(core_id = "core1", depth_cm = 9)

  out <- find_neighbours_by_depth(x, targets, tolerance_cm = 3)

  expect_equal(nrow(out), 1)
  expect_equal(out$depth_cm, 10)
  expect_equal(out$value, 1)
})
