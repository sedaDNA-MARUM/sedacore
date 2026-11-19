make_age_series <- function(core_id, depths, values, parameter_name = "TOC") {
  new_param_long(tibble::tibble(
    core_id = core_id, pangaea_doi = "doi1", depth_cm = depths,
    parameter_group = "toc", parameter_name = parameter_name,
    species = NA_character_, value = values, unit = "%",
    raw_column_name = "raw", is_interpolated = FALSE,
    interp_gap_cm = NA_real_, retrieved_at = Sys.time()
  ))
}

test_that("core_age_model extracts one core's own age_model rows, sorted", {
  x <- dplyr::bind_rows(
    make_age_series("core1", c(20, 0, 10), c(20, 0, 10), parameter_name = "age_model"),
    make_age_series("core2", c(0, 10), c(0, 30), parameter_name = "age_model"),
    make_age_series("core1", c(0, 10), c(1, 2), parameter_name = "TOC")
  )
  out <- core_age_model(x, "core1")

  expect_equal(out$depth_cm, c(0, 10, 20))
  expect_equal(out$age_ka, c(0, 10, 20))
})

test_that("core_age_model returns zero rows when a core has no age_model data", {
  x <- make_age_series("core1", c(0, 10), c(1, 2), parameter_name = "TOC")
  out <- core_age_model(x, "core1")
  expect_equal(nrow(out), 0)
})

test_that("core_age_model excludes NA-valued age_model rows", {
  x <- make_age_series("core1", c(0, 5, 10), c(0, NA, 20), parameter_name = "age_model")
  out <- core_age_model(x, "core1")
  expect_equal(out$depth_cm, c(0, 10))
})

test_that("depth_to_age returns NA for every depth when the age model is empty", {
  empty <- core_age_model(make_age_series("core1", c(0, 10), c(1, 2)), "core1")
  out <- depth_to_age(empty, c(0, 10))
  expect_true(all(is.na(out)))
})

test_that("depth_to_age resolves an exact match", {
  am <- tibble::tibble(depth_cm = c(0, 10, 20), age_ka = c(0, 10, 20))
  expect_equal(depth_to_age(am, 10), 10)
})

test_that("depth_to_age snaps to a close control point within snap_tolerance_cm", {
  am <- tibble::tibble(depth_cm = c(0, 20), age_ka = c(0, 20))
  # 2 cm away from the 0 cm control point -- within the default tolerance.
  expect_equal(depth_to_age(am, 2, snap_tolerance_cm = 4), 0)
})

test_that("depth_to_age linearly interpolates between two bracketing control points", {
  am <- tibble::tibble(depth_cm = c(0, 20), age_ka = c(0, 20))
  expect_equal(depth_to_age(am, 10, snap_tolerance_cm = 0), 10)
})

test_that("depth_to_age never extrapolates beyond the age model's own range", {
  am <- tibble::tibble(depth_cm = c(10, 20), age_ka = c(10, 20))
  out <- depth_to_age(am, c(0, 30))
  expect_true(all(is.na(out)))
})

test_that("depth_to_age respects the gap-size guard between bracketing control points", {
  am <- tibble::tibble(depth_cm = c(0, 10, 200), age_ka = c(0, 10, 200))
  # spacing (median of 10, 190) heavily skewed by the one big gap already
  # captures it -- use a tight gap_multiplier to force a real rejection.
  out <- depth_to_age(am, 105, gap_multiplier = 1, snap_tolerance_cm = 0)
  expect_true(is.na(out))
})

test_that("depth_to_age is vectorized and independent per depth", {
  am <- tibble::tibble(depth_cm = c(0, 10, 20), age_ka = c(0, 10, 20))
  out <- depth_to_age(am, c(0, 5, 20, 100), snap_tolerance_cm = 0)
  expect_equal(out, c(0, 5, 20, NA_real_))
})

test_that("depth_to_age returns NA for an NA depth instead of erroring (age-only series)", {
  am <- tibble::tibble(depth_cm = c(0, 10, 20), age_ka = c(0, 10, 20))
  expect_equal(depth_to_age(am, c(NA, 5), snap_tolerance_cm = 0), c(NA_real_, 5))
})
