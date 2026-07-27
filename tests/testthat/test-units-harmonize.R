base_param_long <- function(...) {
  x <- tibble::tibble(
    core_id = "core1", pangaea_doi = "doi1", depth_cm = 10,
    parameter_group = "susceptibility", parameter_name = "kappa",
    species = NA_character_, value = 5, unit = "10**-6 SI",
    raw_column_name = "kappa [10**-6 SI]", is_interpolated = FALSE,
    interp_gap_cm = NA_real_, retrieved_at = Sys.time()
  )
  overrides <- list(...)
  for (nm in names(overrides)) x[[nm]] <- overrides[[nm]]
  new_param_long(x)
}

test_that("standardize_units is identity for units already at target", {
  x <- base_param_long()
  out <- standardize_units(x)
  expect_equal(out$value, 5)
  expect_equal(out$unit, "10**-6 SI")
})

test_that("standardize_units rescales a non-target unit", {
  conv <- tibble::tibble(
    parameter_name = "kappa", unit_raw = "10**-5 SI",
    unit_target = "10**-6 SI", multiplier = 10, offset = 0
  )
  x <- base_param_long(unit = "10**-5 SI")
  out <- standardize_units(x, conversions = conv)
  expect_equal(out$value, 50)
  expect_equal(out$unit, "10**-6 SI")
})

test_that("standardize_units passes through unmapped combos by default", {
  x <- base_param_long(parameter_name = "unknown_param")
  expect_silent(out <- standardize_units(x, unmapped = "pass_through"))
  expect_equal(out$value, x$value)
})

test_that("standardize_units can warn or error on unmapped combos", {
  x <- base_param_long(parameter_name = "unknown_param")
  expect_warning(standardize_units(x, unmapped = "warn"), class = "sedacore_unmapped_units")
  expect_error(standardize_units(x, unmapped = "error"), class = "sedacore_unmapped_units")
})
