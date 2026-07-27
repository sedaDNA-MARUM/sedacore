test_that("new_sample_registry validates required columns", {
  df <- tibble::tibble(
    core_id = "c1", depth_cm = 1, latitude = 1, longitude = 1,
    water_depth_m = 1, sample_id = "a", age_ka = 1
  )
  x <- new_sample_registry(df)
  expect_true(is_sample_registry(x))
  expect_silent(validate_sample_registry(x))
})

test_that("new_sample_registry errors on missing columns", {
  df <- tibble::tibble(core_id = "c1")
  expect_error(new_sample_registry(df), class = "sedacore_schema_error")
})

test_that("validate_* rejects wrong class", {
  df <- tibble::tibble(
    core_id = "c1", core_label = "c1", latitude = 1, longitude = 1,
    water_depth_m = 1, source_program = "GeoB", n_samples = 1L,
    depth_min_cm = 1, depth_max_cm = 1, pangaea_doi = NA_character_
  )
  expect_error(validate_sample_registry(new_core_registry(df)))
})

test_that("param_long round-trips", {
  df <- tibble::tibble(
    core_id = "c1", pangaea_doi = "doi", depth_cm = 1,
    parameter_group = "isotopes", parameter_name = "d18O", species = "sp",
    value = 1.5, unit = "permil", raw_column_name = "raw",
    is_interpolated = FALSE, interp_gap_cm = NA_real_,
    retrieved_at = Sys.time()
  )
  x <- new_param_long(df)
  expect_true(is_param_long(x))
  expect_silent(validate_param_long(x))
})
