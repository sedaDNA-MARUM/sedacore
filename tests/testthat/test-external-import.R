test_that("import_external_data reads the shipped template into param_long", {
  template <- system.file("extdata", "external_import_template.csv", package = "sedacore")
  out <- import_external_data(template, core_id = "U1478A")

  expect_true(is_param_long(out))
  expect_equal(nrow(out), 2)
  expect_true(all(out$core_id == "U1478A"))
  expect_setequal(out$parameter_name, c("d18O", "d13C"))
})

test_that("import_external_data errors on a malformed CSV", {
  bad <- tempfile(fileext = ".csv")
  readr::write_csv(tibble::tibble(depth_cm = 1, value = 2), bad)
  expect_error(import_external_data(bad, core_id = "U1478A"))
})

test_that("import_external_data records source and doi as passthrough provenance", {
  template <- system.file("extdata", "external_import_template.csv", package = "sedacore")
  out <- import_external_data(
    template, core_id = "U1478A",
    source = "Smith et al. 2020, Table 2", doi = "10.1234/example"
  )

  expect_true(all(out$source == "Smith et al. 2020, Table 2"))
  expect_true(all(out$pangaea_doi == "10.1234/example"))
})

test_that("import_external_data defaults source/doi to NA when not given", {
  template <- system.file("extdata", "external_import_template.csv", package = "sedacore")
  out <- import_external_data(template, core_id = "U1478A")

  expect_true(all(is.na(out$source)))
  expect_true(all(is.na(out$pangaea_doi)))
})
