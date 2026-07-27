vcr::use_cassette("pg_fetch_dataset_single_event", {
  single_event <- pg_fetch_dataset("10.1594/PANGAEA.218244")
})

test_that("pg_fetch_dataset parses a single-event dataset", {
  d <- single_event

  expect_equal(nrow(d), 1)
  expect_equal(d$doi, "10.1594/PANGAEA.218244")
  expect_equal(round(d$latitude, 2), -19.56)
  expect_equal(round(d$longitude, 2), 11.18)
  expect_equal(d$water_depth_m, 980)
  expect_equal(d$campaign, "M20/2")
  expect_true(is.data.frame(d$data[[1]]))
  expect_true("Depth sed [m]" %in% names(d$data[[1]]))
})

test_that("harmonize_columns end-to-end on a real fetched isotope dataset", {
  d <- single_event
  h <- harmonize_columns(d$data[[1]], core_id = "GeoB1706-2", doi = d$doi[1])

  expect_true(is_param_long(h))
  expect_true(nrow(h) > 0)
  expect_equal(unique(h$parameter_group), "isotopes")
  expect_setequal(unique(h$parameter_name), c("d18O", "d13C"))
  expect_equal(unique(h$species), "G. inflata")
  expect_length(attr(h, "fallback_columns"), 0)
  expect_length(attr(h, "dropped_columns"), 0)
})

test_that("pg_fetch_dataset parses a parent/collection DOI into multiple rows", {
  vcr::use_cassette("pg_fetch_dataset_collection", {
    d <- pg_fetch_dataset("10.1594/PANGAEA.680850")
  })

  expect_gt(nrow(d), 1)
  expect_true(all(!is.na(d$doi)))
  expect_true(all(c("latitude", "longitude", "campaign") %in% names(d)))
})
