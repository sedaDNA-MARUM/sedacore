test_that("pg_search_cores returns candidates per core", {
  cr <- new_core_registry(tibble::tibble(
    core_id = "GeoB1706-2", core_label = "GeoB1706-2; 2003",
    latitude = -19.56, longitude = 11.18, water_depth_m = 980,
    source_program = "GeoB", n_samples = 18L,
    depth_min_cm = 5, depth_max_cm = 160, pangaea_doi = NA_character_
  ))

  vcr::use_cassette("pg_search_cores_geob1706", {
    res <- pg_search_cores(cr, count = 5)
  })

  expect_true(all(c("core_id", "score", "doi", "citation") %in% names(res)))
  expect_true(all(res$core_id == "GeoB1706-2"))
  expect_true(nrow(res) > 0)
  expect_true(any(grepl("^10\\.1594/PANGAEA", res$doi)))
})
