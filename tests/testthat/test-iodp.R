test_that("classify_core_source recognizes IODP/ODP/DSDP/GeoB naming conventions", {
  expect_equal(classify_core_source("GeoB1706-2"), "GeoB")
  expect_equal(classify_core_source("U1478A"), "IODP")
  expect_equal(classify_core_source("361-U1478"), "IODP")
  expect_equal(classify_core_source("108-659A"), "ODP")
  expect_equal(classify_core_source("DSDP 502"), "DSDP")
  expect_equal(classify_core_source("some-random-core"), "other")
  expect_true(is.na(classify_core_source(NA_character_)))
})

test_that("search_iodp_campaigns requires expedition or bbox", {
  expect_error(search_iodp_campaigns())
})

test_that("search_iodp_campaigns finds real IODP Expedition 361 datasets", {
  vcr::use_cassette("search_iodp_campaigns_361", {
    res <- search_iodp_campaigns(expedition = "361", count = 5)
  })

  expect_true(nrow(res) > 0)
  expect_true(all(c("score", "doi", "citation") %in% names(res)))
})

test_that("search_iodp_campaigns results feed straight into pg_fetch_dataset/harmonize_columns", {
  vcr::use_cassette("iodp_361_end_to_end", {
    d <- pg_fetch_dataset("10.1594/PANGAEA.923959")
  })

  expect_true(all(grepl("^361-U", d$event_label)))
  expect_equal(classify_core_source(d$event_label[1]), "IODP")
})
