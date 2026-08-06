example_path <- system.file(
  "extdata", "example_core_metadata.tsv",
  package = "sedacore"
)

raw <- if (nzchar(example_path)) {
  readr::read_tsv(example_path, show_col_types = FALSE)
} else {
  readr::read_tsv(
    testthat::test_path("..", "..", "inst", "extdata", "example_core_metadata.tsv"),
    show_col_types = FALSE
  )
}

test_that("as_sample_registry maps raw columns", {
  sr <- as_sample_registry(raw)
  expect_true(is_sample_registry(sr))
  expect_equal(nrow(sr), nrow(raw))
  expect_true("core_label" %in% names(sr))
})

test_that("as_sample_registry errors on missing required columns", {
  expect_error(as_sample_registry(raw, core_id = "nope"))
  expect_error(as_sample_registry(raw, depth_cm = "nope"))
})

test_that("as_sample_registry fills optional columns sensibly when absent", {
  sr <- as_sample_registry(raw, water_depth_m = "nope", sample_id = "nope", age_ka = "nope")
  expect_true(all(is.na(sr$water_depth_m)))
  expect_true(all(is.na(sr$age_ka)))
  expect_false(any(is.na(sr$sample_id)))
  expect_equal(sr$sample_id, paste(sr$core_id, seq_len(nrow(sr)), sep = "_"))
})

test_that("collapse_to_core_registry produces one row per physical core", {
  sr <- as_sample_registry(raw)
  cr <- collapse_to_core_registry(sr)
  expect_true(is_core_registry(cr))
  expect_equal(nrow(cr), dplyr::n_distinct(raw$sample_source))
  expect_setequal(cr$core_id, unique(raw$sample_source))
})

test_that("collapse_to_core_registry classifies source_program", {
  sr <- as_sample_registry(raw)
  cr <- collapse_to_core_registry(sr)
  expect_equal(cr$source_program[cr$core_id == "U1499A"], "IODP")
  expect_true(all(cr$source_program[startsWith(cr$core_id, "EXCORE")] == "other"))
})

test_that("plot_core_map returns a ggplot", {
  cr <- collapse_to_core_registry(as_sample_registry(raw))
  p <- plot_core_map(cr)
  expect_s3_class(p, "ggplot")
})

test_that("leaflet_core_map returns a leaflet widget", {
  cr <- collapse_to_core_registry(as_sample_registry(raw))
  m <- leaflet_core_map(cr)
  expect_s3_class(m, "leaflet")
})

test_that("leaflet_core_map uses distinct, non-blue colors for input vs neighbour cores", {
  expect_false(grepl("^#0", map_colors[["input"]], ignore.case = TRUE))
  expect_false(identical(map_colors[["input"]], map_colors[["neighbour"]]))
})

test_that("leaflet_core_map adds a second layer + legend when neighbours are supplied", {
  cr <- collapse_to_core_registry(as_sample_registry(raw))
  nb <- new_neighbour_registry(tibble::tibble(
    core_id = cr$core_id[1], neighbour_doi = "10.1594/PANGAEA.999999",
    neighbour_core_id = "SOME-OTHER-CORE",
    latitude = cr$latitude[1] + 0.1, longitude = cr$longitude[1] + 0.1,
    distance_km = 12, source_program = "other", date_time = as.POSIXct(NA),
    citation = "Some citation", available_parameters = NA_character_,
    included_in_export = TRUE
  ))

  m <- leaflet_core_map(cr, neighbours = nb)
  calls <- vapply(m$x$calls, function(c) c$method, character(1))

  expect_true(sum(calls == "addCircleMarkers") == 2)
  expect_true("addLegend" %in% calls)
})

test_that("leaflet_core_map drops neighbour rows with missing coordinates without erroring", {
  cr <- collapse_to_core_registry(as_sample_registry(raw))
  nb <- new_neighbour_registry(tibble::tibble(
    core_id = cr$core_id[1], neighbour_doi = "10.1594/PANGAEA.999999",
    neighbour_core_id = "NO-COORDS",
    latitude = NA_real_, longitude = NA_real_,
    distance_km = NA_real_, source_program = "other", date_time = as.POSIXct(NA),
    citation = "Some citation", available_parameters = NA_character_,
    included_in_export = TRUE
  ))

  m <- leaflet_core_map(cr, neighbours = nb)
  calls <- vapply(m$x$calls, function(c) c$method, character(1))
  expect_true(sum(calls == "addCircleMarkers") == 1)
})
