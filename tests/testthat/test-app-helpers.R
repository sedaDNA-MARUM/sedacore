mk_pl <- function(core, name, depths, vals, group = "toc", unit = "u", interp = FALSE) {
  new_param_long(tibble::tibble(
    core_id = core, pangaea_doi = "d", depth_cm = depths, parameter_group = group,
    parameter_name = name, species = NA_character_, value = vals, unit = unit,
    raw_column_name = "r", is_interpolated = interp, interp_gap_cm = NA_real_,
    retrieved_at = Sys.time()
  ))
}

# -- empty tables ---------------------------------------------------------------

test_that("empty constructors are zero-row and valid against the schema", {
  expect_equal(nrow(empty_param_long()), 0)
  expect_no_error(validate_param_long(empty_param_long()))
  expect_equal(names(empty_param_long()), names(schema_spec$param_long))

  nb <- empty_neighbour_registry()
  expect_equal(nrow(nb), 0)
  expect_true(inherits(nb, "neighbour_registry"))
  expect_equal(names(nb), names(schema_spec$neighbour_registry))
})

test_that("empty_depth_matched_long adds target_core_id and nothing else", {
  x <- empty_depth_matched_long()
  expect_equal(names(x), c(names(schema_spec$param_long), "target_core_id"))
  expect_type(x$target_core_id, "character")
})

test_that("empty_from_spec rejects a schema type it has no empty vector for", {
  expect_error(empty_from_spec(list(a = "Date")), "No empty vector")
})

test_that("empty tables bind cleanly with real rows (how the app uses them)", {
  out <- dplyr::bind_rows(empty_param_long(), mk_pl("c1", "TOC", 1, 2))
  expect_equal(nrow(out), 1)
})

# -- keep_selected_measurements ---------------------------------------------------

test_that("keep_selected_measurements keeps everything when nothing is selected", {
  h <- dplyr::bind_rows(mk_pl("c", "TOC", 1, 1), mk_pl("c", "TN", 1, 1))
  expect_equal(nrow(keep_selected_measurements(h, NULL)), 2)
  expect_equal(nrow(keep_selected_measurements(h, character())), 2)
})

test_that("keep_selected_measurements always keeps age_model, filters the rest", {
  h <- dplyr::bind_rows(
    mk_pl("c", "TOC", 1, 1), mk_pl("c", "TN", 1, 1), mk_pl("c", "age_model", 1, 5, "age")
  )
  out <- keep_selected_measurements(h, "TOC")
  expect_setequal(out$parameter_name, c("TOC", "age_model"))
})

# -- harmonize_for_store ------------------------------------------------------------

raw_tbl <- function(...) {
  tibble::tibble(`Depth sed [m]` = c(0.1, 0.2), ...)
}

test_that("harmonize_for_store returns rows plus a 'Harmonized N rows' message", {
  res <- suppressWarnings(harmonize_for_store(
    raw_tbl(`TOC [%]` = c(1, 2)), "c1", "doi1", "Depth sed [m]", load_synonym_map()
  ))
  expect_equal(nrow(res$data), 2)
  expect_equal(res$messages, "Harmonized 2 rows from doi1.")
})

test_that("harmonize_for_store applies the measurement filter after harmonizing", {
  res <- suppressWarnings(harmonize_for_store(
    raw_tbl(`TOC [%]` = c(1, 2), `TN [%]` = c(0.1, 0.2), `Age model [ka]` = c(1, 2)),
    "c1", "doi1", "Depth sed [m]", load_synonym_map(),
    measurement_filter = "TOC"
  ))
  # TOC (selected) and age_model (always kept) survive; TN does not.
  expect_setequal(res$data$parameter_name, c("TOC", "age_model"))
})

test_that("harmonize_for_store reports, rather than appends, when the filter leaves nothing", {
  res <- suppressWarnings(harmonize_for_store(
    raw_tbl(`TOC [%]` = c(1, 2)), "c1", "doi1", "Depth sed [m]", load_synonym_map(),
    measurement_filter = "TN"
  ))
  expect_null(res$data)
  expect_match(res$messages[1], "0 rows from doi1 after restricting")
  expect_match(res$messages[1], "2 rows available before filtering")
})

test_that("harmonize_for_store lists fallback and dropped columns, naming the synonym csv", {
  res <- suppressWarnings(harmonize_for_store(
    raw_tbl(`TOC [%]` = c(1, 2), `Mystery [x]` = c(5, 6), `Comment` = c("a", "b")),
    "c1", "doi1", "Depth sed [m]", load_synonym_map(),
    synonym_csv_name = "my_synonyms.csv"
  ))
  expect_match(res$messages[2], "Not in my_synonyms.csv, imported as-is: Mystery \\[x\\]")
  expect_match(res$messages[3], "Excluded \\(non-numeric.*Comment")
})

test_that("harmonize_for_store turns a harmonize error into a log message, not a crash", {
  res <- harmonize_for_store(
    tibble::tibble(`TOC [%]` = 1), "c1", "doi1", "Depth sed [m]", load_synonym_map()
  )
  expect_null(res$data)
  expect_match(res$messages, "^ERROR harmonizing:")
})

# -- attach_age_ka -----------------------------------------------------------------

test_that("attach_age_ka hides age_model rows and fills age_ka per core, after depth_cm", {
  x <- dplyr::bind_rows(
    mk_pl("c1", "age_model", c(0, 100), c(0, 10), "age"),
    mk_pl("c1", "TOC", c(0, 50, 100), c(1, 2, 3)),
    mk_pl("c2", "age_model", c(0, 100), c(0, 40), "age"),
    mk_pl("c2", "TOC", 50, 9),
    mk_pl("c3", "TOC", 50, 9) # no age model of its own
  )
  out <- attach_age_ka(x)

  expect_false("age_model" %in% out$parameter_name)
  expect_equal(which(names(out) == "age_ka"), which(names(out) == "depth_cm") + 1)
  expect_equal(out$age_ka[out$core_id == "c1"], c(0, 5, 10))
  expect_equal(out$age_ka[out$core_id == "c2"], 20) # c2's own scale, not c1's
  expect_true(is.na(out$age_ka[out$core_id == "c3"]))
})

test_that("attach_age_ka works with no age models at all, and on an empty table", {
  out <- attach_age_ka(mk_pl("c1", "TOC", 1, 1))
  expect_true(is.na(out$age_ka))

  empty <- attach_age_ka(empty_param_long())
  expect_equal(nrow(empty), 0)
  expect_true("age_ka" %in% names(empty))
})

test_that("attach_age_ka keeps an existing age_ka and only fills the gaps", {
  x <- dplyr::bind_rows(
    mk_pl("c1", "age_model", c(0, 100), c(0, 10), "age"),
    mk_pl("c1", "TOC", c(50, 100), c(1, 2))
  )
  x$age_ka <- c(NA, NA, 99, NA) # a pre-existing native age on the first TOC row
  out <- attach_age_ka(x)
  expect_equal(out$age_ka, c(99, 10))
})

test_that("attach_age_ka tolerates NA depths (age-only series)", {
  x <- dplyr::bind_rows(
    mk_pl("c1", "age_model", c(0, 100), c(0, 10), "age"),
    mk_pl("c1", "TOC", c(NA, 50), c(1, 2))
  )
  out <- attach_age_ka(x)
  expect_equal(out$age_ka, c(NA, 5))
})

# -- summarise_measurements -----------------------------------------------------------

test_that("summarise_measurements counts cores and measured vs interpolated rows", {
  x <- dplyr::bind_rows(
    mk_pl("a", "TOC", c(1, 2), c(1, 1)),
    mk_pl("b", "TOC", 1, 1),
    mk_pl("b", "TOC", 3, 1, interp = TRUE)
  )
  out <- summarise_measurements(x)
  expect_equal(nrow(out), 1)
  expect_equal(out$n_cores, 2)
  expect_equal(out$cores, "a, b")
  expect_equal(out$n_measured, 3)
  expect_equal(out$n_interpolated, 1)
})

test_that("summarise_measurements truncates a long core list", {
  x <- dplyr::bind_rows(lapply(sprintf("core%02d", 1:10), function(id) mk_pl(id, "TOC", 1, 1)))
  out <- summarise_measurements(x, max_cores = 3)
  expect_equal(out$n_cores, 10)
  expect_equal(out$cores, "core01, core02, core03 ... (+7 more)")
})

test_that("summarise_measurements handles an empty table", {
  expect_equal(nrow(summarise_measurements(empty_param_long())), 0)
})

# -- interp_candidate_data --------------------------------------------------------------

test_that("interp_candidate_data selects the active core plus its own neighbours only", {
  pl <- dplyr::bind_rows(
    mk_pl("T", "TOC", 1, 1), mk_pl("N1", "TOC", 1, 1), mk_pl("OTHER", "TOC", 1, 1)
  )
  nb <- new_neighbour_registry(tibble::tibble(
    core_id = c("T", "SOMEONE_ELSE"), neighbour_doi = "x", neighbour_core_id = c("N1", "OTHER"),
    latitude = 0, longitude = 0, distance_km = 1, source_program = "other",
    date_time = Sys.time(), citation = "c", available_parameters = "TOC", included_in_export = TRUE
  ))
  res <- interp_candidate_data(pl, nb, "T")
  expect_setequal(res$data$core_id, c("T", "N1"))
  expect_equal(res$neighbour_ids, "N1")
  expect_null(res$message)
})

test_that("interp_candidate_data explains an empty result", {
  res <- interp_candidate_data(empty_param_long(), empty_neighbour_registry(), "T")
  expect_equal(res$message, "No harmonized data yet for the active core; fetch some first.")

  pl <- mk_pl("T", "TOC", 1, 1)
  nb <- new_neighbour_registry(tibble::tibble(
    core_id = "T", neighbour_doi = "x", neighbour_core_id = "N1", latitude = 0, longitude = 0,
    distance_km = 1, source_program = "other", date_time = Sys.time(), citation = "c",
    available_parameters = "TOC", included_in_export = TRUE
  ))
  res2 <- interp_candidate_data(pl, nb, "T", measurement_filter = "TN")
  expect_match(res2$message, "or its fetched neighbours matching the selected measurement filter; fetch some first.")
})

# -- depth_coverage_messages ---------------------------------------------------------------

sr <- function(ages) {
  new_sample_registry(tibble::tibble(
    core_id = "T", depth_cm = seq(50, by = 10, length.out = length(ages)), latitude = 1, longitude = 1,
    water_depth_m = 1, sample_id = paste0("s", seq_along(ages)), age_ka = ages
  ))
}

test_that("depth_coverage_messages tallies a plain same-core report", {
  x <- mk_pl("T", "TOC", c(0, 10, 20), c(1, 2, 3))
  report <- report_depth_coverage(x, c(0, 5), target_core_id = "T")
  msgs <- depth_coverage_messages(report, "T", c(0, 5))
  expect_length(msgs, 1)
  expect_match(msgs, "^Depth coverage for T's 2 sample depth\\(s\\) x 1 parameter\\(s\\): 1 measured, 1 interpolated, 0 out of retrieved range, 0 gaps too large\\.$")
})

test_that("depth_coverage_messages explains a depth mismatch", {
  x <- mk_pl("T", "TOC", c(300, 1000), c(1, 2))
  report <- report_depth_coverage(x, c(40, 200), target_core_id = "T")
  msgs <- depth_coverage_messages(report, "T", c(40, 200))
  expect_match(msgs[2], "your samples are at 40-200 cm, but retrieved data doesn't cover that range for: TOC \\(toc\\): retrieved 300-1000 cm")
})

test_that("depth_coverage_messages does not crash when match_basis is NA alongside a refused neighbour", {
  # The scenario that took the app down: neighbour A has an age model but its
  # data cannot resolve at the target's ages (so those cells have NA
  # match_basis), while neighbour B has no age model at all (refused).
  x <- dplyr::bind_rows(
    mk_pl("A", "TOC", c(0, 20), c(1, 2)), mk_pl("A", "age_model", c(0, 20), c(0, 20), "age"),
    mk_pl("B", "TOC", c(0, 20), c(1, 2))
  )
  report <- report_depth_coverage(x, c(50, 60), target_core_id = "T", sample_registry = sr(c(500, 600)))
  expect_true(anyNA(report$match_basis))

  msgs <- expect_no_error(depth_coverage_messages(report, "T", c(50, 60)))
  expect_true(any(grepl("^Skipped B: no age model available for cross-core matching\\.$", msgs)))
})

test_that("depth_coverage_messages points at XRF data when a refused neighbour has it", {
  x <- dplyr::bind_rows(
    mk_pl("B", "TOC", c(0, 20), c(1, 2)),
    mk_pl("B", "Ca_xrf", c(0, 20), c(5, 6), group = "xrf")
  )
  report <- report_depth_coverage(x, c(50, 60), target_core_id = "T", sample_registry = sr(c(5, 6)))
  msgs <- depth_coverage_messages(report, "T", c(50, 60))
  expect_true(any(grepl("Skipped B: .*XRF data is available for this core", msgs)))
})

test_that("depth_coverage_messages reports sample depths the target has no age for", {
  x <- dplyr::bind_rows(
    mk_pl("N", "TOC", c(0, 20), c(1, 2)), mk_pl("N", "age_model", c(0, 20), c(0, 20), "age")
  )
  report <- report_depth_coverage(x, c(50, 60), target_core_id = "T", sample_registry = sr(c(5, NA)))
  msgs <- depth_coverage_messages(report, "T", c(50, 60))
  expect_true(any(grepl("Skipped cross-core matching for 1 sample depth\\(s\\): no age_ka on file for T", msgs)))
})

# -- small conveniences -----------------------------------------------------------------------

test_that("doi_link builds an anchor, is NA-safe, and honours base_url", {
  expect_match(doi_link("10.1594/PANGAEA.1"), '^<a href="https://doi.pangaea.de/10.1594/PANGAEA.1" target="_blank"')
  expect_match(doi_link("10.5281/zenodo.1", "https://doi.org/"), 'href="https://doi.org/10.5281/zenodo.1"')
  expect_equal(doi_link(c("a", NA))[2], NA_character_)
})

test_that("iodp_expedition_from_core_id parses expedition-prefixed holes only, vectorised", {
  expect_equal(
    iodp_expedition_from_core_id(c("362-U1480A", "U1480A", "GeoB1706-2", "108-659A", "361-U1479")),
    c("362", NA, NA, NA, "361")
  )
})

test_that("export_filename sanitises the core id and falls back to 'export'", {
  d <- as.Date("2026-11-19")
  expect_equal(export_filename("GeoB 18549/1", "harmonized", d), "sedacore_GeoB_18549_1_2026-11-19.tsv")
  expect_equal(export_filename("T1", "depth_matched", d), "sedacore_depth_matched_T1_2026-11-19.tsv")
  expect_equal(export_filename(NULL, "harmonized", d), "sedacore_export_2026-11-19.tsv")
  expect_equal(export_filename("", "harmonized", d), "sedacore_export_2026-11-19.tsv")
  expect_equal(export_filename(NA_character_, "harmonized", d), "sedacore_export_2026-11-19.tsv")
})

test_that("load_core_registries reads the bundled example into both registries", {
  res <- load_core_registries(system.file("extdata", "example_core_metadata.tsv", package = "sedacore"))
  expect_true(inherits(res$sample_registry, "sample_registry"))
  expect_true(inherits(res$core_registry, "core_registry"))
  expect_gt(nrow(res$core_registry), 0)
})
