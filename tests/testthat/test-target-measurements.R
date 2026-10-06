make_row <- function(core_id, parameter_name, depth_cm, is_interpolated = FALSE) {
  new_param_long(tibble::tibble(
    core_id = core_id, pangaea_doi = "doi1", depth_cm = depth_cm,
    parameter_group = "toc", parameter_name = parameter_name,
    species = NA_character_, value = 1, unit = "%",
    raw_column_name = "raw", is_interpolated = is_interpolated,
    interp_gap_cm = NA_real_, retrieved_at = Sys.time()
  ))
}

test_that("load_target_measurements reads the bundled checklist", {
  targets <- load_target_measurements()
  expect_true(all(c("parameter_name", "notes") %in% names(targets)))
  expect_true("TOC" %in% targets$parameter_name)
})

test_that("check_target_measurements marks present parameters as found with correct stats", {
  x <- dplyr::bind_rows(
    make_row("core1", "TOC", 10),
    make_row("core1", "TOC", 20)
  )
  targets <- tibble::tibble(parameter_name = c("TOC", "TN"), notes = c("n1", "n2"))

  out <- check_target_measurements(x, "core1", targets)

  expect_equal(out$status[out$parameter_name == "TOC"], "found")
  expect_equal(out$n_measured[out$parameter_name == "TOC"], 2)
  expect_equal(out$depth_min_cm[out$parameter_name == "TOC"], 10)
  expect_equal(out$depth_max_cm[out$parameter_name == "TOC"], 20)
  expect_equal(out$status[out$parameter_name == "TN"], "missing")
  expect_equal(out$n_measured[out$parameter_name == "TN"], 0)
  expect_true(is.na(out$depth_min_cm[out$parameter_name == "TN"]))
})

test_that("check_target_measurements ignores interpolated rows and other cores' data", {
  x <- dplyr::bind_rows(
    make_row("core1", "TOC", 10, is_interpolated = TRUE),
    make_row("core2", "TOC", 10)
  )
  targets <- tibble::tibble(parameter_name = "TOC", notes = "n1")

  out <- check_target_measurements(x, "core1", targets)
  expect_equal(out$status, "missing")
})

test_that("check_target_measurements returns all-missing for a core with no data at all", {
  x <- make_row("core2", "TOC", 10)
  targets <- tibble::tibble(parameter_name = c("TOC", "TN"), notes = c("n1", "n2"))

  out <- check_target_measurements(x, "core1", targets)
  expect_true(all(out$status == "missing"))
  expect_equal(out$n_measured, c(0, 0))
})

test_that("check_target_measurements falls back to neighbour data when the core has none", {
  x <- dplyr::bind_rows(
    make_row("core1", "TOC", 10),
    make_row("neighbourA", "TN", 15)
  )
  targets <- tibble::tibble(parameter_name = c("TOC", "TN", "CaCO3"), notes = c("n1", "n2", "n3"))

  out <- check_target_measurements(x, "core1", targets, neighbour_core_ids = "neighbourA")

  expect_equal(out$status[out$parameter_name == "TOC"], "found")
  expect_equal(out$source_core_id[out$parameter_name == "TOC"], "core1")
  expect_equal(out$status[out$parameter_name == "TN"], "found in neighbour")
  expect_equal(out$source_core_id[out$parameter_name == "TN"], "neighbourA")
  expect_equal(out$n_measured[out$parameter_name == "TN"], 1)
  expect_equal(out$depth_min_cm[out$parameter_name == "TN"], 15)
  expect_equal(out$status[out$parameter_name == "CaCO3"], "missing")
  expect_true(is.na(out$source_core_id[out$parameter_name == "CaCO3"]))
})

test_that("check_target_measurements prefers the core's own data over a neighbour's for the same parameter", {
  x <- dplyr::bind_rows(
    make_row("core1", "TOC", 10),
    make_row("neighbourA", "TOC", 999)
  )
  targets <- tibble::tibble(parameter_name = "TOC", notes = "n1")

  out <- check_target_measurements(x, "core1", targets, neighbour_core_ids = "neighbourA")
  expect_equal(out$status, "found")
  expect_equal(out$source_core_id, "core1")
  expect_equal(out$depth_min_cm, 10)
})

test_that("check_target_measurements ignores an interpolated row in a neighbour core", {
  x <- make_row("neighbourA", "TOC", 10, is_interpolated = TRUE)
  targets <- tibble::tibble(parameter_name = "TOC", notes = "n1")

  out <- check_target_measurements(x, "core1", targets, neighbour_core_ids = "neighbourA")
  expect_equal(out$status, "missing")
})

test_that("the bundled checklist includes sea-ice proxies, each resolvable by a real synonym", {
  targets <- load_target_measurements()
  expect_true(all(c("IRD_count", "IP25_sed", "PBIP25", "sea_ice_conc", "F_curta_gr") %in% targets$parameter_name))

  # Every name the checklist can match must actually be producible by
  # harmonize_columns(), or the entry would sit at "missing" forever.
  produced <- unique(load_synonym_map()$parameter_name)
  expect_true(all(target_measurement_names(targets) %in% produced))
})

test_that("an also_matches alternative satisfies its checklist entry and reports found_as", {
  x <- make_row("core1", "IRD_per_g", 10)
  targets <- tibble::tibble(
    parameter_name = c("IRD_count", "TOC"), notes = c("n1", "n2"),
    also_matches = c("IRD_per_g|IRD_percent", NA)
  )

  out <- check_target_measurements(x, "core1", targets)

  expect_equal(out$status[out$parameter_name == "IRD_count"], "found")
  expect_equal(out$found_as[out$parameter_name == "IRD_count"], "IRD_per_g")
  expect_equal(out$n_measured[out$parameter_name == "IRD_count"], 1)
  expect_equal(out$status[out$parameter_name == "TOC"], "missing")
  expect_true(is.na(out$found_as[out$parameter_name == "TOC"]))
})

test_that("several variants of one checklist entry are pooled and all listed in found_as", {
  x <- dplyr::bind_rows(
    make_row("core1", "IRD_per_g", 10),
    make_row("core1", "IRD_percent", 30)
  )
  targets <- tibble::tibble(
    parameter_name = "IRD_count", notes = "n1", also_matches = "IRD_per_g|IRD_percent"
  )

  out <- check_target_measurements(x, "core1", targets)

  expect_equal(out$n_measured, 2)
  expect_equal(out$depth_min_cm, 10)
  expect_equal(out$depth_max_cm, 30)
  # Order is locale-dependent (sort()), so compare as a set.
  expect_setequal(strsplit(out$found_as, ", ")[[1]], c("IRD_per_g", "IRD_percent"))
})

test_that("also_matches also works for the neighbour fallback", {
  x <- make_row("neighbourA", "IP25_toc", 15)
  targets <- tibble::tibble(parameter_name = "IP25_sed", notes = "n1", also_matches = "IP25_toc")

  out <- check_target_measurements(x, "core1", targets, neighbour_core_ids = "neighbourA")

  expect_equal(out$status, "found in neighbour")
  expect_equal(out$found_as, "IP25_toc")
})

test_that("target_measurement_names lists headline names plus alternatives, without duplicates", {
  targets <- tibble::tibble(
    parameter_name = c("IRD_count", "TOC"), notes = c("a", "b"),
    also_matches = c("IRD_per_g|IRD_percent", NA)
  )
  expect_setequal(target_measurement_names(targets), c("IRD_count", "IRD_per_g", "IRD_percent", "TOC"))
  # A checklist with no also_matches column at all still works.
  expect_setequal(target_measurement_names(targets[, c("parameter_name", "notes")]), c("IRD_count", "TOC"))
})
