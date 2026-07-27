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
