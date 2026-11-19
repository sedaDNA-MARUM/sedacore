mk_sites <- function() {
  tibble::tibble(
    site_id = c("A", "B"), label = c("A", "B"), latitude = c(10, -20), longitude = c(-30, 5),
    water_depth_m = c(1000, NA), basin = c("North Atlantic Ocean", "South Atlantic Ocean"),
    program = c("GeoB", "other"), sources = "test"
  )
}

mk_series <- function() {
  tibble::tibble(
    site_id = c("A", "A", "A", "B"), series_id = c("s1", "s1", "s2", "s3"), source = "test",
    parameter_name = c("TOC", "TOC", "d18O_benthic", "TOC"),
    proxy_class = c("organic carbon", "organic carbon", "d18O benthic", "organic carbon"),
    species = NA_character_, material = NA_character_,
    role = c("central", "central", "central", "lower"), is_recalc = FALSE,
    depth_cm = c(10, 20, 10, 5), age_ka = c(1, 2, 1, 3), age_lo_ka = NA_real_, age_hi_ka = NA_real_,
    value = c(1.5, 2.5, 3.5, 4.5), value_lo = NA_real_, value_hi = NA_real_, unit = "percent"
  )
}

new_test_bundle <- function(dir = withr::local_tempdir(.local_envir = parent.frame()), ...) {
  write_bundle(
    list(sites = mk_sites(), series = mk_series()), dir,
    sources = list(list(id = "test", name = "Test source", version = "1", doi = "10.0/x", license = "CC0", citation = "c")),
    ...
  )
  dir
}

test_that("a bundle round-trips: tables, row counts and manifest survive", {
  dir <- new_test_bundle()
  b <- sedacore_open_bundle(dir)
  withr::defer(close_bundle(b))

  expect_s3_class(b, "sedacore_bundle")
  expect_equal(b$manifest$bundle_format, 1)
  expect_equal(b$manifest$tables$series$rows, 4)
  expect_equal(b$manifest$sources[[1]]$license, "CC0")
  expect_equal(b$manifest$package, "sedacore")

  sites <- bundle_table(b, "sites")
  expect_equal(sites$site_id, c("A", "B"))
  expect_true(is.na(sites$water_depth_m[2]))
  expect_equal(nrow(bundle_table(b, "series")), 4)
})

test_that("bundle_series filters by site, parameter, age range and role, in age order", {
  dir <- new_test_bundle()
  b <- sedacore_open_bundle(dir)
  withr::defer(close_bundle(b))

  expect_equal(nrow(bundle_series(b)), 3) # default role = central: B's 'lower' row excluded
  expect_equal(nrow(bundle_series(b, role = NULL)), 4)
  expect_equal(unique(bundle_series(b, site_id = "A")$site_id), "A")
  expect_equal(bundle_series(b, parameter_name = "d18O_benthic")$value, 3.5)
  expect_equal(bundle_series(b, age_range = c(1.5, 5))$value, 2.5)
  expect_equal(bundle_series(b, site_id = "A", parameter_name = "TOC")$age_ka, c(1, 2))
  expect_equal(nrow(bundle_series(b, site_id = c("A", "B"), role = c("central", "lower"))), 4)
})

test_that("bundle queries are parameterised, not string-pasted (no injection via a site id)", {
  dir <- new_test_bundle()
  b <- sedacore_open_bundle(dir)
  withr::defer(close_bundle(b))

  nasty <- "A'; DROP VIEW series; --"
  expect_equal(nrow(bundle_series(b, site_id = nasty)), 0)
  expect_equal(nrow(bundle_series(b)), 3) # the view is still there
})

test_that("bundle_table refuses a table that is not in the bundle", {
  b <- sedacore_open_bundle(new_test_bundle())
  withr::defer(close_bundle(b))
  expect_error(bundle_table(b, "nope"), "not in this bundle")
})

test_that("write_bundle validates against the schema", {
  d <- withr::local_tempdir()
  expect_error(write_bundle(list(sites = mk_sites()[, -1]), d), "missing required column\\(s\\): site_id")
  bad <- mk_sites()
  bad$latitude <- as.character(bad$latitude)
  expect_error(write_bundle(list(sites = bad), d), "wrong type: latitude \\(expected numeric\\)")
  expect_error(write_bundle(list(sitez = mk_sites()), d), "Unknown bundle table")
  expect_error(write_bundle(list(mk_sites()), d), "named list")
})

test_that("extra columns are allowed (how a source adds detail without a schema change)", {
  s <- mk_sites()
  s$palmod_file <- c("a.RDS", "b.RDS")
  d <- withr::local_tempdir()
  write_bundle(list(sites = s), d)
  b <- sedacore_open_bundle(d)
  withr::defer(close_bundle(b))
  expect_true("palmod_file" %in% names(bundle_table(b, "sites")))
})

test_that("write_bundle will not silently overwrite an existing bundle", {
  dir <- new_test_bundle()
  expect_error(write_bundle(list(sites = mk_sites()), dir), "already exists")
  expect_no_error(write_bundle(list(sites = mk_sites()), dir, overwrite = TRUE))
})

test_that("opening verifies checksums, so an edited or truncated file is caught", {
  dir <- new_test_bundle()
  f <- file.path(dir, "series.parquet")
  bytes <- readBin(f, "raw", file.size(f))
  writeBin(c(bytes, as.raw(0)), f)
  expect_error(sedacore_open_bundle(dir), "Checksum mismatch for series.parquet")
  # Even without the checksum, DuckDB refuses a damaged parquet file.
  expect_error(sedacore_open_bundle(dir, verify = FALSE), "magic bytes|parquet")
})

test_that("opening reports a missing table file and an unsupported format clearly", {
  dir <- new_test_bundle()
  file.remove(file.path(dir, "sites.parquet"))
  expect_error(sedacore_open_bundle(dir), "table file is missing")

  dir2 <- new_test_bundle()
  m <- jsonlite::read_json(file.path(dir2, "manifest.json"))
  m$bundle_format <- 99L
  jsonlite::write_json(m, file.path(dir2, "manifest.json"), auto_unbox = TRUE)
  expect_error(sedacore_open_bundle(dir2), "Unsupported bundle format 99")
})

test_that("sedacore_bundle_path resolves explicit path, option, env var, and fails helpfully", {
  dir <- new_test_bundle()
  expect_equal(sedacore_bundle_path(dir), normalizePath(dir))

  withr::with_options(list(sedacore.bundle = dir), expect_equal(sedacore_bundle_path(), normalizePath(dir)))
  withr::with_envvar(c(SEDACORE_BUNDLE = dir), expect_equal(sedacore_bundle_path(), normalizePath(dir)))

  empty <- withr::local_tempdir()
  withr::with_options(list(sedacore.bundle = empty), withr::with_envvar(c(SEDACORE_BUNDLE = ""), {
    # (falls through to the cache dir, which has no bundle in the test environment)
    expect_error(sedacore_bundle_path(empty), "No sedacore data bundle found")
  }))
})

test_that("the cache dir is derived from the package name (so a rename needs no edit)", {
  expect_match(basename(sedacore_cache_dir()), "^sedacore$")
})

test_that("printing a bundle summarises it", {
  b <- sedacore_open_bundle(new_test_bundle())
  withr::defer(close_bundle(b))
  out <- paste(capture.output(print(b)), collapse = "\n")
  expect_match(out, "<sedacore_bundle>")
  expect_match(out, "series +4 rows")
  expect_match(out, "sources: test 1")
})

test_that("closing a bundle twice is a quiet no-op (it is also closed by the finalizer)", {
  b <- sedacore_open_bundle(new_test_bundle())
  expect_no_warning(close_bundle(b))
  expect_no_warning(close_bundle(b))
  expect_false(DBI::dbIsValid(b$con))
})
