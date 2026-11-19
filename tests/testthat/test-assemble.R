palmod_fixtures <- function(files = c("GeoB3938_1_age.RDS", "MD95_2024_age.RDS")) {
  build_palmod(test_path("..", "fixtures", "palmod"), files)
}

# A pretend second source: the same GeoB3938 core spelled the way your own core
# list spells it (with a hyphen), plus one site only it has. Series ids are
# prefixed, as a real adapter must do, so they stay unique across sources.
fake_second_source <- function(p) {
  s <- p$tables
  keep <- s$sites$site_id == "GeoB3938_1"
  sites <- s$sites[keep, ]
  sites$site_id <- sites$label <- "GeoB3938-1"
  sites$latitude <- sites$latitude + 0.001
  sites$sources <- "wpa_test"
  sites <- dplyr::bind_rows(sites, tibble::tibble(
    site_id = "WPAONLY-1", label = "WPAONLY-1", latitude = -10, longitude = 20, water_depth_m = 3000,
    basin = "South Atlantic Ocean", program = "other", sources = "wpa_test"
  ))
  recode <- function(x) {
    x <- x[x$site_id == "GeoB3938_1", ]
    x$site_id <- "GeoB3938-1"
    if ("series_id" %in% names(x)) x$series_id <- paste0("wpa:", x$series_id)
    if ("source" %in% names(x)) x$source <- "wpa_test"
    x
  }
  list(
    tables = list(
      sites = sites, series = recode(s$series), series_catalog = recode(s$series_catalog),
      age_control = recode(s$age_control)[0, ], site_links = recode(s$site_links)
    ),
    source = list(id = "wpa_test", name = "Fake atlas", version = "0", doi = "10.0/fake", license = "CC0", citation = "c")
  )
}

test_that("assemble_bundle folds two sources onto canonical sites and re-keys every table", {
  p <- palmod_fixtures()
  a <- assemble_bundle(list(p, fake_second_source(p)), priority = c("palmod_v2", "wpa_test"))
  t <- a$tables

  expect_setequal(t$sites$site_id, c("GeoB3938-1", "MD95-2024", "WPAONLY-1"))
  expect_false("n_aliases" %in% names(t$sites))
  expect_true(all(c("GeoB3938-1", "MD95-2024") %in% unique(t$series$site_id)))
  expect_false(any(grepl("_", c(t$series$site_id, t$series_catalog$site_id, t$age_control$site_id, t$site_links$site_id))))
  expect_equal(t$site_aliases$site_id[t$site_aliases$alias == "GeoB3938_1"], "GeoB3938-1")
  expect_equal(t$site_aliases$site_id[t$site_aliases$alias == "GeoB3938-1"], "GeoB3938-1")
  for (nm in names(t)) expect_no_error(check_bundle_table(t[[nm]], nm))
})

test_that("a measurement held by both sources is counted once: PALMOD's copy is preferred", {
  p <- palmod_fixtures()
  a <- assemble_bundle(list(p, fake_second_source(p)), priority = c("palmod_v2", "wpa_test"))
  cat <- a$tables$series_catalog
  g <- cat[cat$site_id == "GeoB3938-1" & cat$parameter_name == "d18O_benthic", ]
  expect_equal(nrow(g), 2)
  expect_equal(g$source[g$preferred], "palmod_v2")
  expect_equal(sum(g$preferred), 1)

  # the flag is carried onto the observations
  s <- a$tables$series
  expect_equal(sort(unique(s$source[s$site_id == "GeoB3938-1" & s$preferred])), "palmod_v2")
  expect_true(any(!s$preferred))
  # exactly one preferred series per (site, parameter, species) among central series
  central <- cat[cat$role == "central", ]
  key <- paste(central$site_id, central$parameter_name, normalize_species(central$species))
  expect_equal(as.integer(tapply(central$preferred, key, sum)), rep(1L, length(unique(key))))

  # swap the priority and the WPA copy wins instead
  b <- assemble_bundle(list(p, fake_second_source(p)), priority = c("wpa_test", "palmod_v2"))
  g2 <- b$tables$series_catalog[b$tables$series_catalog$site_id == "GeoB3938-1" &
    b$tables$series_catalog$parameter_name == "d18O_benthic", ]
  expect_equal(g2$source[g2$preferred], "wpa_test")
})

test_that("duplicate links that appear only after canonical re-keying are dropped", {
  p <- palmod_fixtures()
  a <- assemble_bundle(list(p, fake_second_source(p)), priority = c("palmod_v2", "wpa_test"))
  l <- a$tables$site_links
  expect_equal(anyDuplicated(l[, c("site_id", "link_type", "url")]), 0)
})

test_that("assemble_bundle refuses a source missing from `priority`, and duplicate source ids", {
  p <- palmod_fixtures()
  expect_error(assemble_bundle(list(p, fake_second_source(p)), priority = "palmod_v2"), "missing from `priority`")
  expect_error(assemble_bundle(list(p, p), priority = "palmod_v2"), "must be unique")
})

test_that("series ids that collide across sources are refused", {
  p <- palmod_fixtures()
  w <- fake_second_source(p)
  w$tables$series_catalog$series_id <- sub("^wpa:", "", w$tables$series_catalog$series_id)
  w$tables$series$series_id <- sub("^wpa:", "", w$tables$series$series_id)
  expect_error(assemble_bundle(list(p, w), priority = c("palmod_v2", "wpa_test")), "must be unique")
})

test_that("an assembled bundle round-trips and serves only preferred series by default", {
  p <- palmod_fixtures()
  a <- assemble_bundle(list(p, fake_second_source(p)), priority = c("palmod_v2", "wpa_test"))
  dir <- withr::local_tempdir()
  write_bundle(a$tables, dir, sources = a$sources, crosswalked = TRUE)
  b <- sedacore_open_bundle(dir)
  withr::defer(close_bundle(b))

  expect_true(b$manifest$crosswalked)
  expect_equal(nrow(bundle_table(b, "site_aliases")), nrow(a$tables$site_aliases))
  expect_equal(vapply(b$manifest$sources, `[[`, "", "id"), c("palmod_v2", "wpa_test"))

  s_pref <- bundle_series(b, site_id = "GeoB3938-1", parameter_name = "d18O_benthic")
  s_all <- bundle_series(b, site_id = "GeoB3938-1", parameter_name = "d18O_benthic", preferred_only = FALSE)
  expect_equal(unique(s_pref$source), "palmod_v2")
  expect_equal(nrow(s_all), 2 * nrow(s_pref))
})

test_that("preferred_only is refused on a bundle that was not built with a crosswalk", {
  dir <- withr::local_tempdir()
  write_bundle(list(sites = tibble::tibble(
    site_id = "A", label = "A", latitude = 1, longitude = 1, water_depth_m = 1, basin = "x", program = "y", sources = "z"
  )), dir)
  b <- sedacore_open_bundle(dir)
  withr::defer(close_bundle(b))
  expect_false(isTRUE(b$manifest$crosswalked))
  expect_error(bundle_series(b, preferred_only = TRUE), "no `preferred` flag")
})

test_that("on PALMOD alone, holes of one drilling site fold together (real fixtures have none; synthetic does)", {
  p <- palmod_fixtures("GeoB3938_1_age.RDS")
  hole_b <- p
  hole_b$tables$sites$site_id <- hole_b$tables$sites$label <- "108_658C"
  hole_b$tables$sites$latitude <- 20.7493; hole_b$tables$sites$longitude <- -18.5809
  hole_a <- p
  hole_a$tables$sites$site_id <- hole_a$tables$sites$label <- "108_658"
  hole_a$tables$sites$latitude <- 20.7492; hole_a$tables$sites$longitude <- -18.5808
  for (nm in c("series", "series_catalog", "age_control", "site_links")) {
    hole_a$tables[[nm]]$site_id <- "108_658"
    hole_b$tables[[nm]]$site_id <- "108_658C"
    if (nm %in% c("series", "series_catalog")) hole_b$tables[[nm]]$series_id <- paste0("c:", hole_b$tables[[nm]]$series_id)
  }
  hole_a$source$id <- "palmod_a"; hole_b$source$id <- "palmod_b"
  a <- assemble_bundle(list(hole_a, hole_b), c("palmod_a", "palmod_b"))
  expect_equal(a$tables$sites$site_id, "ODP-108-658")
  expect_equal(nrow(a$tables$site_aliases), 2)
})
