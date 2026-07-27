test_that("bbox_from_radius returns a plausible box around the centre", {
  bb <- bbox_from_radius(latitude = -19.561667, longitude = 11.175, radius_km = 100)

  expect_length(bb, 4)
  expect_true(bb[1] < 11.175 && bb[3] > 11.175)
  expect_true(bb[2] < -19.561667 && bb[4] > -19.561667)
})

test_that("circle_filter keeps only points within radius and drops missing coords", {
  x <- tibble::tibble(
    id = c("near", "far", "missing"),
    latitude = c(-19.6, 10, NA),
    longitude = c(11.2, 10, NA)
  )

  out <- circle_filter(x, center_lat = -19.561667, center_lon = 11.175, radius_km = 100)

  expect_equal(out$id, "near")
  expect_true(out$distance_km < 100)
})

test_that("find_neighbours finds real nearby ODP sites around GeoB1706-2, excluding self", {
  cr <- new_core_registry(tibble::tibble(
    core_id = "GeoB1706-2", core_label = "GeoB1706-2; 2003",
    latitude = -19.561667, longitude = 11.175, water_depth_m = 980,
    source_program = "GeoB", n_samples = 18L,
    depth_min_cm = 5, depth_max_cm = 160, pangaea_doi = NA_character_
  ))

  vcr::use_cassette("find_neighbours_geob1706", {
    nb <- find_neighbours(cr, "GeoB1706-2", radius_km = 100, count = 30)
  })

  expect_true(is_neighbour_registry(nb))
  expect_true(nrow(nb) > 0)
  expect_true(all(nb$distance_km <= 100))
  expect_false(any(nb$neighbour_core_id == "GeoB1706-2"))
  expect_true(any(grepl("^175-1081|^75-532", nb$neighbour_core_id)))
  expect_true(all(nb$source_program %in% c("ODP", "other", "GeoB", "IODP", "DSDP")))
  expect_true(all(c("latitude", "longitude") %in% names(nb)))
  expect_true(all(!is.na(nb$latitude) & !is.na(nb$longitude)))
  expect_true("citation" %in% names(nb))
  expect_true(all(!is.na(nb$citation) & nzchar(nb$citation)))
})

test_that("normalize_core_label matches a plain core_id against a PANGAEA label with a station suffix", {
  # Real bug: PANGAEA's event_label for a core already in the user's file
  # can carry a cruise station suffix (e.g. "GeoB18549-2 (528-2)") that an
  # exact-string comparison against the plain core_id ("GeoB18549-2")
  # fails to recognize as the same core, so it wasn't actually excluded.
  expect_equal(normalize_core_label("GeoB18549-2"), normalize_core_label("GeoB18549-2 (528-2)"))
  expect_equal(normalize_core_label("175-1081A"), "175-1081a")
  expect_true(is.na(normalize_core_label(NA_character_)))
})

test_that("find_neighbours excludes cores already present anywhere in x, not just core_id itself", {
  # Same core registry as above, but with a second row for a core we know
  # (from the cassette above) shows up as a real candidate near
  # GeoB1706-2 -- confirms it's excluded because it's already in the
  # user's own dataset, mirroring the GeoB18549-1/-2 case.
  cr <- new_core_registry(tibble::tibble(
    core_id = c("GeoB1706-2", "175-1081A"),
    core_label = c("GeoB1706-2; 2003", "175-1081A"),
    latitude = c(-19.561667, -19.4), longitude = c(11.175, 11.2),
    water_depth_m = c(980, 4000), source_program = c("GeoB", "IODP"),
    n_samples = c(18L, 1L), depth_min_cm = c(5, 0), depth_max_cm = c(160, 0),
    pangaea_doi = NA_character_
  ))

  vcr::use_cassette("find_neighbours_geob1706", {
    nb <- find_neighbours(cr, "GeoB1706-2", radius_km = 100, count = 30)
  })

  expect_false(any(nb$neighbour_core_id == "175-1081A"))
})

test_that("find_neighbours calls the progress callback once per candidate", {
  cr <- new_core_registry(tibble::tibble(
    core_id = "GeoB1706-2", core_label = "GeoB1706-2; 2003",
    latitude = -19.561667, longitude = 11.175, water_depth_m = 980,
    source_program = "GeoB", n_samples = 18L,
    depth_min_cm = 5, depth_max_cm = 160, pangaea_doi = NA_character_
  ))
  calls <- list()

  vcr::use_cassette("find_neighbours_geob1706", {
    nb <- find_neighbours(
      cr, "GeoB1706-2", radius_km = 100, count = 30,
      progress = function(i, n, doi) calls[[length(calls) + 1]] <<- list(i = i, n = n, doi = doi)
    )
  })

  expect_equal(length(calls), 30)
  expect_equal(calls[[1]]$i, 1)
  expect_equal(calls[[30]]$i, 30)
  expect_true(all(vapply(calls, function(x) x$n == 30, logical(1))))
})

test_that("find_neighbours filters results to the requested source_programs", {
  cr <- new_core_registry(tibble::tibble(
    core_id = "GeoB1706-2", core_label = "GeoB1706-2; 2003",
    latitude = -19.561667, longitude = 11.175, water_depth_m = 980,
    source_program = "GeoB", n_samples = 18L,
    depth_min_cm = 5, depth_max_cm = 160, pangaea_doi = NA_character_
  ))

  vcr::use_cassette("find_neighbours_geob1706", {
    odp_only <- find_neighbours(cr, "GeoB1706-2", radius_km = 100, count = 30, source_programs = "ODP")
  })
  expect_true(is_neighbour_registry(odp_only))
  expect_true(nrow(odp_only) > 0)
  expect_true(all(odp_only$source_program == "ODP"))

  vcr::use_cassette("find_neighbours_geob1706", {
    none <- find_neighbours(cr, "GeoB1706-2", radius_km = 100, count = 30, source_programs = "DSDP")
  })
  expect_true(is_neighbour_registry(none))
  expect_equal(nrow(none), 0)
})

test_that("order_neighbours ranks GeoB first, then IODP, then others, by recency within each", {
  x <- tibble::tibble(
    neighbour_core_id = c("old-odp", "new-other", "old-geob", "new-geob", "iodp-only"),
    source_program = c("ODP", "other", "GeoB", "GeoB", "IODP"),
    date_time = as.POSIXct(c("2001-01-01", "2020-01-01", "1995-01-01", "2015-01-01", NA), tz = "UTC")
  )

  out <- order_neighbours(x)

  expect_equal(out$neighbour_core_id, c("new-geob", "old-geob", "iodp-only", "old-odp", "new-other"))
})

test_that("find_neighbours errors on an unknown core_id", {
  cr <- new_core_registry(tibble::tibble(
    core_id = "core1", core_label = "core1", latitude = 1, longitude = 1,
    water_depth_m = 1, source_program = "other", n_samples = 1L,
    depth_min_cm = 1, depth_max_cm = 1, pangaea_doi = NA_character_
  ))
  expect_error(find_neighbours(cr, "not-a-core", radius_km = 10))
})
