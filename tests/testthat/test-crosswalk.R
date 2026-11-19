mk_site <- function(id, lat, lon, source, depth = NA_real_, basin = NA_character_) {
  tibble::tibble(
    site_id = id, label = id, latitude = lat, longitude = lon, water_depth_m = depth,
    basin = basin, program = NA_character_, sources = source, source = source
  )
}

test_that("haversine_km gives known distances", {
  expect_equal(haversine_km(0, 0, 0, 1), 111.195, tolerance = 1e-3)
  expect_equal(haversine_km(0, 0, 0, 0), 0)
  expect_equal(haversine_km(50.205667, -45.6858, 50.21, -45.69), 0.57, tolerance = 0.03) # GeoB18549 vs MD95-2024
  expect_equal(haversine_km(c(0, 0), c(0, 0), c(0, 10), c(1, 0))[2], 1111.95, tolerance = 1e-3) # vectorised
})

test_that("site_name_key ignores case and separators, and keeps NA", {
  expect_equal(site_name_key(c("GeoB3938_1", "GeoB3938-1", "geob3938 1", "MD95_2024", NA)),
    c("geob39381", "geob39381", "geob39381", "md952024", NA))
})

test_that("iodp_site_key recognises ocean-drilling sites in every spelling, and nothing else", {
  expect_equal(
    iodp_site_key(c("108_658C", "108-658", "ODP 658B", "ODP Site 658", "DSDP Site 502", "94_502", "658")),
    c("dsdp-odp-658", "dsdp-odp-658", "dsdp-odp-658", "dsdp-odp-658", "dsdp-odp-502", "dsdp-odp-502", NA)
  )
  expect_equal(
    iodp_site_key(c("303_U1308", "U1308A", "IODP Site U1308", "361-U1479", "u1308")),
    rep(c("iodp-U1308", "iodp-U1479")[c(1, 1, 1, 2, 1)])
  )
  # not drilling sites: ordinary core names, and "leg" numbers beyond ODP
  expect_true(all(is.na(iodp_site_key(c("GeoB3938_1", "MD95_2024", "GIK15666_6", "KNR166_2_2", "400_1234", NA)))))
})

test_that("canonical_site_id builds readable ids", {
  expect_equal(canonical_site_id(c("108_658C", "ODP 658B")), "ODP-108-658")
  expect_equal(canonical_site_id(c("ODP 658B", "108_658C")), "ODP-108-658") # the leg is used whenever any member has it
  expect_equal(canonical_site_id("ODP 658B"), "ODP-658")
  expect_equal(canonical_site_id("94_502"), "DSDP-94-502") # legs up to 96 are DSDP
  expect_equal(canonical_site_id(c("303_U1308", "U1308A")), "IODP-U1308")
  expect_equal(canonical_site_id(c("GeoB3938_1", "GeoB3938-1")), "GeoB3938-1")
  expect_equal(canonical_site_id("KNR166_2_2"), "KNR166-2-2")
})

test_that("the same core under different spellings and sources becomes one site", {
  sites <- dplyr::bind_rows(
    mk_site("GeoB3938_1", 12.25833, -58.33, "palmod", 1972, "North Atlantic Ocean"),
    mk_site("MD95_2024", 50.21, -45.69, "palmod", 3448, "North Atlantic Ocean"),
    mk_site("GeoB3938-1", 12.2584, -58.3302, "wpa", NA, NA),
    mk_site("MD95-2024", 50.2104, -45.6905, "wpa", 3450, "North Atlantic Ocean")
  )
  cw <- build_crosswalk(sites, priority = c("palmod", "wpa"))
  expect_equal(sort(cw$sites$site_id), c("GeoB3938-1", "MD95-2024"))
  expect_equal(cw$sites$n_aliases, c(2, 2))
  g <- cw$sites[cw$sites$site_id == "GeoB3938-1", ]
  expect_equal(g$sources, "palmod,wpa")
  expect_equal(g$latitude, 12.25833) # position from the higher-priority source
  expect_equal(g$water_depth_m, 1972)
  expect_equal(nrow(cw$review), 0)

  a <- cw$site_aliases
  expect_equal(a$site_id[a$alias == "GeoB3938_1" & a$source == "palmod"], "GeoB3938-1")
  expect_equal(a$site_id[a$alias == "GeoB3938-1" & a$source == "wpa"], "GeoB3938-1")
  expect_equal(a$match_method[a$alias == "GeoB3938_1"], "representative")
  expect_equal(a$match_method[a$alias == "GeoB3938-1" & a$source == "wpa"], "name_key")
})

test_that("a missing field in the top-priority record is filled from the next one", {
  sites <- dplyr::bind_rows(
    mk_site("X_1", 1, 1, "palmod", NA, NA),
    mk_site("X-1", 1.001, 1, "wpa", 2500, "South Atlantic Ocean")
  )
  s <- build_crosswalk(sites, c("palmod", "wpa"))$sites
  expect_equal(s$water_depth_m, 2500)
  expect_equal(s$basin, "South Atlantic Ocean")
})

test_that("IODP holes and spellings of one drilling site merge, across sources", {
  sites <- dplyr::bind_rows(
    mk_site("108_658", 20.7492, -18.5808, "palmod"),
    mk_site("108_658C", 20.7493, -18.5809, "palmod"),
    mk_site("ODP 658B", 20.7490, -18.5810, "wpa"),
    mk_site("303_U1308", 49.88, -24.24, "palmod"),
    mk_site("U1308A", 49.8801, -24.2401, "wpa")
  )
  cw <- build_crosswalk(sites, c("palmod", "wpa"))
  expect_setequal(cw$sites$site_id, c("ODP-108-658", "IODP-U1308"))
  expect_equal(cw$sites$n_aliases[cw$sites$site_id == "ODP-108-658"], 3)
  expect_equal(unique(cw$site_aliases$site_id[cw$site_aliases$alias %in% c("108_658", "108_658C", "ODP 658B")]), "ODP-108-658")
  expect_equal(cw$sites$program[cw$sites$site_id == "ODP-108-658"], "ODP")
})

test_that("two different cores a few hundred metres apart are NOT merged -- they are neighbours, and are listed for review", {
  sites <- dplyr::bind_rows(
    mk_site("MD95_2024", 50.21, -45.69, "palmod"),
    mk_site("GeoB18549-1", 50.205667, -45.6858, "mine")
  )
  cw <- build_crosswalk(sites, c("palmod", "mine"))
  expect_equal(nrow(cw$sites), 2)
  expect_equal(cw$review$type, "near_different_name")
  expect_match(cw$review$detail, "0\\.5\\d km apart")
})

test_that("the same name far apart stays two sites and is flagged (separator-insensitive names can collide)", {
  sites <- dplyr::bind_rows(
    mk_site("GeoB1-23", 10, -20, "palmod"),
    mk_site("GeoB12-3", -40, 10, "wpa") # reduces to the same key, 5000+ km away
  )
  cw <- build_crosswalk(sites, c("palmod", "wpa"))
  expect_equal(nrow(cw$sites), 2)
  expect_equal(anyDuplicated(cw$sites$site_id), 0)
  expect_true("same_name_far_apart" %in% cw$review$type)
})

test_that("canonical site ids are unique and every input record has exactly one alias", {
  sites <- dplyr::bind_rows(
    mk_site("A_1", 1, 1, "palmod"), mk_site("A-1", 1.01, 1, "wpa"), mk_site("B_2", 5, 5, "palmod"),
    mk_site("108_658", 20, -18, "palmod"), mk_site("ODP 658A", 20, -18, "wpa")
  )
  cw <- build_crosswalk(sites, c("palmod", "wpa"))
  expect_equal(anyDuplicated(cw$sites$site_id), 0)
  expect_equal(nrow(cw$site_aliases), nrow(sites))
  expect_equal(anyDuplicated(cw$site_aliases[, c("source", "alias")]), 0)
  expect_true(all(cw$site_aliases$site_id %in% cw$sites$site_id))
})

test_that("source priority decides the position and the order labels are offered", {
  sites <- dplyr::bind_rows(mk_site("X-1", 1, 1, "wpa"), mk_site("X_1", 1.01, 1, "palmod"))
  s1 <- build_crosswalk(sites, c("palmod", "wpa"))$sites
  s2 <- build_crosswalk(sites, c("wpa", "palmod"))$sites
  expect_equal(s1$latitude, 1.01)
  expect_equal(s2$latitude, 1)
})

test_that("apply_crosswalk re-keys a source's tables, and refuses an id it cannot map", {
  sites <- dplyr::bind_rows(mk_site("GeoB3938_1", 12, -58, "palmod"), mk_site("GeoB3938-1", 12, -58, "wpa"))
  cw <- build_crosswalk(sites, c("palmod", "wpa"))
  tables <- list(
    sites = sites[1, ],
    series = tibble::tibble(site_id = c("GeoB3938_1", "GeoB3938_1"), series_id = c("a", "b")),
    site_links = tibble::tibble(site_id = "GeoB3938_1", url = "u")
  )
  out <- apply_crosswalk(tables, cw$site_aliases, "palmod")
  expect_false("sites" %in% names(out))
  expect_equal(unique(out$series$site_id), "GeoB3938-1")
  expect_equal(out$site_links$site_id, "GeoB3938-1")

  tables$series$site_id[1] <- "NOT_A_SITE"
  expect_error(apply_crosswalk(tables, cw$site_aliases, "palmod"), "no crosswalk entry.*NOT_A_SITE")
})

test_that("resolve_site_labels finds your cores by name, ignoring position, and never guesses", {
  sites <- build_crosswalk(dplyr::bind_rows(
    mk_site("GeoB3938_1", 12, -58, "palmod"), mk_site("MD95_2024", 50, -45, "palmod"),
    mk_site("108_658C", 20, -18, "palmod"), mk_site("GeoB1-23", 10, -20, "palmod"), mk_site("GeoB12-3", -40, 10, "palmod")
  ), "palmod")
  r <- resolve_site_labels(c("GeoB3938-1", "MD95-2024", "ODP 658B", "GeoB18549-1", "GeoB123"), sites$sites, sites$site_aliases)

  expect_equal(r$site_id[1:3], c("GeoB3938-1", "MD95-2024", "ODP-108-658"))
  expect_equal(r$match[1:3], c("name", "name", "drilling_site"))
  expect_equal(r$match[4], "none") # your aDNA core is not in the compilation
  expect_equal(r$match[5], "ambiguous") # two different sites share that key: no guess
  expect_true(is.na(r$site_id[5]))
})

test_that("a source's own spelling resolves through the aliases", {
  cw <- build_crosswalk(dplyr::bind_rows(mk_site("GIK15666_6", 35, -7, "palmod")), "palmod")
  expect_equal(resolve_site_labels("GIK15666_6", cw$sites, cw$site_aliases)$site_id, "GIK15666-6")
})

# -- golden test: your 11 aDNA cores against the real PALMOD sites (skipped unless available) ---------------

test_that("the real PALMOD sites resolve your aDNA cores the way your own matching found them", {
  dir <- Sys.getenv("PALMOD_DIR")
  tsv <- Sys.getenv("SEDACORE_GOLDEN_MATCH")
  skip_if(dir == "" || tsv == "" || !file.exists(tsv), "PALMOD_DIR / SEDACORE_GOLDEN_MATCH not set")

  palmod <- build_palmod(dir)$tables$sites
  palmod$source <- "palmod_v2"
  cw <- build_crosswalk(palmod, "palmod_v2")
  golden <- readr::read_tsv(tsv, show_col_types = FALSE)

  r <- resolve_site_labels(golden$core, cw$sites, cw$site_aliases)
  # Only GeoB3938-1 is present verbatim in PALMOD (match_type "direct"); none of the others should resolve.
  direct <- golden$match_type == "direct"
  expect_equal(sum(direct), 1)
  expect_equal(r$site_id[direct], "GeoB3938-1")
  expect_true(all(is.na(r$site_id[!direct])))

  # Every neighbour you matched is a real PALMOD site and resolves to a canonical one.
  matched <- resolve_site_labels(golden$matched_site, cw$sites, cw$site_aliases)
  expect_false(anyNA(matched$site_id))

  # ...and sits where your own search put it (distance_km from your match table).
  sites <- cw$sites
  core_pos <- golden[, c("core", "latitude", "longitude")]
  d <- haversine_km(
    golden$latitude, golden$longitude,
    sites$latitude[match(matched$site_id, sites$site_id)], sites$longitude[match(matched$site_id, sites$site_id)]
  )
  expect_equal(d, golden$distance_km, tolerance = 0.02)
})

test_that("the crosswalk at real scale: every PALMOD site gets exactly one canonical id and nothing real is merged away", {
  dir <- Sys.getenv("PALMOD_DIR")
  skip_if(dir == "" || !dir.exists(dir), "PALMOD_DIR not set")
  palmod <- build_palmod(dir)
  sites <- palmod$tables$sites
  sites$source <- "palmod_v2"
  cw <- build_crosswalk(sites, "palmod_v2")

  expect_equal(nrow(cw$site_aliases), 475)
  expect_equal(anyDuplicated(cw$sites$site_id), 0)
  # Merges are only holes of one drilling site (legs/sites), never ordinary cores.
  merged <- cw$sites[cw$sites$n_aliases > 1, ]
  expect_true(all(grepl("^(ODP|DSDP|IODP)-", merged$site_id)))
  expect_true(all(cw$site_aliases$distance_km <= 25))
})

# -- overrides ---------------------------------------------------------------------------------------------------

mk_overrides <- function(...) {
  tibble::tibble(source_a = "palmod", alias_a = "MD97_2138", source_b = "wpa", alias_b = "MD97-2138", note = "latitude sign error in PALMOD", ...)
}

test_that("a name collision is flagged with each member's position and depth", {
  sites <- dplyr::bind_rows(
    mk_site("MD97_2138", -1.25, 146.23, "palmod", 1900),
    mk_site("MD97-2138", 1.25, 146.14, "wpa", 1960)
  )
  cw <- build_crosswalk(sites, c("palmod", "wpa"))
  expect_equal(nrow(cw$sites), 2)
  expect_equal(unique(cw$review$type), "same_name_far_apart")
  d <- cw$review$detail[1]
  expect_match(d, "palmod: MD97_2138 \\(-1.250, 146.230, 1900 m\\)")
  expect_match(d, "wpa: MD97-2138 \\(1.250, 146.140, 1960 m\\)")
})

test_that("an override merges records whose positions disagree, and says so in the aliases", {
  sites <- dplyr::bind_rows(
    mk_site("MD97_2138", -1.25, 146.23, "palmod", 1900),
    mk_site("MD97-2138", 1.25, 146.14, "wpa", 1960)
  )
  cw <- build_crosswalk(sites, c("palmod", "wpa"), overrides = mk_overrides())
  expect_equal(nrow(cw$sites), 1)
  expect_equal(cw$sites$site_id, "MD97-2138")
  expect_equal(cw$sites$latitude, -1.25) # priority source's position
  expect_equal(cw$site_aliases$match_method[cw$site_aliases$source == "wpa"], "override")
  expect_equal(nrow(cw$review), 0) # no longer something to look at
})

test_that("an override for a source that is not in this build is silently skipped", {
  sites <- mk_site("MD97_2138", -1.25, 146.23, "palmod", 1900)
  expect_no_warning(cw <- build_crosswalk(sites, "palmod", overrides = mk_overrides()))
  expect_equal(nrow(cw$sites), 1)
})

test_that("an override that names a record that does not exist warns instead of failing quietly", {
  sites <- dplyr::bind_rows(mk_site("MD97_2138", -1.25, 146.23, "palmod"), mk_site("OTHER", 5, 5, "wpa"))
  expect_warning(
    build_crosswalk(sites, c("palmod", "wpa"), overrides = mk_overrides()),
    "matches no record"
  )
})

test_that("the shipped overrides file loads, is empty, and a malformed one is refused", {
  o <- load_crosswalk_overrides()
  expect_equal(names(o), c("source_a", "alias_a", "source_b", "alias_b", "note")[names(o) %in% names(o)])
  expect_true(all(c("source_a", "alias_a", "source_b", "alias_b") %in% names(o)))
  expect_equal(nrow(o), 0)
  bad <- withr::local_tempfile(fileext = ".csv")
  writeLines("a,b\n1,2", bad)
  expect_error(load_crosswalk_overrides(bad), "must have columns")
})

test_that("assemble_bundle applies overrides to the whole bundle, not just the site list", {
  p <- build_palmod(test_path("..", "fixtures", "palmod"), "GeoB3938_1_age.RDS")
  w <- read_wpa_sources(test_path("..", "fixtures", "wpa", "WA_Stable_Isotopes_2022_Sources_sample.txt"))
  # pretend the atlas spelled this core differently enough to miss the name match
  w$tables$sites$site_id[w$tables$sites$site_id == "GeoB3938-1"] <- "GB3938/1"
  w$tables$sites$label[w$tables$sites$site_id == "GB3938/1"] <- "GB3938/1"
  w$tables$indexed_records$site_id[w$tables$indexed_records$site_id == "GeoB3938-1"] <- "GB3938/1"
  w$tables$site_links$site_id[w$tables$site_links$site_id == "GeoB3938-1"] <- "GB3938/1"
  ov <- tibble::tibble(source_a = "palmod_v2", alias_a = "GeoB3938_1", source_b = "wpa_2022", alias_b = "GB3938/1", note = "checked")

  without <- assemble_bundle(list(p, w), c("palmod_v2", "wpa_2022"), overrides = NULL)
  with <- assemble_bundle(list(p, w), c("palmod_v2", "wpa_2022"), overrides = ov)
  expect_equal(nrow(without$tables$sites) - nrow(with$tables$sites), 1)
  # the atlas's record (still labelled with the atlas's own spelling) now sits on the merged site
  rec <- with$tables$indexed_records[with$tables$indexed_records$label == "GeoB3938-1", ]
  expect_equal(rec$site_id, "GeoB3938-1")
  expect_equal(without$tables$indexed_records$site_id[without$tables$indexed_records$label == "GeoB3938-1"], "GB3938/1")
})
