fx_wpa <- function() test_path("..", "fixtures", "wpa", "WA_Stable_Isotopes_2022_Sources_sample.txt")

test_that("extract_dois handles DOIs with parentheses, punctuation and surrounding parentheses", {
  expect_equal(extract_dois("Sikes: Paleoceanography, 9, 31-45, doi:10.1029/93PA02198, 1994."), "10.1029/93PA02198")
  expect_equal(extract_dois("doi:10.1016/S0967-0637(99)00076-X, 1999."), "10.1016/S0967-0637(99)00076-X")
  expect_equal(extract_dois("(see doi:10.1029/93PA02198)"), "10.1029/93PA02198")
  expect_equal(extract_dois("(10.1016/S0031-0182(96)00136-8)."), "10.1016/S0031-0182(96)00136-8")
  expect_equal(extract_dois("two: 10.1594/PANGAEA.1 and 10.1594/PANGAEA.2"), c("10.1594/PANGAEA.1", "10.1594/PANGAEA.2"))
  expect_equal(extract_dois("no doi here"), character())
})

test_that("read_wpa_sources builds one site per core, with the atlas's own positions", {
  r <- read_wpa_sources(fx_wpa())
  s <- r$tables$sites
  expect_equal(nrow(s), 10) # 11 records on 10 cores (64PE-174P13 has two species)
  expect_equal(s$n_records[s$site_id == "64PE-174P13"], 2)
  g <- s[s$site_id == "GeoB3938-1", ]
  expect_equal(c(g$latitude, g$longitude, g$water_depth_m), c(12.258, -58.33, 1972))
  expect_equal(g$program, "GeoB")
  expect_true(all(is.na(s$basin))) # the atlas has no basin; derived later
  expect_equal(unique(s$sources), "wpa_2022")
  expect_true(all(s$position_spread_km < 1))
  expect_no_error(check_bundle_table(s, "sites"))
})

test_that("the atlas holds aDNA cores PALMOD lacks (GeoB1706-2, GeoB9064-1)", {
  s <- read_wpa_sources(fx_wpa())$tables$sites
  expect_true(all(c("GeoB1706-2", "GeoB9064-1", "GeoB3938-1") %in% s$site_id))
})

test_that("each record becomes an indexed record with its species, source DOI and radiocarbon flag", {
  r <- read_wpa_sources(fx_wpa())$tables$indexed_records
  expect_equal(nrow(r), 11)
  expect_equal(anyDuplicated(r$record_id), 0)
  expect_equal(unique(r$proxy_class), "foram isotopes")
  g <- r[r$site_id == "GeoB3938-1", ]
  expect_equal(g$species, "Cibicidoides wuellerstorfi")
  expect_equal(g$doi, "https://doi.org/10.1594/PANGAEA.730467")
  expect_true(g$has_radiocarbon)
  expect_false(r$has_radiocarbon[r$site_id == "GeoB1706-2"])
  expect_equal(r$doi[r$site_id == "12PC51"], "ftp://ftp.ncdc.noaa.gov/pub/data/paleo/paleocean/by_contributor/sikes1994/")
  # sources with no URL (laboratory, personal communication, a table in a paper) have no doi
  expect_true(all(is.na(r$doi[r$site_id %in% c("BS88-6-10B", "GeoB1028-5", "AA_GC5", "EN32-PC6")])))
  expect_true(all(nzchar(r$citation)))
  expect_no_error(check_bundle_table(r, "indexed_records"))
})

test_that("links separate publications, data and chronology, with canonical DOI URLs", {
  l <- read_wpa_sources(fx_wpa())$tables$site_links
  g <- l[l$site_id == "GeoB3938-1", ]
  expect_true("https://doi.org/10.1594/PANGAEA.730467" %in% g$url[g$link_type == "data"])
  expect_true("https://doi.org/10.5194/cp-10-843-2014" %in% g$url[g$link_type == "chronology_publication"])
  expect_true("https://doi.org/10.1016/S0967-0637(99)00076-X" %in% g$url) # parenthesised DOI intact
  # a PANGAEA DOI quoted in a reference is data, not a publication
  expect_false(any(grepl("PANGAEA", l$url[l$link_type == "primary_publication"])))
  expect_equal(anyDuplicated(l[, c("site_id", "link_type", "url")]), 0)
  expect_no_error(check_bundle_table(l, "site_links"))
})

test_that("a malformed source table is an error, not a guess", {
  bad <- withr::local_tempfile(fileext = ".txt")
  writeLines("Core\tLongitude\n12PC51\t1", bad)
  expect_error(read_wpa_sources(bad), "missing column")
  nopos <- withr::local_tempfile(fileext = ".txt")
  writeLines(c("Core\tLongitude\tLatitude\tDepth\tSpecies\tReference Stable Isotopes\tSource Stable Isotopes", "X\t\t\t100\tsp\tr\ts"), nopos)
  expect_error(read_wpa_sources(nopos), "without a numeric position")
})

test_that("the source id, version and licence are recorded for the manifest", {
  s <- read_wpa_sources(fx_wpa())$source
  expect_equal(s$id, "wpa_2022")
  expect_equal(s$doi, "10.5281/zenodo.6337519")
  expect_equal(s$license, "CC-BY-4.0")
  expect_equal(s$n_records, 11)
})

# -- joining the atlas to PALMOD ----------------------------------------------------------------------------

palmod_fixture <- function(files = c("GeoB3938_1_age.RDS", "MD95_2024_age.RDS")) {
  build_palmod(test_path("..", "fixtures", "palmod"), files)
}

test_that("the atlas and PALMOD fold onto one site per core, and the atlas adds its own cores", {
  a <- assemble_bundle(list(palmod_fixture(), read_wpa_sources(fx_wpa())), priority = c("palmod_v2", "wpa_2022"))
  t <- a$tables
  # GeoB3938_1 (PALMOD) and GeoB3938-1 (atlas) are the same core
  expect_equal(sum(t$sites$site_id == "GeoB3938-1"), 1)
  expect_equal(t$sites$sources[t$sites$site_id == "GeoB3938-1"], "palmod_v2,wpa_2022")
  expect_true(all(c("GeoB1706-2", "GeoB9064-1", "MD95-2024") %in% t$sites$site_id))
  expect_equal(nrow(t$sites), 10 + 2 - 1) # 10 atlas cores + 2 PALMOD sites - 1 shared
  # positions and depth come from PALMOD where both have the site
  expect_equal(t$sites$water_depth_m[t$sites$site_id == "GeoB3938-1"], 1972)
  expect_equal(t$sites$basin[t$sites$site_id == "GeoB3938-1"], "North Atlantic Ocean")
  # the atlas's own core still has no basin
  expect_true(is.na(t$sites$basin[t$sites$site_id == "GeoB1706-2"]))
})

test_that("index records are re-keyed onto canonical sites and carry through to the bundle tables", {
  a <- assemble_bundle(list(palmod_fixture(), read_wpa_sources(fx_wpa())), priority = c("palmod_v2", "wpa_2022"))
  r <- a$tables$indexed_records
  expect_equal(nrow(r), 11)
  expect_true(all(r$site_id %in% a$tables$sites$site_id))
  expect_equal(r$site_id[r$label == "GeoB3938-1"], "GeoB3938-1")
  expect_false("GeoB3938_1" %in% r$site_id)
  # links from both sources end up on the one site
  l <- a$tables$site_links
  expect_setequal(unique(l$source[l$site_id == "GeoB3938-1"]), c("palmod_v2", "wpa_2022"))
})

test_that("an assembled PALMOD + atlas bundle round-trips, with the index records queryable", {
  a <- assemble_bundle(list(palmod_fixture(), read_wpa_sources(fx_wpa())), priority = c("palmod_v2", "wpa_2022"))
  dir <- withr::local_tempdir()
  write_bundle(a$tables, dir, sources = a$sources, crosswalked = TRUE)
  b <- sedacore_open_bundle(dir)
  withr::defer(close_bundle(b))
  expect_equal(nrow(bundle_table(b, "indexed_records")), 11)
  expect_equal(vapply(b$manifest$sources, `[[`, "", "id"), c("palmod_v2", "wpa_2022"))
  q <- bundle_query(b, "SELECT count(DISTINCT site_id) n FROM indexed_records WHERE has_radiocarbon")
  expect_equal(q$n, 4)
})

test_that("an index record with no site yet keeps a NA site_id through the crosswalk", {
  p <- palmod_fixture("GeoB3938_1_age.RDS")
  w <- read_wpa_sources(fx_wpa())
  w$tables$indexed_records$site_id[1] <- NA_character_
  a <- assemble_bundle(list(p, w), priority = c("palmod_v2", "wpa_2022"))
  expect_equal(sum(is.na(a$tables$indexed_records$site_id)), 1)
})
