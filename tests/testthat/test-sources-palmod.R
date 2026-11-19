fx <- function(f) test_path("..", "fixtures", "palmod", f)

# Read a fixture, change it, write it to a temp file: how the edge cases below
# are produced without needing more real data.
modified_fixture <- function(f, fn) {
  x <- readRDS(fx(f))
  x <- fn(x)
  path <- withr::local_tempfile(fileext = ".RDS", .local_envir = parent.frame())
  saveRDS(x, path)
  path
}

# -- name resolution -----------------------------------------------------------------

test_that("palmod_resolve_name picks the most specific map row", {
  map <- load_palmod_parameter_map()
  r <- palmod_resolve_name
  expect_equal(r("surface.temp", "MgCa", NA, map), "SST_MgCa")
  expect_equal(r("surface.temp", "UK37", NA, map), "SST_UK37")
  expect_equal(r("surface.temp", "never heard of it", NA, map), "SST_other") # catch-all
  expect_equal(r("surface.temp", NA, "degC", map), "SST_other")
  expect_equal(r("IRD", NA, "percent", map), "IRD_percent")
  expect_equal(r("IRD", NA, "no_g", map), "IRD_per_g")
  expect_equal(r("IRD", NA, "n", map), "IRD_count")
  expect_equal(r("IRD", NA, "n_cm3", map), "IRD_per_cm3")
  expect_equal(r("IRD", NA, NA, map), "IRD_other") # catch-all, but a listed unit always wins
  expect_equal(r("IRD", NA, "percent", map), "IRD_percent")
  expect_equal(r("benthic.d18O", "benthic foraminifera calcite", "permil VPDB", map), "d18O_benthic")
  expect_true(is.na(r("unknown.parameter", NA, NA, map)))
})

test_that("planktonic and benthic isotopes get distinct names (they must never be pooled)", {
  map <- load_palmod_parameter_map()
  expect_equal(
    palmod_resolve_name(c("benthic.d18O", "planktonic.d18O"), c(NA, NA), c(NA, NA), map),
    c("d18O_benthic", "d18O_planktonic")
  )
})

test_that("every canonical name the map can produce is in the proxy vocabulary", {
  map <- load_palmod_parameter_map()
  vocab <- load_proxy_vocabulary()
  expect_setequal(setdiff(map$parameter_name, vocab$parameter_name), character())
})

test_that("palmod_clean_type fixes the source's typos", {
  expect_equal(palmod_clean_type(c("measured", "meaasured", "measurd", "inferred", NA)),
    c("measured", "measured", "measured", "inferred", NA))
})

# -- a real site: GeoB3938_1 (an aDNA core present verbatim in PALMOD) ----------------------

test_that("read_palmod_site reads site metadata, with water depth as a positive number", {
  r <- read_palmod_site(fx("GeoB3938_1_age.RDS"))
  expect_equal(r$site$site_id, "GeoB3938_1")
  expect_equal(r$site$water_depth_m, 1972) # PALMOD stores -1972
  expect_equal(r$site$basin, "North Atlantic Ocean")
  expect_equal(r$site$program, "GeoB")
  expect_equal(r$site$latitude, 12.25833, tolerance = 1e-5)
})

test_that("read_palmod_site maps proxies to canonical names and classes", {
  r <- read_palmod_site(fx("GeoB3938_1_age.RDS"))
  cat <- r$series_catalog
  expect_setequal(cat$parameter_name, c("d18O_benthic", "TOC", "CaCO3"))
  expect_equal(cat$proxy_class[cat$parameter_name == "d18O_benthic"], "d18O benthic")
  expect_equal(cat$species[cat$parameter_name == "d18O_benthic"], "C. wuellerstorfi")
  expect_equal(cat$unit[cat$parameter_name == "d18O_benthic"], "permil PDB") # raw, not converted
  expect_true(all(cat$role == "central"))
  expect_equal(nrow(r$skipped), 0)
})

test_that("series rows carry the recalibrated Bacon age, never a missing age or value", {
  r <- read_palmod_site(fx("GeoB3938_1_age.RDS"))
  s <- r$series
  expect_false(anyNA(s$age_ka))
  expect_false(anyNA(s$value))
  expect_true(all(s$age_lo_ka <= s$age_ka & s$age_ka <= s$age_hi_ka))
  # ages are exactly PALMOD's median age at the row's depth
  x <- readRDS(fx("GeoB3938_1_age.RDS"))
  toc <- s[s$parameter_name == "TOC", ]
  expect_equal(sort(toc$age_ka), sort(x$data$medianAge_kaBP[!is.na(x$data$medianAge_kaBP) & !is.na(x$data$TOC_percent)]))
  # depth_merged is metres in PALMOD; the bundle uses centimetres
  expect_equal(max(toc$depth_cm), max(x$data$depth_merged[!is.na(x$data$medianAge_kaBP) & !is.na(x$data$TOC_percent)]) * 100)
})

test_that("catalog accounts for every observation: aged + dropped = total", {
  for (f in c("GeoB3938_1_age.RDS", "MD01_2412_age.RDS", "MD95_2024_age.RDS")) {
    cat <- read_palmod_site(fx(f))$series_catalog
    expect_equal(cat$n_obs_aged + cat$n_dropped_no_age, cat$n_obs, info = f)
  }
})

test_that("series rows per series equal the catalog's aged count", {
  r <- read_palmod_site(fx("MD95_2024_age.RDS"))
  n <- table(r$series$series_id)
  cat <- r$series_catalog
  expect_equal(as.integer(n[cat$series_id]), cat$n_obs_aged)
})

# -- MD95_2024: most Bacon ages missing (the GeoB18549 neighbour) ----------------------------------

test_that("MD95_2024 loses most bulk rows to missing Bacon ages, and says so", {
  cat <- read_palmod_site(fx("MD95_2024_age.RDS"))$series_catalog
  toc <- cat[cat$parameter_name == "TOC", ]
  expect_gt(toc$n_dropped_no_age, 5 * toc$n_obs_aged)
  # PALMOD's recalibrated ages for this core stop around 12 ka: nothing here
  # covers the LGM, even though the underlying data do.
  expect_lt(max(cat$age_max_ka, na.rm = TRUE), 12.5)
})

# -- MD01_2412: recalculated temperature with bounds ---------------------------------------------------

test_that("an adjacent .lo/central/.hi triple is folded into value_lo/value_hi", {
  r <- read_palmod_site(fx("MD01_2412_age.RDS"))
  cat <- r$series_catalog
  expect_false(any(cat$role %in% c("lower", "upper"))) # folded, not separate series
  sst <- cat[cat$parameter_name == "SST_UK37", ]
  expect_equal(nrow(sst), 2) # the original and the recalculated series
  expect_equal(sst$has_band[order(sst$is_recalc)], c(FALSE, TRUE))

  s <- r$series[r$series$is_recalc & r$series$parameter_name == "SST_UK37", ]
  expect_false(anyNA(s$value_lo))
  expect_true(all(s$value_lo <= s$value & s$value <= s$value_hi))
  orig <- r$series[!r$series$is_recalc & r$series$parameter_name == "SST_UK37", ]
  expect_true(all(is.na(orig$value_lo)))
})

test_that("a triple whose bounds do not bracket the value is NOT folded (it stays separate)", {
  path <- modified_fixture("MD01_2412_age.RDS", function(x) {
    j <- which(x$meta$Parameter == "surface.temp.lo")
    k <- which(x$meta$Parameter == "surface.temp.hi")
    tmp <- x$data[[j]]; x$data[[j]] <- x$data[[k]]; x$data[[k]] <- tmp # swap lo/hi
    x
  })
  r <- read_palmod_site(path)
  expect_true(all(c("lower", "upper") %in% r$series_catalog$role))
  expect_false(any(r$series_catalog$has_band))
})

test_that("a lone .lo/.hi (no matching neighbours) is kept as its own series, not lost", {
  path <- modified_fixture("MD01_2412_age.RDS", function(x) {
    x$meta$Material[x$meta$Parameter == "surface.temp.lo"] <- "TEX86" # breaks the Material match
    x
  })
  r <- read_palmod_site(path)
  expect_true("lower" %in% r$series_catalog$role)
})

# -- skips, typos, structural errors -----------------------------------------------------------------------

test_that("an unmapped proxy is reported in `skipped`, not silently dropped", {
  path <- modified_fixture("GeoB3938_1_age.RDS", function(x) {
    x$meta$Parameter[x$meta$Parameter == "CaCO3"] <- "mystery.proxy"
    x
  })
  r <- read_palmod_site(path)
  expect_equal(r$skipped$palmod_parameter, "mystery.proxy")
  expect_false("CaCO3" %in% r$series_catalog$parameter_name)
})

test_that("ParameterType typos are normalised in the catalog", {
  path <- modified_fixture("GeoB3938_1_age.RDS", function(x) {
    x$meta$ParameterType[x$meta$Parameter == "TOC"] <- "meaasured"
    x
  })
  cat <- read_palmod_site(path)$series_catalog
  expect_equal(cat$parameter_type[cat$parameter_name == "TOC"], "measured")
})

test_that("a meta/data misalignment is an error, not a guess", {
  path <- modified_fixture("GeoB3938_1_age.RDS", function(x) {
    x$meta <- x$meta[-1, ]
    x
  })
  expect_error(read_palmod_site(path), "expected to align")
})

test_that("a site with no depth column still reads (age-native), with NA depth", {
  path <- modified_fixture("GeoB3938_1_age.RDS", function(x) {
    x$meta$Parameter[x$meta$Parameter == "depth_merged"] <- "label"
    x
  })
  r <- read_palmod_site(path)
  expect_true(all(is.na(r$series$depth_cm)))
  expect_gt(nrow(r$series), 0)
})

# -- age control ------------------------------------------------------------------------------------------------------

test_that("age_control carries the dated levels with typed columns", {
  ac <- read_palmod_site(fx("MD01_2412_age.RDS"))$age_control
  expect_equal(nrow(ac), 16)
  # radiocarbon dates plus a tephra level (which has no calendar-age bounds)
  expect_setequal(unique(ac$control_type), c("14C", "tephra"))
  expect_type(ac$rejected, "logical")
  expect_type(ac$depth_cm, "double")
  dated <- ac[!is.na(ac$cal_age_min_ka) & !is.na(ac$cal_age_max_ka), ]
  expect_gt(nrow(dated), 0)
  expect_true(all(dated$cal_age_min_ka <= dated$cal_age_ka & dated$cal_age_ka <= dated$cal_age_max_ka))
})

# -- links ---------------------------------------------------------------------------------------------------------------

test_that("links are canonical DOI URLs, deduplicated across resolver forms", {
  l <- read_palmod_site(fx("MD95_2024_age.RDS"))$site_links
  expect_false(any(grepl("doi.pangaea.de", l$url)))
  expect_equal(anyDuplicated(l[, c("link_type", "url")]), 0)
  expect_true(all(c("primary_publication", "data", "chronology_publication") %in% l$link_type))
  expect_true("https://doi.org/10.1594/PANGAEA.778963" %in% l$url) # the dataset used earlier in the project
})

# -- building the bundle ---------------------------------------------------------------------------------------------------

test_that("build_palmod combines sites and every table passes the bundle schema", {
  files <- c("GeoB3938_1_age.RDS", "MD01_2412_age.RDS", "MD95_2024_age.RDS")
  b <- build_palmod(test_path("..", "fixtures", "palmod"), files)
  expect_equal(nrow(b$tables$sites), 3)
  expect_equal(b$source$id, "palmod_v2")
  expect_equal(b$source$version, "2.0.0")
  expect_equal(b$source$doi, "10.1594/PANGAEA.984602")
  expect_equal(b$source$n_sites, 3)
  for (nm in names(b$tables)) expect_no_error(check_bundle_table(b$tables[[nm]], nm))
  expect_equal(nrow(b$skipped), 0)
})

test_that("build_palmod refuses duplicate sites and an empty directory", {
  expect_error(
    build_palmod(test_path("..", "fixtures", "palmod"), c("GeoB3938_1_age.RDS", "GeoB3938_1_age.RDS")),
    "Duplicate PALMOD site ids"
  )
  expect_error(build_palmod(withr::local_tempdir()), "No PALMOD .RDS files")
})

test_that("a built PALMOD bundle round-trips through the bundle layer", {
  files <- c("GeoB3938_1_age.RDS", "MD95_2024_age.RDS")
  built <- build_palmod(test_path("..", "fixtures", "palmod"), files)
  dir <- withr::local_tempdir()
  write_bundle(built$tables, dir, sources = list(built$source), region = "global")
  b <- sedacore_open_bundle(dir)
  withr::defer(close_bundle(b))

  expect_equal(nrow(bundle_table(b, "sites")), 2)
  s <- bundle_series(b, site_id = "GeoB3938_1", parameter_name = "d18O_benthic")
  expect_equal(nrow(s), built$tables$series_catalog$n_obs_aged[built$tables$series_catalog$site_id == "GeoB3938_1" &
    built$tables$series_catalog$parameter_name == "d18O_benthic"])
  expect_true(all(diff(s$age_ka) >= 0)) # age order
  expect_equal(b$manifest$sources[[1]]$license, "CC-BY-4.0")
})

# -- reconciliation against the full download (skipped unless available) ------------------------------------------------------

test_that("the full PALMOD V2 download reconciles with what was profiled by hand", {
  dir <- Sys.getenv("PALMOD_DIR")
  skip_if(dir == "" || !dir.exists(dir), "PALMOD_DIR not set (the 1.2 GB download is not part of the repo)")

  b <- build_palmod(dir)
  sites <- b$tables$sites
  cat <- b$tables$series_catalog
  expect_equal(nrow(sites), 475)
  expect_equal(sum(sites$basin == "North Atlantic Ocean"), 154)
  expect_equal(sum(sites$basin == "South Atlantic Ocean"), 61)

  sites_with <- function(nm) length(unique(cat$site_id[cat$parameter_name %in% nm & cat$role == "central"]))
  expect_equal(sites_with("d18O_benthic"), 345)
  expect_equal(sites_with("d18O_planktonic"), 310)
  expect_equal(sites_with(c("IRD_count", "IRD_per_g", "IRD_percent", "IRD_per_cm3", "IRD_other")), 37)
  # Every proxy column maps to something: nothing is silently unmapped
  expect_equal(nrow(b$skipped), 0)
  for (nm in names(b$tables)) check_bundle_table(b$tables[[nm]], nm)
})
