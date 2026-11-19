# Build the sedacore data bundle from the raw sources.
#
#   PALMOD_DIR=~/Downloads/PALMOD_marine_beta_2_0_0_070825 Rscript data-raw/build_bundle.R
#
# Environment:
#   PALMOD_DIR  directory of PALMOD V2 <site>_age.RDS files (required)
#   WPA_SOURCES optional: path to the atlas source table (see below)
#   BUNDLE_DIR  output directory (default: <cache dir>/bundle, where the app
#               looks first)
#
# Sources are added one section at a time; each adapter returns the same
# tables (see `bundle_schema` in R/bundle.R) so that they can simply be
# row-bound before writing. So far: PALMOD. Planned: the Mulitza World
# Paleoceanographic Atlas, the PANGAEA metadata index, the IODP holes feed.

suppressMessages(devtools::load_all(quiet = TRUE))

palmod_dir <- Sys.getenv("PALMOD_DIR", unset = path.expand("~/Downloads/PALMOD_marine_beta_2_0_0_070825"))
out_dir <- Sys.getenv("BUNDLE_DIR", unset = file.path(sedacore_cache_dir(), "bundle"))

# -- PALMOD -------------------------------------------------------------------
t0 <- Sys.time()
palmod <- build_palmod(palmod_dir, progress = function(i, n, f) {
  if (i %% 100 == 0 || i == n) message(sprintf("  PALMOD %d / %d (%s)", i, n, f))
})
message(sprintf("PALMOD read in %.0f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

if (nrow(palmod$skipped) > 0) {
  message("\nProxy columns with no entry in palmod_parameter_map.csv (not in the bundle):")
  print(dplyr::count(palmod$skipped, palmod_parameter, material, sort = TRUE))
}

# -- World Atlas of foraminiferal isotopes: source table (optional) ------------
# WPA_SOURCES = path to WA_Stable_Isotopes_2022_Sources.txt (Zenodo
# 10.5281/zenodo.6337519). Its age-controlled time series (PANGAEA.974598) are
# a separate, later adapter.
sources <- list(palmod)
priority <- "palmod_v2"
wpa_sources <- Sys.getenv("WPA_SOURCES", unset = "")
if (nzchar(wpa_sources)) {
  wpa <- read_wpa_sources(wpa_sources)
  message(sprintf("WPA source table: %d cores, %d records", nrow(wpa$tables$sites), nrow(wpa$tables$indexed_records)))
  sources <- c(sources, list(wpa))
  priority <- c(priority, "wpa_2022")
}

# -- assemble: crosswalk, canonical sites, preferred series -------------------
assembled <- assemble_bundle(sources, priority = priority)
if (nrow(assembled$review) > 0) {
  message("\nCrosswalk items for a human to look at:")
  print(assembled$review, n = 30)
}
message(sprintf(
  "\n%d source sites -> %d canonical sites; %d preferred of %d central series",
  sum(vapply(sources, function(x) nrow(x$tables$sites), 1L)), nrow(assembled$tables$sites),
  sum(assembled$tables$series_catalog$preferred), sum(assembled$tables$series_catalog$role == "central")
))

# -- write ---------------------------------------------------------------------
write_bundle(assembled$tables, out_dir,
  sources = assembled$sources, region = "global",
  crosswalked = TRUE, overwrite = TRUE
)

b <- sedacore_open_bundle(out_dir)
on.exit(close_bundle(b), add = TRUE)
print(b)
message(sprintf("\nBundle size on disk: %.1f MB", sum(file.size(list.files(out_dir, full.names = TRUE))) / 1024^2))
