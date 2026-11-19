# Build the sedacore data bundle from the raw sources.
#
#   PALMOD_DIR=~/Downloads/PALMOD_marine_beta_2_0_0_070825 Rscript data-raw/build_bundle.R
#
# Environment:
#   PALMOD_DIR  directory of PALMOD V2 <site>_age.RDS files (required)
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

tables <- palmod$tables
sources <- list(palmod$source)

# -- write ---------------------------------------------------------------------
write_bundle(tables, out_dir, sources = sources, region = "global (all PALMOD V2 sites)", overwrite = TRUE)

b <- sedacore_open_bundle(out_dir)
on.exit(close_bundle(b), add = TRUE)
print(b)
message(sprintf("\nBundle size on disk: %.1f MB", sum(file.size(list.files(out_dir, full.names = TRUE))) / 1024^2))
