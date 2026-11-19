# Cut small test fixtures from three real PALMOD V2 site files.
#
# The adapter (R/sources-palmod.R) never reads the Bacon age ensembles, which
# are almost the whole size of a site file (1000 columns each), so they are
# dropped here. Everything the adapter does read is left untouched.
#
#   PALMOD_DIR=~/Downloads/PALMOD_marine_beta_2_0_0_070825 Rscript data-raw/make_palmod_fixtures.R
#
# Data: PALMOD 130k marine palaeoclimate synthesis V2.0.0, Jonkers et al.
# (2026), https://doi.org/10.1594/PANGAEA.984602, CC BY 4.0.

dir <- Sys.getenv("PALMOD_DIR", unset = path.expand("~/Downloads/PALMOD_marine_beta_2_0_0_070825"))
out <- file.path("tests", "fixtures", "palmod")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

for (f in c("GeoB3938_1_age.RDS", "MD01_2412_age.RDS", "MD95_2024_age.RDS")) {
  x <- readRDS(file.path(dir, f))
  x$AgeEnsemblesBACON <- NULL
  x$SurfaceTempEnsembles <- NULL
  saveRDS(x, file.path(out, f), compress = "xz")
  cat(sprintf("%-22s %7.1f KB\n", f, file.size(file.path(out, f)) / 1024))
}
