# Escape hatch for one-off PANGAEA datasets whose raw columns can't be
# handled generically by inst/extdata/column_synonyms.csv (e.g. a file that
# mixes multiple species' measurements in one set of columns, needing a
# coalesce() before harmonize_columns() can pivot it correctly). Keyed by
# DOI. Each entry is a function(raw_data) -> raw_data, applied by
# apply_dataset_override() (R/pangaea-parse.R) before harmonize_columns().
#
# Example (not real, illustrative only):
# dataset_overrides <- list(
#   "10.1594/PANGAEA.XXXXXX" = function(raw_data) {
#     dplyr::mutate(
#       raw_data,
#       species = ifelse(!is.na(`Cibicides sp. d18O [permil PDB]`), "Cibicides sp.", "Melonis sp.")
#     )
#   }
# )

dataset_overrides <- list()
