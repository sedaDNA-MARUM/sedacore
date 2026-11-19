# De-duplicating the same measurement across sources.
#
# PALMOD and the Mulitza atlas both compile G. bulloides d18O for the same
# core. Both are kept (they differ in age model, and one may have rows the
# other lacks), but coverage must count the record once. Each group of series
# that measure the same thing at the same site gets exactly one `preferred`
# series; everything downstream that counts or plots uses only those.

#' Canonical form of a foraminiferal species name
#'
#' Sources write the same species as `Cibicidoides wuellerstorfi`,
#' `C. wuellerstorfi` or `Planulina wuellerstorfi`. This lower-cases, drops
#' qualifiers in brackets, abbreviates the genus to its initial, and applies a
#' small table of known synonyms.
#'
#' @param x Character vector of species names.
#' @return Character vector (`NA` stays `NA`; empty strings become `NA`).
#' @keywords internal
#' @export
normalize_species <- function(x) {
  y <- tolower(trimws(gsub("\\(.*?\\)", "", x)))
  y <- gsub("\\s+", " ", y)
  # genus synonyms first, so they abbreviate to the accepted genus
  y <- sub("^planulina wuellerstorfi", "cibicidoides wuellerstorfi", y)
  y <- sub("^cibicides wuellerstorfi", "cibicidoides wuellerstorfi", y)
  y <- sub("^globigerina bulloides", "globigerina bulloides", y)
  # abbreviate "Genus species ..." -> "g. species ..."; leave "g. species" alone
  y <- sub("^([a-z])[a-z]+ ([a-z])", "\\1. \\2", y)
  y[is.na(x) | !nzchar(y)] <- NA_character_
  y
}

#' Choose one preferred series per (site, proxy, species)
#'
#' Among `central` series with the same `site_id`, `parameter_name` and
#' normalised species, the winner is, in order: the source listed first in
#' `priority`; a recalculated series over the original (PALMOD ships both for
#' temperature); the one with more age-controlled observations; then
#' `series_id`, so the choice is deterministic. Non-central series (bounds,
#' errors) are never preferred.
#'
#' @param catalog A `series_catalog` (needs `source`, `is_recalc`,
#'   `n_obs_aged`, `role`, `species`; `series_id` must be unique).
#' @param priority Source ids, highest priority first.
#' @return `catalog` with a logical `preferred` column.
#' @keywords internal
#' @export
mark_preferred_series <- function(catalog, priority) {
  if (anyDuplicated(catalog$series_id)) {
    rlang::abort("`series_id` must be unique across sources before choosing preferred series.")
  }
  catalog$.species <- normalize_species(catalog$species)
  ord <- order(
    catalog$site_id, catalog$parameter_name, catalog$.species,
    match(catalog$source, priority), !catalog$is_recalc, -catalog$n_obs_aged, catalog$series_id,
    na.last = TRUE
  )
  key <- paste(catalog$site_id, catalog$parameter_name, ifelse(is.na(catalog$.species), "", catalog$.species), sep = "\r")
  first_in_group <- !duplicated(key[ord])
  winner <- rep(FALSE, nrow(catalog))
  winner[ord[first_in_group]] <- TRUE
  catalog$preferred <- winner & catalog$role == "central"
  # a group whose best candidate is non-central has no central winner: pick
  # among the central ones only
  has_winner <- stats::ave(catalog$preferred, key, FUN = any)
  need <- !has_winner & catalog$role == "central"
  if (any(need)) {
    sub <- which(need)
    o2 <- sub[order(
      key[sub], match(catalog$source[sub], priority), !catalog$is_recalc[sub],
      -catalog$n_obs_aged[sub], catalog$series_id[sub]
    )]
    catalog$preferred[o2[!duplicated(key[o2])]] <- TRUE
  }
  catalog$.species <- NULL
  catalog
}

#' Carry `preferred` from the catalog onto the observations
#'
#' @param series A `series` table.
#' @param catalog A catalog from [mark_preferred_series()].
#' @return `series` with `preferred`.
#' @keywords internal
#' @export
add_preferred_to_series <- function(series, catalog) {
  series$preferred <- catalog$preferred[match(series$series_id, catalog$series_id)]
  if (anyNA(series$preferred)) {
    rlang::abort("Some `series` rows belong to a series that is not in the catalog.")
  }
  series
}
