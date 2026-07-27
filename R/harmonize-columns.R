#' Load the column-synonym lookup table
#'
#' @param path Path to `column_synonyms.csv`.
#' @return A tibble: `raw_pattern, match_type, parameter_group,
#'   parameter_name, unit_raw, species_capture_group, notes`.
#' @export
load_synonym_map <- function(
    path = system.file("extdata", "column_synonyms.csv", package = "sedacore")) {
  readr::read_csv(
    path,
    show_col_types = FALSE,
    col_types = readr::cols(.default = "c")
  )
}

#' @noRd
match_column_to_synonym <- function(col, synonym_map) {
  for (i in seq_len(nrow(synonym_map))) {
    row <- synonym_map[i, ]
    is_hit <- if (identical(row$match_type, "exact")) {
      identical(col, row$raw_pattern)
    } else {
      isTRUE(stringr::str_detect(col, row$raw_pattern))
    }
    if (!isTRUE(is_hit)) next

    species <- NA_character_
    if (!is.na(row$species_capture_group)) {
      gi <- as.integer(row$species_capture_group) + 1L
      m <- stringr::str_match(col, row$raw_pattern)
      if (ncol(m) >= gi) species <- m[1, gi]
    }

    return(list(
      parameter_group = row$parameter_group,
      parameter_name = row$parameter_name,
      unit_raw = row$unit_raw,
      species = species
    ))
  }
  NULL
}

#' Best-effort parameter info for a column not in the synonym table
#'
#' PANGAEA has thousands of distinct measurement types; a synonym table
#' can never be fully exhaustive, and new real-world datasets keep
#' surfacing columns it doesn't yet cover (trace element/nutrient
#' pore-water chemistry -- Fe, PO4, Si(OH)4 -- was the case that prompted
#' this). Rather than drop genuinely numeric data just because it isn't
#' yet catalogued, parse a reasonable `parameter_name`/`unit` straight out
#' of the raw column header (splitting off a trailing `[unit]`, if any)
#' and tag it `parameter_group = "uncategorized"` so it's still visible
#' and exportable -- add a proper synonym-table row later once it's
#' apparent the same column recurs.
#'
#' @param col A raw column name, e.g. `"Fe [µmol/l]"`.
#' @return A list shaped like `match_column_to_synonym()`'s return value.
#' @noRd
fallback_column_info <- function(col) {
  unit_match <- stringr::str_match(col, "\\[([^\\]]+)\\]\\s*$")
  unit <- if (!is.na(unit_match[1, 2])) unit_match[1, 2] else NA_character_
  name <- stringr::str_trim(stringr::str_remove(col, "\\[[^\\]]+\\]\\s*$"))
  if (!nzchar(name)) name <- col
  list(parameter_group = "uncategorized", parameter_name = name, unit_raw = unit, species = NA_character_)
}

#' Harmonize a raw PANGAEA parameter table into `param_long`
#'
#' Pivots a raw wide PANGAEA data table (as returned in a row of
#' [pg_fetch_dataset()]'s `data` list-column) into the long
#' [new_param_long()] schema. Every numeric column is kept -- the whole
#' point of this package is to gather as much contextual data as possible,
#' so nothing measurable is silently dropped just because it isn't yet in
#' `column_synonyms.csv`. Columns matching `synonym_map` (see
#' [load_synonym_map()]) get their canonical `parameter_group`/
#' `parameter_name`/`unit`/`species`; anything else falls back to
#' `fallback_column_info()` (raw column name/unit, `parameter_group =
#' "uncategorized"`) -- add a synonym-table row later once a fallback
#' column turns out to recur across datasets. A column is only ever
#' excluded if it's genuinely non-numeric (e.g. free-text comment fields
#' coerce entirely to `NA` and carry no measurement to represent).
#'
#' @param raw_data A raw wide data frame from a PANGAEA dataset.
#' @param core_id The core this dataset belongs to.
#' @param doi The dataset's PANGAEA DOI (for provenance).
#' @param depth_col Name of the depth column in `raw_data`, before
#'   harmonization (PANGAEA convention: metres).
#' @param depth_unit_to_cm Multiplier converting `depth_col` to
#'   centimetres (default `100`, i.e. metres to cm).
#' @param synonym_map A synonym table, from [load_synonym_map()].
#' @param retrieved_at Timestamp to stamp on every row.
#' @param unmatched One of `"warn"` (default), `"silent"`, or `"error"`:
#'   what to do about columns not in `synonym_map` that fall back to
#'   `fallback_column_info()` -- they're included either way, this only
#'   controls whether/how loudly it's flagged (useful for noticing when a
#'   fallback column recurs often enough to deserve a proper synonym-table
#'   entry).
#' @return A `param_long` tibble. Column names that used the fallback
#'   parser are attached as `attr(x, "fallback_columns")`; genuinely
#'   non-numeric columns that were excluded are attached as
#'   `attr(x, "dropped_columns")`.
#' @export
harmonize_columns <- function(raw_data,
                               core_id,
                               doi = NA_character_,
                               depth_col = "Depth sed [m]",
                               depth_unit_to_cm = 100,
                               synonym_map = load_synonym_map(),
                               retrieved_at = Sys.time(),
                               unmatched = c("warn", "silent", "error")) {
  unmatched <- match.arg(unmatched)

  if (!depth_col %in% names(raw_data)) {
    rlang::abort(sprintf("`depth_col` (\"%s\") not found in `raw_data`.", depth_col))
  }
  # round() to clear floating-point noise from the unit conversion (e.g.
  # 0.28 * 100 == 28.000000000000004, not 28) -- left unrounded, this both
  # pollutes exported TSVs with spurious digits and silently breaks exact-
  # depth matches in interpolate_depth_series() against clean sample
  # depths. 6 dp is far finer than any real sediment sampling resolution.
  depth_cm <- round(as.numeric(raw_data[[depth_col]]) * depth_unit_to_cm, 6)

  candidate_cols <- setdiff(names(raw_data), depth_col)
  matched <- stats::setNames(
    lapply(candidate_cols, match_column_to_synonym, synonym_map = synonym_map),
    candidate_cols
  )
  is_matched <- !vapply(matched, is.null, logical(1))
  fallback_cols <- candidate_cols[!is_matched]

  if (length(fallback_cols) > 0) {
    msg <- sprintf(
      "Column(s) in dataset %s not in column_synonyms.csv, imported using their raw name: %s",
      doi %||% "?", paste(fallback_cols, collapse = ", ")
    )
    if (unmatched == "warn") rlang::warn(msg, class = "sedacore_unmatched_columns")
    if (unmatched == "error") rlang::abort(msg, class = "sedacore_unmatched_columns")
  }
  for (col in fallback_cols) {
    matched[[col]] <- fallback_column_info(col)
  }

  dropped_cols <- character()
  rows <- purrr::imap_dfr(matched, function(m, col) {
    vals <- suppressWarnings(as.numeric(raw_data[[col]]))
    if (all(is.na(vals))) {
      dropped_cols <<- c(dropped_cols, col)
      return(tibble::tibble())
    }
    tibble::tibble(
      core_id = core_id,
      pangaea_doi = doi,
      depth_cm = depth_cm,
      parameter_group = m$parameter_group,
      parameter_name = m$parameter_name,
      species = m$species,
      value = vals,
      unit = m$unit_raw,
      raw_column_name = col,
      is_interpolated = FALSE,
      interp_gap_cm = NA_real_,
      retrieved_at = retrieved_at
    )
  })

  if (nrow(rows) == 0) {
    rows <- tibble::tibble(
      core_id = character(), pangaea_doi = character(), depth_cm = numeric(),
      parameter_group = character(), parameter_name = character(),
      species = character(), value = numeric(), unit = character(),
      raw_column_name = character(), is_interpolated = logical(),
      interp_gap_cm = numeric(), retrieved_at = as.POSIXct(character())
    )
  }

  out <- new_param_long(rows)
  attr(out, "fallback_columns") <- setdiff(fallback_cols, dropped_cols)
  attr(out, "dropped_columns") <- dropped_cols
  out
}
