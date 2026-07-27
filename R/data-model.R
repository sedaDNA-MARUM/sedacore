#' Required columns and types for each sedacore schema
#' @noRd
schema_spec <- list(
  sample_registry = list(
    core_id        = "character",
    depth_cm       = "numeric",
    latitude       = "numeric",
    longitude      = "numeric",
    water_depth_m  = "numeric",
    sample_id      = "character",
    age_ka         = "numeric"
  ),
  core_registry = list(
    core_id        = "character",
    core_label     = "character",
    latitude       = "numeric",
    longitude      = "numeric",
    water_depth_m  = "numeric",
    source_program = "character",
    n_samples      = "integer",
    depth_min_cm   = "numeric",
    depth_max_cm   = "numeric",
    pangaea_doi    = "character"
  ),
  param_long = list(
    core_id          = "character",
    pangaea_doi      = "character",
    depth_cm         = "numeric",
    parameter_group  = "character",
    parameter_name   = "character",
    species          = "character",
    value            = "numeric",
    unit             = "character",
    raw_column_name  = "character",
    is_interpolated  = "logical",
    interp_gap_cm    = "numeric",
    retrieved_at     = "POSIXct"
  ),
  neighbour_registry = list(
    core_id             = "character",
    neighbour_doi       = "character",
    neighbour_core_id   = "character",
    latitude            = "numeric",
    longitude           = "numeric",
    distance_km         = "numeric",
    source_program      = "character",
    date_time           = "POSIXct",
    citation            = "character",
    included_in_export  = "logical"
  )
)

#' @noRd
check_type <- function(x, type) {
  switch(type,
    POSIXct   = inherits(x, "POSIXct"),
    integer   = is.numeric(x),
    numeric   = is.numeric(x),
    character = is.character(x),
    logical   = is.logical(x),
    inherits(x, type)
  )
}

#' @noRd
new_sedacore_tibble <- function(x, subclass, spec) {
  missing_cols <- setdiff(names(spec), names(x))
  if (length(missing_cols) > 0) {
    rlang::abort(
      c(
        sprintf("`%s` is missing required column(s):", subclass),
        rlang::set_names(missing_cols, rep("*", length(missing_cols)))
      ),
      class = "sedacore_schema_error"
    )
  }
  x <- tibble::as_tibble(x)
  x <- x[, unique(c(names(spec), setdiff(names(x), names(spec))))]
  tibble::new_tibble(x, nrow = nrow(x), class = subclass)
}

#' @noRd
validate_sedacore_tibble <- function(x, subclass, spec) {
  if (!inherits(x, subclass)) {
    rlang::abort(sprintf("`x` must have class \"%s\".", subclass))
  }
  missing_cols <- setdiff(names(spec), names(x))
  if (length(missing_cols) > 0) {
    rlang::abort(sprintf(
      "`%s` is missing required column(s): %s",
      subclass, paste(missing_cols, collapse = ", ")
    ))
  }
  bad_types <- purrr::imap_chr(spec, function(type, col) {
    if (check_type(x[[col]], type)) NA_character_ else type
  })
  bad_types <- bad_types[!is.na(bad_types)]
  if (length(bad_types) > 0) {
    rlang::abort(sprintf(
      "`%s` has columns with the wrong type: %s",
      subclass,
      paste(sprintf("%s (expected %s)", names(bad_types), bad_types), collapse = ", ")
    ))
  }
  invisible(x)
}

#' Construct a `sample_registry`
#'
#' One row per subsample: the formalized, typed shape of a raw
#' `core_metadata.tsv`-style input file.
#'
#' @param x A data frame with at least the columns `core_id`, `depth_cm`,
#'   `latitude`, `longitude`, `water_depth_m`, `sample_id`, `age_ka`.
#' @return A `sample_registry` tibble.
#' @export
new_sample_registry <- function(x) {
  new_sedacore_tibble(x, "sample_registry", schema_spec$sample_registry)
}

#' @rdname new_sample_registry
#' @export
validate_sample_registry <- function(x) {
  validate_sedacore_tibble(x, "sample_registry", schema_spec$sample_registry)
}

#' @rdname new_sample_registry
#' @export
is_sample_registry <- function(x) inherits(x, "sample_registry")

#' Construct a `core_registry`
#'
#' One row per physical core: the map-visualization and retrieval unit.
#'
#' @param x A data frame with at least the columns `core_id`, `core_label`,
#'   `latitude`, `longitude`, `water_depth_m`, `source_program`,
#'   `n_samples`, `depth_min_cm`, `depth_max_cm`, `pangaea_doi`.
#' @return A `core_registry` tibble.
#' @export
new_core_registry <- function(x) {
  new_sedacore_tibble(x, "core_registry", schema_spec$core_registry)
}

#' @rdname new_core_registry
#' @export
validate_core_registry <- function(x) {
  validate_sedacore_tibble(x, "core_registry", schema_spec$core_registry)
}

#' @rdname new_core_registry
#' @export
is_core_registry <- function(x) inherits(x, "core_registry")

#' Construct a `param_long` table
#'
#' Long-format contextual-data table: one row per (core, depth, parameter)
#' measurement. This is the universal schema that retrieval, harmonization,
#' neighbour-joining, unit conversion, interpolation, and export all read
#' and write, replacing what used to be five separate wide tables
#' (water density, TOC, isotopes, susceptibility, XRF).
#'
#' @param x A data frame with at least the columns `core_id`, `pangaea_doi`,
#'   `depth_cm`, `parameter_group`, `parameter_name`, `species`, `value`,
#'   `unit`, `raw_column_name`, `is_interpolated`, `interp_gap_cm`,
#'   `retrieved_at`.
#' @return A `param_long` tibble.
#' @export
new_param_long <- function(x) {
  new_sedacore_tibble(x, "param_long", schema_spec$param_long)
}

#' @rdname new_param_long
#' @export
validate_param_long <- function(x) {
  validate_sedacore_tibble(x, "param_long", schema_spec$param_long)
}

#' @rdname new_param_long
#' @export
is_param_long <- function(x) inherits(x, "param_long")

#' Construct a `neighbour_registry`
#'
#' Candidate neighbouring cores/datasets found within a search radius of an
#' input core, before or after user curation in the Shiny app.
#'
#' @param x A data frame with at least the columns `core_id`,
#'   `neighbour_doi`, `neighbour_core_id`, `latitude`, `longitude`,
#'   `distance_km`, `source_program`, `date_time`, `citation`,
#'   `included_in_export`.
#' @return A `neighbour_registry` tibble.
#' @export
new_neighbour_registry <- function(x) {
  new_sedacore_tibble(x, "neighbour_registry", schema_spec$neighbour_registry)
}

#' @rdname new_neighbour_registry
#' @export
validate_neighbour_registry <- function(x) {
  validate_sedacore_tibble(x, "neighbour_registry", schema_spec$neighbour_registry)
}

#' @rdname new_neighbour_registry
#' @export
is_neighbour_registry <- function(x) inherits(x, "neighbour_registry")
