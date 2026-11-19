# The data bundle: a versioned directory of parquet tables plus a manifest.
#
# The atlas cannot be built by fetching thousands of datasets live, so the
# compiled sources (PALMOD, the Mulitza atlas, the PANGAEA/IODP indexes)
# are turned *offline* into a small set of tables that the app only reads.
# A bundle is just
#
#   <dir>/manifest.json      what is in it, from which sources/versions, under
#                            which licences, with checksums
#   <dir>/<table>.parquet    one file per table (see `bundle_schema`)
#
# DuckDB reads the parquet files in place (no full load into memory) and can
# aggregate them in SQL, which is what the gap-map queries need. Nothing here
# is specific to the Atlantic: the region is a filter at query time.

# Required columns and types per table. Extra columns are allowed (and are
# how a source adds detail without a schema change), mirroring `schema_spec`.
#' @noRd
bundle_schema <- list(
  sites = c(
    site_id = "character", label = "character", latitude = "numeric", longitude = "numeric",
    water_depth_m = "numeric", basin = "character", program = "character", sources = "character"
  ),
  series = c(
    site_id = "character", series_id = "character", source = "character",
    parameter_name = "character", proxy_class = "character", species = "character",
    material = "character", role = "character", is_recalc = "logical",
    depth_cm = "numeric", age_ka = "numeric", age_lo_ka = "numeric", age_hi_ka = "numeric",
    value = "numeric", value_lo = "numeric", value_hi = "numeric", unit = "character"
  ),
  series_catalog = c(
    site_id = "character", series_id = "character", source = "character",
    parameter_name = "character", proxy_class = "character", species = "character",
    material = "character", parameter_type = "character", role = "character",
    is_recalc = "logical", unit = "character", n_obs = "integer", n_obs_aged = "integer",
    n_dropped_no_age = "integer", age_min_ka = "numeric", age_max_ka = "numeric"
  ),
  age_control = c(
    site_id = "character", source = "character", control_type = "character",
    depth_cm = "numeric", age_14c_ka = "numeric", cal_age_ka = "numeric",
    cal_age_min_ka = "numeric", cal_age_max_ka = "numeric", rejected = "logical"
  ),
  site_links = c(
    site_id = "character", link_type = "character", label = "character",
    url = "character", source = "character"
  )
)

#' Where bundles live by default
#'
#' The one place that decides the cache location (derived from the package
#' name, so a rename needs no edit here).
#'
#' @return A directory path (not created).
#' @keywords internal
#' @export
sedacore_cache_dir <- function() {
  tools::R_user_dir(utils::packageName() %||% "sedacore", "cache")
}

#' Locate the bundle to use
#'
#' In order: an explicit `path`, the `sedacore.bundle` option, the
#' `SEDACORE_BUNDLE` environment variable, then `bundle/` under
#' [sedacore_cache_dir()].
#'
#' @param path Optional explicit bundle directory.
#' @return A directory path that contains a `manifest.json`.
#' @keywords internal
#' @export
sedacore_bundle_path <- function(path = NULL) {
  candidates <- c(
    path,
    getOption("sedacore.bundle"),
    Sys.getenv("SEDACORE_BUNDLE", unset = NA_character_),
    file.path(sedacore_cache_dir(), "bundle")
  )
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  found <- candidates[file.exists(file.path(candidates, "manifest.json"))]
  if (length(found) == 0) {
    rlang::abort(c(
      "No sedacore data bundle found.",
      i = sprintf("Looked in: %s", paste(candidates, collapse = ", ")),
      i = "Build one with the scripts in data-raw/, or point option `sedacore.bundle` or env var `SEDACORE_BUNDLE` at an existing bundle directory."
    ))
  }
  normalizePath(found[1])
}

#' @noRd
check_bundle_table <- function(x, name) {
  spec <- bundle_schema[[name]]
  if (is.null(spec)) {
    rlang::abort(sprintf("Unknown bundle table \"%s\" (known: %s).", name, paste(names(bundle_schema), collapse = ", ")))
  }
  missing_cols <- setdiff(names(spec), names(x))
  if (length(missing_cols) > 0) {
    rlang::abort(sprintf("Bundle table `%s` is missing required column(s): %s", name, paste(missing_cols, collapse = ", ")))
  }
  bad <- names(spec)[!vapply(names(spec), function(col) check_type(x[[col]], spec[[col]]), logical(1))]
  if (length(bad) > 0) {
    rlang::abort(sprintf(
      "Bundle table `%s` has columns of the wrong type: %s",
      name, paste(sprintf("%s (expected %s)", bad, spec[bad]), collapse = ", ")
    ))
  }
  invisible(x)
}

#' Write a data bundle
#'
#' Validates each table against the bundle schema, writes it as a zstd-
#' compressed parquet file, and writes `manifest.json` recording the format
#' version, what was built and when, the sources (with versions, DOIs and
#' licences) and each file's row count and checksum.
#'
#' @param tables A named list of data frames; names must be bundle tables
#'   (`sites`, `series`, `series_catalog`, `age_control`, `site_links`).
#' @param dir Output directory (created if needed).
#' @param sources A list of source descriptions, each a list with at least
#'   `id`, `name`, `version`, `doi`, `license` and `citation`.
#' @param region Free-text note on what region/scope the bundle covers.
#' @param overwrite Replace an existing bundle in `dir`.
#' @return `dir`, invisibly.
#' @keywords internal
#' @export
write_bundle <- function(tables, dir, sources = list(), region = "global", overwrite = FALSE) {
  if (is.null(names(tables)) || any(!nzchar(names(tables)))) {
    rlang::abort("`tables` must be a named list.")
  }
  for (nm in names(tables)) check_bundle_table(tables[[nm]], nm)
  if (file.exists(file.path(dir, "manifest.json")) && !overwrite) {
    rlang::abort(sprintf("A bundle already exists in \"%s\"; use `overwrite = TRUE` to replace it.", dir))
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  files <- lapply(names(tables), function(nm) {
    path <- file.path(dir, paste0(nm, ".parquet"))
    duckdb::duckdb_register(con, "bundle_tmp", tables[[nm]])
    DBI::dbExecute(con, sprintf(
      "COPY (SELECT * FROM bundle_tmp) TO %s (FORMAT parquet, COMPRESSION zstd)",
      DBI::dbQuoteString(con, path)
    ))
    duckdb::duckdb_unregister(con, "bundle_tmp")
    list(
      file = basename(path), rows = nrow(tables[[nm]]),
      md5 = unname(tools::md5sum(path)), bytes = file.size(path)
    )
  })
  names(files) <- names(tables)

  manifest <- list(
    bundle_format = 1L,
    created = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    package = utils::packageName() %||% "sedacore",
    package_version = as.character(utils::packageVersion(utils::packageName() %||% "sedacore")),
    region = region,
    tables = files,
    sources = sources
  )
  jsonlite::write_json(manifest, file.path(dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
  invisible(dir)
}

#' Open a data bundle
#'
#' Returns a handle whose tables are DuckDB views over the parquet files, so
#' large tables are queried in place rather than loaded. Checksums are
#' verified against the manifest, so a truncated or edited file is caught
#' here rather than producing quietly wrong maps.
#'
#' @param path Bundle directory; default via [sedacore_bundle_path()].
#' @param verify Check file checksums against the manifest.
#' @return A `sedacore_bundle`. Close it with [close_bundle()].
#' @keywords internal
#' @export
sedacore_open_bundle <- function(path = NULL, verify = TRUE) {
  path <- sedacore_bundle_path(path)
  manifest <- jsonlite::read_json(file.path(path, "manifest.json"), simplifyVector = FALSE)

  if (!identical(manifest$bundle_format, 1L)) {
    rlang::abort(sprintf(
      "Unsupported bundle format %s (this version reads format 1).",
      format(manifest$bundle_format %||% "unknown")
    ))
  }
  for (nm in names(manifest$tables)) {
    f <- file.path(path, manifest$tables[[nm]]$file)
    if (!file.exists(f)) rlang::abort(sprintf("Bundle table file is missing: %s", f))
    if (verify && !identical(unname(tools::md5sum(f)), manifest$tables[[nm]]$md5)) {
      rlang::abort(sprintf("Checksum mismatch for %s: the file differs from what the manifest recorded.", manifest$tables[[nm]]$file))
    }
  }

  con <- DBI::dbConnect(duckdb::duckdb())
  # If any table fails to open (e.g. a corrupt parquet file), do not leak the
  # connection on the way out.
  opened <- FALSE
  on.exit(if (!opened) disconnect_if_open(con), add = TRUE)
  for (nm in names(manifest$tables)) {
    DBI::dbExecute(con, sprintf(
      "CREATE VIEW %s AS SELECT * FROM read_parquet(%s)",
      DBI::dbQuoteIdentifier(con, nm),
      DBI::dbQuoteString(con, file.path(path, manifest$tables[[nm]]$file))
    ))
  }

  opened <- TRUE
  b <- new.env(parent = emptyenv())
  b$con <- con
  b$path <- path
  b$manifest <- manifest
  reg.finalizer(b, function(e) disconnect_if_open(e$con), onexit = TRUE)
  structure(b, class = "sedacore_bundle")
}

# Closing is idempotent: a bundle is closed explicitly *and* by the finalizer
# at garbage collection, and the second disconnect must be a quiet no-op.
#' @noRd
disconnect_if_open <- function(con) {
  if (DBI::dbIsValid(con)) DBI::dbDisconnect(con, shutdown = TRUE)
  invisible(NULL)
}

#' @rdname sedacore_open_bundle
#' @param b A `sedacore_bundle`.
#' @export
close_bundle <- function(b) {
  disconnect_if_open(b$con)
}

#' @export
print.sedacore_bundle <- function(x, ...) {
  m <- x$manifest
  cat("<sedacore_bundle> ", x$path, "\n", sep = "")
  cat("  built ", m$created, " by ", m$package, " ", m$package_version, "; region: ", m$region, "\n", sep = "")
  for (nm in names(m$tables)) cat(sprintf("  %-15s %9s rows\n", nm, format(m$tables[[nm]]$rows, big.mark = ",")))
  cat("  sources: ", paste(vapply(m$sources, function(s) paste0(s$id, " ", s$version %||% ""), ""), collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' Run a SQL query against a bundle
#'
#' Always pass values through `params` (positional `?` placeholders) rather
#' than pasting them into `sql`.
#'
#' @param b A `sedacore_bundle`.
#' @param sql A SQL string over the bundle's tables.
#' @param params A list of values for the `?` placeholders.
#' @return A tibble.
#' @keywords internal
#' @export
bundle_query <- function(b, sql, params = list()) {
  tibble::as_tibble(DBI::dbGetQuery(b$con, sql, params = unname(params)))
}

#' Read a whole bundle table
#'
#' Only sensible for the small tables (`sites`, `series_catalog`,
#' `age_control`, `site_links`); use [bundle_series()] for `series`.
#'
#' @param b A `sedacore_bundle`.
#' @param name Table name.
#' @return A tibble.
#' @keywords internal
#' @export
bundle_table <- function(b, name) {
  if (!name %in% names(b$manifest$tables)) {
    rlang::abort(sprintf("Table \"%s\" is not in this bundle (has: %s).", name, paste(names(b$manifest$tables), collapse = ", ")))
  }
  bundle_query(b, sprintf("SELECT * FROM %s", DBI::dbQuoteIdentifier(b$con, name)))
}

#' Read observations from the `series` table
#'
#' @param b A `sedacore_bundle`.
#' @param site_id,parameter_name Optional character vectors to restrict to.
#' @param age_range Optional length-2 numeric, ka BP (inclusive).
#' @param role Which roles to return; default only the central estimate
#'   (uncertainty companions are folded into `value_lo`/`value_hi`).
#' @return A tibble, ordered by site, series and age.
#' @keywords internal
#' @export
bundle_series <- function(b, site_id = NULL, parameter_name = NULL, age_range = NULL, role = "central") {
  where <- character()
  params <- list()
  add_in <- function(col, vals) {
    where <<- c(where, sprintf("%s IN (%s)", col, paste(rep("?", length(vals)), collapse = ", ")))
    params <<- c(params, as.list(vals))
  }
  if (!is.null(site_id)) add_in("site_id", site_id)
  if (!is.null(parameter_name)) add_in("parameter_name", parameter_name)
  if (!is.null(role)) add_in("role", role)
  if (!is.null(age_range)) {
    where <- c(where, "age_ka BETWEEN ? AND ?")
    params <- c(params, as.list(age_range))
  }
  sql <- paste(
    "SELECT * FROM series",
    if (length(where) > 0) paste("WHERE", paste(where, collapse = " AND ")),
    "ORDER BY site_id, series_id, age_ka"
  )
  bundle_query(b, sql, params)
}
