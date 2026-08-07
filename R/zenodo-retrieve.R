#' @noRd
zenodo_base_url <- "https://zenodo.org/api/records"

#' Load the Zenodo/IODP column-synonym lookup table
#'
#' Same shape and purpose as [load_synonym_map()], but for
#' Zenodo/IODP-convention column headers (e.g. `"Organic carbon (wt%)"`,
#' `"Top depth CSF-A (m)"`), which don't match PANGAEA's own conventions
#' (`"TOC [%]"`, `"Depth sed [m]"`) -- pass this as [harmonize_columns()]'s
#' `synonym_map` argument when harmonizing Zenodo-sourced data. Only
#' CHNS/elemental-analysis columns are confirmed and seeded so far (see
#' `inst/extdata/zenodo_column_synonyms.csv`); other Zenodo IODP
#' measurement types will need their own header mappings added
#' incrementally, the same way `column_synonyms.csv` grew over time.
#'
#' @param path Path to `zenodo_column_synonyms.csv`.
#' @return A tibble: `raw_pattern, match_type, parameter_group,
#'   parameter_name, unit_raw, species_capture_group, notes`.
#' @export
load_zenodo_synonym_map <- function(
    path = system.file("extdata", "zenodo_column_synonyms.csv", package = "sedacore")) {
  readr::read_csv(
    path,
    show_col_types = FALSE,
    col_types = readr::cols(.default = "c")
  )
}

#' @noRd
strip_html <- function(x) {
  if (is.null(x) || is.na(x)) {
    return(NA_character_)
  }
  stringr::str_squish(stringr::str_remove_all(x, "<[^>]+>"))
}

#' Search Zenodo's IODP community for datasets
#'
#' Confirmed live to matter: recent JOIDES Resolution expeditions (e.g.
#' 390-395) aren't yet mirrored to PANGAEA at all
#' (`pangaear::pg_search(query = "campaign:392")` returns zero results),
#' but real core-level geochemistry for them is archived here -- e.g.
#' Expedition 362's "Elemental analysis (CHNS)" record. Uses Zenodo's own
#' documented REST API directly (`crul`, matching how `pangaear` itself
#' makes requests) -- no auth needed for public reads.
#'
#' @param query Free-text search query (matched against title/description/
#'   etc.), or `NULL` to search by `expedition` alone.
#' @param expedition Optional IODP expedition number (e.g. `"362"`) --
#'   appended to the query as `"Expedition <n>"`, since Zenodo has no
#'   dedicated expedition-number metadata field to filter on directly.
#' @param count Max results to request.
#' @return A tibble: `doi`, `title`, `description`, `n_files`.
#' @export
zenodo_search_iodp <- function(query = NULL, expedition = NULL, count = 20) {
  if (is.null(query) && is.null(expedition)) {
    rlang::abort("Provide at least one of `query` or `expedition`.")
  }
  q <- if (!is.null(expedition)) {
    paste(c(sprintf("Expedition %s", expedition), query), collapse = " ")
  } else {
    query
  }

  # No `sort` param -- Zenodo's default is relevance ("bestmatch"), which
  # matters a lot: confirmed live that `sort = "mostrecent"` completely
  # overrides text-match quality, burying an exact title match (e.g. the
  # real Expedition 362 CHNS record) under newer, unrelated results
  # (bathymetry/navigation from 390/391) that merely share the query
  # terms' recency-community membership.
  cli <- crul::HttpClient$new(url = zenodo_base_url)
  res <- tryCatch(
    cli$get(query = list(communities = "iodp", q = q, size = count)),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    rlang::abort(sprintf("Failed to search Zenodo: %s", conditionMessage(res)))
  }
  res$raise_for_status()

  parsed <- jsonlite::fromJSON(res$parse("UTF-8"), simplifyVector = FALSE)
  hits <- parsed$hits$hits
  if (length(hits) == 0) {
    return(tibble::tibble(
      doi = character(), title = character(), description = character(), n_files = integer()
    ))
  }

  purrr::map_dfr(hits, function(h) {
    tibble::tibble(
      doi = h$doi %||% NA_character_,
      title = h$metadata$title %||% NA_character_,
      description = strip_html(h$metadata$description),
      n_files = length(h$files %||% list())
    )
  })
}

#' @noRd
zenodo_record_id_from_doi <- function(doi) {
  # Zenodo DOIs look like "10.5281/zenodo.19867114" -- the record ID is
  # the numeric suffix after the last dot.
  m <- stringr::str_match(doi, "zenodo\\.([0-9]+)$")
  if (is.na(m[1, 2])) {
    rlang::abort(sprintf("Could not parse a Zenodo record ID out of DOI \"%s\".", doi))
  }
  m[1, 2]
}

#' @noRd
looks_like_tabular_file <- function(filename) {
  grepl("\\.(csv|tsv|txt)$", filename, ignore.case = TRUE)
}

#' @noRd
read_zenodo_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  # suppressWarnings(): real IODP LIMS export columns are often sparse
  # (e.g. "Organic carbon (wt%)" can be entirely blank in a file that
  # otherwise has real Carbon/Nitrogen/Sulfur data -- confirmed live),
  # which readr flags as a "parsing issue" (it guessed logical instead of
  # double for an all-NA column) even though nothing is actually wrong;
  # harmonize_columns() handles the resulting NAs correctly either way.
  first <- if (identical(ext, "csv")) {
    tryCatch(suppressWarnings(readr::read_csv(path, show_col_types = FALSE, guess_max = 10000)), error = function(e) NULL)
  } else {
    tryCatch(suppressWarnings(readr::read_tsv(path, show_col_types = FALSE, guess_max = 10000)), error = function(e) NULL)
  }
  # A single resulting column is a strong sign the delimiter guessed from
  # the file extension was wrong (some Zenodo IODP files use ".txt" for
  # what's actually comma-delimited, or vice versa) -- try the other one
  # before giving up.
  if (is.null(first) || ncol(first) <= 1) {
    second <- if (identical(ext, "csv")) {
      tryCatch(suppressWarnings(readr::read_tsv(path, show_col_types = FALSE, guess_max = 10000)), error = function(e) NULL)
    } else {
      tryCatch(suppressWarnings(readr::read_csv(path, show_col_types = FALSE, guess_max = 10000)), error = function(e) NULL)
    }
    if (!is.null(second) && ncol(second) > 1) {
      first <- second
    }
  }
  if (is.null(first) || ncol(first) <= 1) {
    return(NULL)
  }
  tibble::as_tibble(first, .name_repair = "minimal")
}

#' Fetch a Zenodo IODP record's data
#'
#' Downloads a record's file(s) -- there can be several per record (e.g.
#' `data_by_hole.zip`, `data_by_splice.zip`, a README, a site-locations
#' HTML file, all on the same record) -- unzips any `.zip` files, and
#' reads every real tabular file found (`.csv`/`.tsv`/`.txt`, delimiter
#' auto-detected) into a raw data frame, skipping anything else (HTML,
#' READMEs, checksums) without erroring. Mirrors [pg_fetch_dataset()]'s
#' multi-child-row shape: one row per data file found, so the same
#' per-row [harmonize_columns()] loop pattern used for PANGAEA data
#' (`fetch_and_harmonize()` in the Shiny app) works here unchanged.
#'
#' @param doi A Zenodo DOI (e.g. `"10.5281/zenodo.19867114"`).
#' @param retries Number of retry attempts on failure fetching the record
#'   metadata itself (not per-file downloads, which are best-effort and
#'   silently skipped on failure -- one bad file within a many-file
#'   record shouldn't fail the whole fetch).
#' @return A tibble with one row per data file found: `doi`, `title`,
#'   `source_file` (the file it came from, for provenance), and `data`
#'   (a list-column holding each file's raw parsed table). Zero rows if
#'   no tabular file was found in the record at all.
#' @export
zenodo_fetch_dataset <- function(doi, retries = 2) {
  check_utf8_locale()
  record_id <- zenodo_record_id_from_doi(doi)

  cli <- crul::HttpClient$new(url = paste0(zenodo_base_url, "/", record_id))
  attempt <- 0
  repeat {
    attempt <- attempt + 1
    res <- tryCatch(cli$get(), error = function(e) e)
    if (!inherits(res, "error")) break
    if (attempt > retries) {
      rlang::abort(sprintf(
        "Failed to fetch Zenodo record %s after %d attempt(s): %s",
        doi, attempt, conditionMessage(res)
      ))
    }
    Sys.sleep(1)
  }
  res$raise_for_status()
  record <- jsonlite::fromJSON(res$parse("UTF-8"), simplifyVector = FALSE)

  title <- record$metadata$title %||% NA_character_
  files <- record$files %||% list()

  tmp_dir <- tempfile("zenodo_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  data_files <- list()
  for (f in files) {
    fname <- f$key
    url <- f$links$self
    if (is.null(fname) || is.null(url)) {
      next
    }
    dest <- file.path(tmp_dir, fname)
    # Fetch into memory and write the bytes ourselves, rather than
    # crul's disk = dest streaming mode -- confirmed live that
    # vcr/webmockr's request mocking doesn't reliably intercept/replay
    # disk-streamed downloads (the file silently ends up missing/empty
    # under a cassette even though the call itself doesn't error),
    # whereas an in-memory response body is exactly what a cassette is
    # built to capture.
    ok <- tryCatch(
      {
        file_res <- crul::HttpClient$new(url = url)$get()
        file_res$raise_for_status()
        writeBin(file_res$content, dest)
        TRUE
      },
      error = function(e) FALSE
    )
    if (!isTRUE(ok)) {
      next
    }

    if (grepl("\\.zip$", fname, ignore.case = TRUE)) {
      extract_dir <- file.path(tmp_dir, tools::file_path_sans_ext(fname))
      dir.create(extract_dir, showWarnings = FALSE)
      tryCatch(utils::unzip(dest, exdir = extract_dir), error = function(e) NULL)
      extracted <- list.files(extract_dir, full.names = TRUE, recursive = TRUE)
      for (tf in extracted[looks_like_tabular_file(extracted)]) {
        parsed <- read_zenodo_table(tf)
        if (!is.null(parsed)) {
          data_files[[length(data_files) + 1]] <- list(source_file = basename(tf), data = parsed)
        }
      }
    } else if (looks_like_tabular_file(fname)) {
      parsed <- read_zenodo_table(dest)
      if (!is.null(parsed)) {
        data_files[[length(data_files) + 1]] <- list(source_file = fname, data = parsed)
      }
    }
  }

  if (length(data_files) == 0) {
    return(tibble::tibble(
      doi = character(), title = character(), source_file = character(), data = list()
    ))
  }

  tibble::tibble(
    doi = doi,
    title = title,
    source_file = vapply(data_files, function(x) x$source_file, character(1)),
    data = lapply(data_files, function(x) x$data)
  )
}
