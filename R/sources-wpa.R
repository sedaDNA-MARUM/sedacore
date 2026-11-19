# Adapter: World Atlas of late Quaternary Foraminiferal Oxygen and Carbon
# Isotope Ratios -- source table (Mulitza et al. 2022).
#
# This is the *record catalogue* of the atlas, not its data. For each of the
# atlas's ~2,100 isotope records it gives the core, position, water depth,
# species, the full reference, and where the data came from (usually a
# PANGAEA DOI). It has no values and no ages -- in the atlas itself the ages
# are applied by PaleoDataView from separate radiocarbon files, so only the
# 2024.0 release's processed time series carry them. What this table does give
# is cheap, citable and complete:
#
#   * the atlas's cores, which the crosswalk folds onto the canonical sites
#     (it contains cores PALMOD lacks, e.g. GeoB1706-2 and GeoB9064-1);
#   * `indexed_records`: records that exist and are compiled somewhere, but
#     whose age-controlled values are not in the bundle -- the "exists, not
#     (yet) usable with ages" state of the gap map;
#   * links: reference and data source per core.

#' Pull DOIs out of free text
#'
#' Elsevier DOIs contain parentheses (`10.1016/S0967-0637(99)00076-X`), so a
#' parenthesis cannot simply end a DOI; but a DOI cited *inside* parentheses
#' (`(doi:10.1029/93PA02198)`) must not swallow the closing one. Trailing
#' punctuation and any *unbalanced* trailing `)` are trimmed.
#'
#' @param text One string.
#' @return Unique DOIs (without a resolver prefix).
#' @noRd
extract_dois <- function(text) {
  d <- stringr::str_extract_all(text, "10\\.[0-9]{4,9}/[^\\s,;]+")[[1]]
  d <- sub("[.,;:]+$", "", d)
  repeat {
    unbalanced <- grepl("\\)$", d) & stringr::str_count(d, "\\)") > stringr::str_count(d, "\\(")
    if (!any(unbalanced)) break
    d[unbalanced] <- sub("[.,;:]*\\)$", "", d[unbalanced])
    d <- sub("[.,;:]+$", "", d)
  }
  unique(d)
}

#' Read the atlas source table
#'
#' @param path Path to `WA_Stable_Isotopes_2022_Sources.txt` (Zenodo
#'   10.5281/zenodo.6337519): tab-separated, Latin-1, CRLF.
#' @return A list: `tables` (`sites`, `site_links`, `indexed_records`) and
#'   `source` (the manifest entry), the same shape as [build_palmod()]. Site
#'   ids are the atlas's own core names.
#' @keywords internal
#' @export
read_wpa_sources <- function(path) {
  x <- readr::read_tsv(
    path,
    col_types = readr::cols(.default = "c"), locale = readr::locale(encoding = "latin1"),
    quote = "", show_col_types = FALSE, progress = FALSE
  )
  required <- c("Core", "Longitude", "Latitude", "Depth", "Species", "Reference Stable Isotopes", "Source Stable Isotopes")
  missing_cols <- setdiff(required, names(x))
  if (length(missing_cols) > 0) {
    rlang::abort(sprintf("Source table is missing column(s): %s", paste(missing_cols, collapse = ", ")))
  }
  x <- x[!is.na(x$Core) & nzchar(trimws(x$Core)), ]
  x$Core <- trimws(x$Core)
  x$lat <- suppressWarnings(as.numeric(x$Latitude))
  x$lon <- suppressWarnings(as.numeric(x$Longitude))
  x$depth <- suppressWarnings(as.numeric(x$Depth))
  if (anyNA(x$lat) | anyNA(x$lon)) {
    rlang::abort("Source table has records without a numeric position.")
  }
  rc_ref <- col_or(x, "Reference Radiocarbon", NA_character_)
  rc_src <- col_or(x, "Source Radiocarbon", NA_character_)
  blank_to_na <- function(v) ifelse(is.na(v) | !nzchar(trimws(v)), NA_character_, trimws(v))

  # One site per core. A core is listed once per species and the position
  # should agree between rows; if it does not, the first is used and the
  # spread is returned so it can be seen.
  by_core <- split(seq_len(nrow(x)), x$Core)
  sites <- dplyr::bind_rows(lapply(by_core, function(i) {
    tibble::tibble(
      site_id = x$Core[i[1]], label = x$Core[i[1]], latitude = x$lat[i[1]], longitude = x$lon[i[1]],
      water_depth_m = x$depth[i[1]], basin = NA_character_,
      program = classify_core_source(gsub("_", "-", x$Core[i[1]], fixed = TRUE)),
      sources = "wpa_2022", n_records = length(i),
      position_spread_km = max(haversine_km(x$lat[i[1]], x$lon[i[1]], x$lat[i], x$lon[i]))
    )
  }))

  link_rows <- function(core, text, type, label) {
    text <- blank_to_na(text)
    if (is.na(text)) {
      return(NULL)
    }
    doi <- extract_dois(text)
    if (length(doi) > 0) {
      is_data <- grepl("^10\\.1594/PANGAEA", doi)
      return(tibble::tibble(
        site_id = core, link_type = ifelse(is_data, "data", type), label = ifelse(is_data, "Data", label),
        url = paste0("https://doi.org/", doi), source = "wpa_2022"
      ))
    }
    u <- stringr::str_extract(text, "(https?|ftp)://[^\\s,;]+")
    if (is.na(u)) {
      return(NULL)
    }
    tibble::tibble(site_id = core, link_type = "data", label = "Data", url = u, source = "wpa_2022")
  }
  links <- dplyr::bind_rows(lapply(seq_len(nrow(x)), function(i) {
    dplyr::bind_rows(
      link_rows(x$Core[i], x[["Reference Stable Isotopes"]][i], "primary_publication", "Publication"),
      link_rows(x$Core[i], x[["Source Stable Isotopes"]][i], "data", "Data"),
      link_rows(x$Core[i], rc_ref[i], "chronology_publication", "Radiocarbon"),
      link_rows(x$Core[i], rc_src[i], "data", "Radiocarbon data")
    )
  }))
  links <- dplyr::distinct(links, .data$site_id, .data$link_type, .data$url, .keep_all = TRUE)

  first_url <- function(text) {
    vapply(blank_to_na(text), function(t) {
      if (is.na(t)) {
        return(NA_character_)
      }
      d <- extract_dois(t)
      if (length(d) > 0) paste0("https://doi.org/", d[1]) else stringr::str_extract(t, "(https?|ftp)://[^\\s,;]+")
    }, character(1), USE.NAMES = FALSE)
  }
  records <- tibble::tibble(
    record_id = sprintf("wpa_2022:%04d", seq_len(nrow(x))),
    source = "wpa_2022", site_id = x$Core, label = x$Core,
    latitude = x$lat, longitude = x$lon,
    proxy_class = "foram isotopes", species = blank_to_na(x$Species),
    doi = first_url(x[["Source Stable Isotopes"]]),
    citation = blank_to_na(x[["Reference Stable Isotopes"]]),
    data_source = blank_to_na(x[["Source Stable Isotopes"]]),
    has_radiocarbon = !is.na(blank_to_na(rc_ref)) | !is.na(blank_to_na(rc_src))
  )

  list(
    tables = list(sites = sites, site_links = links, indexed_records = records),
    source = list(
      id = "wpa_2022",
      name = "World Atlas of late Quaternary Foraminiferal Oxygen and Carbon Isotope Ratios (source table)",
      version = "2022.1",
      doi = "10.5281/zenodo.6337519",
      license = "CC-BY-4.0",
      citation = paste(
        "Mulitza, S. et al. (2022): Data Sources for the World Atlas of late Quaternary Foraminiferal Oxygen",
        "and Carbon Isotope Ratios 2021. Zenodo, https://doi.org/10.5281/zenodo.6337519; atlas: Mulitza et al.",
        "(2022), Earth System Science Data 14, 2553-2611, https://doi.org/10.5194/essd-14-2553-2022"
      ),
      n_sites = nrow(sites), n_records = nrow(records)
    )
  )
}
