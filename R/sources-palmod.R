# Adapter: PALMOD 130k marine palaeoclimate synthesis, V2.0.0
# (Jonkers et al. 2026, PANGAEA.984602, CC BY 4.0).
#
# PALMOD ships one R list per site: `site`, `chron` (dated levels), `meta`
# (one row per column of `data`, in the same order), `data` (a wide table on
# a merged depth axis), and Bacon age ensembles (the bulk of the download,
# not needed here). This turns it into the bundle's tables.
#
# Decisions, each driven by what the data actually looks like:
#
# * AGE. Every series gets the *recalibrated Bacon age* (median and 2.5/97.5 %
#   percentiles) of its merged-depth row -- never the original publication's
#   age column, so that everything compiled from PALMOD shares one chronology
#   method (IntCal20/Marine20). The Bacon age is missing outside the age
#   model's range (about a third of rows overall; 80 % at MD95_2024), and
#   those rows are dropped from `series`. They are not hidden: each series'
#   `n_dropped_no_age` in `series_catalog` says how many were lost.
# * NAMES. `Parameter` alone is not enough: `surface.temp` is MgCa, UK37,
#   TEX86 or assemblage-based depending on `Material`, and IRD comes in
#   counts, per gram and percent. A map file decides the canonical name from
#   (Parameter, Material, unit), most specific match first.
# * UNCERTAINTY. Recalculated temperatures come as an adjacent
#   `.lo` / central / `.hi` triple. When the neighbours match in Material and
#   Species and the bounds actually bracket the central value, they are folded
#   into `value_lo`/`value_hi`; otherwise they stay separate series with a
#   `lower`/`upper` role. `.error` columns keep an `error` role. Only
#   `central` series count towards coverage.

#' @noRd
read_csv_chr <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
}

#' Load the PALMOD-to-canonical parameter map
#'
#' @param path Path to `palmod_parameter_map.csv`.
#' @return A tibble: `palmod_parameter, material, unit_raw, parameter_name, notes`.
#' @keywords internal
#' @export
load_palmod_parameter_map <- function(
    path = system.file("extdata", "palmod_parameter_map.csv", package = "sedacore")) {
  read_csv_chr(path)
}

#' Load the proxy vocabulary
#'
#' Maps a canonical `parameter_name` to a coarser `proxy_class` and a
#' `family`, which is what the coverage and gap views group by. Names not in
#' the vocabulary are classed `"other"`, never dropped.
#'
#' @param path Path to `proxy_vocabulary.csv`.
#' @return A tibble: `parameter_name, proxy_class, family, description`.
#' @keywords internal
#' @export
load_proxy_vocabulary <- function(
    path = system.file("extdata", "proxy_vocabulary.csv", package = "sedacore")) {
  read_csv_chr(path)
}

#' @noRd
palmod_resolve_name <- function(parameter, material, unit, map) {
  blank <- function(x) is.na(x) | !nzchar(x)
  vapply(seq_along(parameter), function(i) {
    hit <- map$palmod_parameter == parameter[i] &
      (blank(map$material) | (!is.na(material[i]) & map$material == material[i])) &
      (blank(map$unit_raw) | (!is.na(unit[i]) & map$unit_raw == unit[i]))
    if (!any(hit)) {
      return(NA_character_)
    }
    # most specific first: a row that pins material and/or unit beats a catch-all
    specificity <- (!blank(map$material[hit])) + (!blank(map$unit_raw[hit]))
    map$parameter_name[hit][which.max(specificity)]
  }, character(1))
}

#' @noRd
palmod_clean_type <- function(x) {
  # The source has typos ("meaasured", "measurd"); there are only two states.
  out <- ifelse(is.na(x), NA_character_, ifelse(grepl("^me", x), "measured", x))
  as.character(out)
}

#' @noRd
col_or <- function(df, name, default) {
  if (name %in% names(df)) df[[name]] else rep(default, nrow(df))
}

# Parameters that are structure, not proxies.
#' @noRd
palmod_structural <- "^(meanAge_kaBP|medianAge_kaBP|Age_kaBP_.*percentile|age|depth.*|label|note|notes)$"

#' @noRd
palmod_links <- function(site_id, meta, chron) {
  pick <- function(x, type, label) {
    x <- unique(unlist(strsplit(as.character(x[!is.na(x)]), ";", fixed = TRUE)))
    x <- trimws(sub("^doi:\\s*", "", trimws(x), ignore.case = TRUE))
    x <- x[nzchar(x)]
    if (length(x) == 0) {
      return(NULL)
    }
    # Anything carrying a DOI -- bare, or inside a doi.pangaea.de / doi.org
    # URL -- becomes one canonical https://doi.org/<doi> link, so the same
    # dataset is not listed twice under two resolvers.
    doi <- stringr::str_extract(x, "10\\.[0-9]{4,9}/[^\\s]+")
    url <- ifelse(!is.na(doi), paste0("https://doi.org/", doi),
      ifelse(grepl("^https?://", x), x, NA_character_)
    )
    keep <- !is.na(url)
    if (!any(keep)) {
      return(NULL)
    }
    tibble::tibble(site_id = site_id, link_type = type, label = label, url = url[keep], source = "palmod_v2")
  }
  dplyr::bind_rows(
    pick(col_or(meta, "PublicationDOI", NA_character_), "primary_publication", "Publication"),
    pick(col_or(meta, "PublicationLink", NA_character_), "primary_publication", "Publication"),
    pick(col_or(meta, "DataDOI", NA_character_), "data", "Data"),
    pick(col_or(meta, "DataLink", NA_character_), "data", "Data"),
    pick(col_or(chron, "ChronDOI", NA_character_), "chronology_publication", "Chronology")
  ) |>
    dplyr::distinct(.data$link_type, .data$url, .keep_all = TRUE)
}

#' Read one PALMOD site file into bundle tables
#'
#' @param path Path to a `<site>_age.RDS` file.
#' @param map,vocab The parameter map and proxy vocabulary.
#' @return A list: `site` (one row), `series`, `series_catalog`,
#'   `age_control`, `site_links`, and `skipped` (proxy columns with no entry
#'   in the parameter map, so nothing is dropped silently).
#' @keywords internal
#' @export
read_palmod_site <- function(path, map = load_palmod_parameter_map(), vocab = load_proxy_vocabulary()) {
  x <- readRDS(path)
  meta <- as.data.frame(x$meta)
  d <- x$data
  if (nrow(meta) != ncol(d)) {
    rlang::abort(sprintf(
      "%s: `meta` has %d rows but `data` has %d columns; they are expected to align.",
      basename(path), nrow(meta), ncol(d)
    ))
  }
  s <- x$site
  site_id <- s$SiteName[1]
  param <- meta$Parameter
  nrows <- nrow(d)

  col_of <- function(p) {
    j <- which(param == p)
    if (length(j) == 0) rep(NA_real_, nrows) else as.numeric(d[[j[1]]])
  }
  age <- col_of("medianAge_kaBP")
  age_lo <- col_of("Age_kaBP_2.5percentile")
  age_hi <- col_of("Age_kaBP_97.5percentile")

  # Depth: first of the known depth columns present; composite (mcd) depth is
  # also metres. A site with no depth column simply has no depth_cm.
  depth_cm <- rep(NA_real_, nrows)
  for (p in c("depth_merged", "depth_m", "depth", "depth.composite")) {
    j <- which(param == p)
    if (length(j) > 0) {
      scale <- if (identical(meta$ParameterUnit[j[1]], "cm")) 1 else 100
      depth_cm <- as.numeric(d[[j[1]]]) * scale
      break
    }
  }

  base_param <- sub("\\.(lo|hi|error)$", "", param)
  role <- ifelse(grepl("\\.lo$", param), "lower",
    ifelse(grepl("\\.hi$", param), "upper", ifelse(grepl("\\.error$", param), "error", "central"))
  )
  is_proxy <- !grepl(palmod_structural, param)
  material <- col_or(meta, "Material", NA_character_)
  species <- col_or(meta, "Species", NA_character_)
  unit <- col_or(meta, "ParameterUnit", NA_character_)
  ptype <- palmod_clean_type(col_or(meta, "ParameterType", NA_character_))
  canonical <- palmod_resolve_name(base_param, material, unit, map)

  same <- function(a, b) identical(as.character(a), as.character(b))
  vocab_class <- function(nm) {
    cl <- vocab$proxy_class[match(nm, vocab$parameter_name)]
    ifelse(is.na(cl), "other", cl)
  }

  # Fold an adjacent .lo / central / .hi triple into the central series, but
  # only if it is genuinely that: same Material and Species, and the bounds
  # bracket the central value on (nearly) every row where all three exist.
  band <- list()
  folded <- integer()
  for (j in which(is_proxy & role == "central" & param == "surface.temp")) {
    if (j > 1 && j < length(param) && param[j - 1] == "surface.temp.lo" && param[j + 1] == "surface.temp.hi" &&
      same(material[j - 1], material[j]) && same(material[j + 1], material[j]) &&
      same(species[j - 1], species[j]) && same(species[j + 1], species[j])) {
      lo <- as.numeric(d[[j - 1]])
      mid <- as.numeric(d[[j]])
      hi <- as.numeric(d[[j + 1]])
      ok <- !is.na(lo) & !is.na(mid) & !is.na(hi)
      if (any(ok) && mean(lo[ok] <= mid[ok] & mid[ok] <= hi[ok]) >= 0.95) {
        band[[as.character(j)]] <- list(lo = lo, hi = hi)
        folded <- c(folded, j - 1L, j + 1L)
      }
    }
  }

  series <- list()
  catalog <- list()
  skipped <- list()
  for (j in which(is_proxy)) {
    if (j %in% folded) next
    name <- canonical[j]
    if (is.na(name)) {
      skipped[[length(skipped) + 1]] <- tibble::tibble(site_id = site_id, palmod_parameter = param[j], material = material[j])
      next
    }
    vals <- suppressWarnings(as.numeric(d[[j]]))
    keep <- !is.na(vals) & !is.na(age)
    sid <- meta$TSuID[j]
    if (is.na(sid)) sid <- sprintf("%s_col%d", site_id, j)
    b <- band[[as.character(j)]]
    recalc <- grepl("^recalc", names(d)[j], ignore.case = TRUE)

    catalog[[length(catalog) + 1]] <- tibble::tibble(
      site_id = site_id, series_id = sid, source = "palmod_v2",
      parameter_name = name, proxy_class = vocab_class(name),
      species = species[j], material = material[j], parameter_type = ptype[j],
      role = role[j], is_recalc = recalc, unit = unit[j],
      n_obs = as.integer(sum(!is.na(vals))), n_obs_aged = as.integer(sum(keep)),
      n_dropped_no_age = as.integer(sum(!is.na(vals) & is.na(age))),
      age_min_ka = if (any(keep)) min(age[keep]) else NA_real_,
      age_max_ka = if (any(keep)) max(age[keep]) else NA_real_,
      has_band = !is.null(b)
    )
    if (any(keep)) {
      series[[length(series) + 1]] <- tibble::tibble(
        site_id = site_id, series_id = sid, source = "palmod_v2",
        parameter_name = name, proxy_class = vocab_class(name),
        species = species[j], material = material[j], role = role[j], is_recalc = recalc,
        depth_cm = depth_cm[keep], age_ka = age[keep], age_lo_ka = age_lo[keep], age_hi_ka = age_hi[keep],
        value = vals[keep],
        value_lo = if (is.null(b)) NA_real_ else b$lo[keep],
        value_hi = if (is.null(b)) NA_real_ else b$hi[keep],
        unit = unit[j]
      )
    }
  }

  ch <- as.data.frame(x$chron)
  age_control <- if (nrow(ch) > 0) {
    tibble::tibble(
      site_id = site_id, source = "palmod_v2",
      control_type = as.character(col_or(ch, "ChronType", NA_character_)),
      depth_cm = as.numeric(col_or(ch, "ChronDepthMid_cm", NA_real_)),
      age_14c_ka = as.numeric(col_or(ch, "ChronAge14C_kaBP", NA_real_)),
      age_14c_err_up_ka = as.numeric(col_or(ch, "ChronAge14CErrorUp_ka", NA_real_)),
      age_14c_err_down_ka = as.numeric(col_or(ch, "ChronAge14CErrorDown_ka", NA_real_)),
      reservoir_age_ka = as.numeric(col_or(ch, "ChronReservoirAge_ka", NA_real_)),
      cal_age_ka = as.numeric(col_or(ch, "ChronCalendarAge_kaBP", NA_real_)),
      cal_age_min_ka = as.numeric(col_or(ch, "ChronCalendarAgeMin_kaBP", NA_real_)),
      cal_age_max_ka = as.numeric(col_or(ch, "ChronCalendarAgeMax_kaBP", NA_real_)),
      material = as.character(col_or(ch, "ChronDatedMaterial", NA_character_)),
      species = as.character(col_or(ch, "ChronDatedSpecies", NA_character_)),
      rejected = as.logical(as.vector(col_or(ch, "ChronAgeRejected", NA))),
      reference_doi = as.character(col_or(ch, "ChronDOI", NA_character_))
    )
  } else {
    NULL
  }

  list(
    site = tibble::tibble(
      site_id = site_id, label = site_id, latitude = as.numeric(s$SiteLat[1]),
      longitude = as.numeric(s$SiteLon[1]),
      # PALMOD stores depth as a negative elevation (-1225); the atlas wants
      # water depth as a positive number of metres.
      water_depth_m = abs(as.numeric(s$SiteDepth_m[1])),
      basin = as.character(s$OceanBasin[1]),
      program = classify_core_source(gsub("_", "-", site_id, fixed = TRUE)),
      sources = "palmod_v2"
    ),
    series = dplyr::bind_rows(series),
    series_catalog = dplyr::bind_rows(catalog),
    age_control = age_control,
    site_links = palmod_links(site_id, meta, ch),
    skipped = dplyr::bind_rows(skipped)
  )
}

#' Build the PALMOD part of a data bundle
#'
#' Reads every `*.RDS` site file in `dir` and combines them into the bundle
#' tables, plus a source description for the manifest.
#'
#' @param dir Directory of PALMOD `<site>_age.RDS` files.
#' @param files Optionally, only these file names (for tests).
#' @param progress Optional `function(i, n, file)`.
#' @return A list: `tables` (named list ready for [write_bundle()]), `source`
#'   (manifest entry) and `skipped` (proxy columns with no map entry).
#' @keywords internal
#' @export
build_palmod <- function(dir, files = NULL, progress = NULL) {
  if (is.null(files)) files <- list.files(dir, pattern = "RDS$")
  if (length(files) == 0) rlang::abort(sprintf("No PALMOD .RDS files found in \"%s\".", dir))
  map <- load_palmod_parameter_map()
  vocab <- load_proxy_vocabulary()

  parts <- lapply(seq_along(files), function(i) {
    if (!is.null(progress)) progress(i, length(files), files[i])
    read_palmod_site(file.path(dir, files[i]), map, vocab)
  })
  bind <- function(nm) dplyr::bind_rows(lapply(parts, `[[`, nm))

  sites <- bind("site")
  if (anyDuplicated(sites$site_id)) {
    rlang::abort(sprintf("Duplicate PALMOD site ids: %s", paste(unique(sites$site_id[duplicated(sites$site_id)]), collapse = ", ")))
  }
  version <- readRDS(file.path(dir, files[1]))$PALMODversion %||% "unknown"

  list(
    tables = list(
      sites = sites, series = bind("series"), series_catalog = bind("series_catalog"),
      age_control = bind("age_control"), site_links = bind("site_links")
    ),
    source = list(
      id = "palmod_v2",
      name = "PALMOD 130k marine palaeoclimate data synthesis",
      version = version,
      doi = "10.1594/PANGAEA.984602",
      license = "CC-BY-4.0",
      citation = paste(
        "Jonkers, L., Hollstein, M., Siccha, M., Kucera, M. (2026): PALMOD 130k marine",
        "palaeoclimate data synthesis V2_0_0. PANGAEA, https://doi.org/10.1594/PANGAEA.984602;",
        "Earth System Science Data 18(5), 3013-3068, https://doi.org/10.5194/essd-18-3013-2026"
      ),
      n_sites = nrow(sites)
    ),
    skipped = bind("skipped")
  )
}
