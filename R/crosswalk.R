# Crosswalk: deciding when records from different sources are the same site.
#
# The same core turns up under different spellings across sources
# (`GeoB3938_1` in PALMOD, `GeoB3938-1` in your core list, `MD95_2024` vs
# `MD95-2024`), and IODP holes under several forms (`108_658C`, `ODP 658B`,
# `108-658`). Counting each spelling as a new site would double-count
# coverage, so every source's sites are folded onto one canonical site.
#
# The rule, deliberately conservative:
#
#   * MERGE on name evidence AND proximity. Two records merge only if their
#     names are equivalent (separator-insensitive; IODP holes by site) and
#     they lie within `max_km` of each other.
#   * NEVER merge on proximity alone. GeoB18549-1 and MD95-2024 are 0.57 km
#     apart but are different cores -- finding exactly such neighbours is what
#     the app is for. Such pairs are listed for review instead.
#   * Name collisions are not merged either: the same name more than
#     `max_km` apart stays two sites (and is flagged), because
#     separator-insensitive names can collide by accident (`GeoB1-23` and
#     `GeoB12-3` both reduce to `geob123`).

#' Great-circle distance in kilometres
#'
#' @param lat1,lon1,lat2,lon2 Degrees; vectorised (recycled).
#' @return Distances in km.
#' @keywords internal
#' @export
haversine_km <- function(lat1, lon1, lat2, lon2) {
  rad <- pi / 180
  a <- sin((lat2 - lat1) * rad / 2)^2 +
    cos(lat1 * rad) * cos(lat2 * rad) * sin((lon2 - lon1) * rad / 2)^2
  2 * 6371.0088 * asin(pmin(1, sqrt(a)))
}

#' Separator-insensitive key for a core name
#'
#' Lower-cases and drops everything that is not a letter or digit, so
#' `GeoB3938_1`, `GeoB3938-1` and `geob3938 1` compare equal. Too blunt to
#' decide a merge on its own (see the notes at the top of this file).
#'
#' @param label Character vector.
#' @return Character vector of keys (`NA` stays `NA`).
#' @keywords internal
#' @export
site_name_key <- function(label) {
  ifelse(is.na(label), NA_character_, gsub("[^a-z0-9]", "", tolower(label)))
}

#' Key identifying a DSDP/ODP/IODP *site*, whatever the hole or spelling
#'
#' Holes of one site are metres apart and belong to one site. Recognises
#' leg/site/hole (`108_658C`, `108-658`), `ODP 658B`, `DSDP Site 502`,
#' expedition + U-site (`303_U1308`, `U1308A`, `IODP Site U1308`). DSDP and
#' ODP share one site numbering, so their keys share a namespace.
#'
#' @param label Character vector.
#' @return Character vector: `"dsdp-odp-658"`, `"iodp-U1308"`, or `NA` where
#'   the label is not an ocean-drilling site.
#' @keywords internal
#' @export
iodp_site_key <- function(label) {
  x <- trimws(label)
  out <- rep(NA_character_, length(x))

  # expedition + U-site, or a bare U-site: U1308, 303_U1308, IODP Site U1308A
  u <- stringr::str_match(x, stringr::regex("^(?:IODP[\\s_-]*)?(?:Site[\\s_-]*)?(?:\\d{2,3}[\\s_-])?(U\\d{4})[A-Z]?$", ignore_case = TRUE))[, 2]
  out[!is.na(u)] <- paste0("iodp-", toupper(u[!is.na(u)]))

  # "ODP 658B" / "DSDP Site 502"
  named <- stringr::str_match(x, stringr::regex("^(?:ODP|DSDP)[\\s_-]*(?:Site[\\s_-]*)?(\\d{1,4})[A-Z]?$", ignore_case = TRUE))[, 2]
  out[is.na(out) & !is.na(named)] <- paste0("dsdp-odp-", named[is.na(out) & !is.na(named)])

  # leg-site(-hole): 108_658C, 108-658, 94 502 (legs up to 210 are DSDP/ODP)
  ls <- stringr::str_match(x, "^(\\d{1,3})[\\s_-](\\d{1,4})[A-Z]?$")
  ok <- is.na(out) & !is.na(ls[, 2]) & suppressWarnings(as.integer(ls[, 2])) <= 210
  out[ok] <- paste0("dsdp-odp-", ls[ok, 3])
  out
}

#' @noRd
drilling_program <- function(leg) ifelse(leg <= 96, "DSDP", "ODP")

#' Canonical id for a group of equivalent labels
#'
#' IODP sites: `ODP-108-658` / `DSDP-94-502` when any member carries the leg,
#' else `ODP-658`; `IODP-U1308`. Everything else: the first label (callers pass
#' them in priority order) with `_` turned into `-`, which is how the core
#' names are actually written (PALMOD replaces hyphens with underscores).
#'
#' @param labels Labels in priority order.
#' @return One string.
#' @noRd
canonical_site_id <- function(labels) {
  keys <- iodp_site_key(labels)
  key <- keys[!is.na(keys)][1]
  if (!is.na(key)) {
    if (startsWith(key, "iodp-")) {
      return(paste0("IODP-", sub("^iodp-", "", key)))
    }
    site <- sub("^dsdp-odp-", "", key)
    ls <- stringr::str_match(trimws(labels), "^(\\d{1,3})[\\s_-](\\d{1,4})[A-Z]?$")
    has_leg <- !is.na(ls[, 2]) & !is.na(keys)
    if (any(has_leg)) {
      leg <- as.integer(ls[which(has_leg)[1], 2])
      return(sprintf("%s-%d-%s", drilling_program(leg), leg, site))
    }
    return(paste0("ODP-", site))
  }
  gsub("_", "-", labels[1], fixed = TRUE)
}

#' Fold sites from several sources onto canonical sites
#'
#' @param sites A data frame of every source's sites stacked together, with
#'   the bundle `sites` columns plus `source` (which source each row came
#'   from) -- see [build_palmod()].
#' @param priority Source ids, highest priority first. Decides which record
#'   supplies a canonical site's position, depth and basin, and the order
#'   labels are offered to the canonical id.
#' @param max_km Names must also agree on position within this distance to
#'   merge.
#' @param overrides Optional hand-made decisions, from
#'   [load_crosswalk_overrides()]: pairs of records to treat as one site
#'   whatever their positions say (e.g. a sign error in one source's
#'   coordinates, which a person has checked). An override whose two sources
#'   are both in `sites` but whose record is not there raises a warning.
#' @return A list: `sites` (canonical, bundle `sites` columns plus
#'   `n_aliases`), `site_aliases` (`source`, `alias`, `site_id`,
#'   `match_method`, `distance_km`), and `review` (what a human should look
#'   at: `same_name_far_apart` and `near_different_name`).
#' @keywords internal
#' @export
build_crosswalk <- function(sites, priority, max_km = 25, overrides = NULL) {
  stopifnot(all(c("site_id", "latitude", "longitude", "source") %in% names(sites)))
  sites <- sites[order(match(sites$source, priority), sites$site_id), ]
  sites$.iodp <- iodp_site_key(sites$label)
  sites$.key <- ifelse(!is.na(sites$.iodp), sites$.iodp, paste0("name:", site_name_key(sites$label)))

  # Group by key, then split each group into spatial clusters (greedy, in
  # priority order, against each cluster's representative).
  cluster_id <- integer(nrow(sites))
  next_id <- 0L
  for (idx in split(seq_len(nrow(sites)), sites$.key)) {
    reps <- integer()
    for (i in idx) {
      d <- if (length(reps)) {
        haversine_km(sites$latitude[reps], sites$longitude[reps], sites$latitude[i], sites$longitude[i])
      } else {
        numeric()
      }
      hit <- which(!is.na(d) & d <= max_km)
      if (length(hit) > 0) {
        cluster_id[i] <- cluster_id[reps[hit[1]]]
      } else {
        next_id <- next_id + 1L
        cluster_id[i] <- next_id
        reps <- c(reps, i)
      }
    }
  }
  forced <- integer()
  if (!is.null(overrides) && nrow(overrides) > 0) {
    for (k in seq_len(nrow(overrides))) {
      o <- overrides[k, ]
      if (!all(c(o$source_a, o$source_b) %in% sites$source)) next # a source not in this build
      ia <- which(sites$source == o$source_a & sites$site_id == o$alias_a)
      ib <- which(sites$source == o$source_b & sites$site_id == o$alias_b)
      if (length(ia) != 1 || length(ib) != 1) {
        rlang::warn(sprintf(
          "Crosswalk override %s/%s + %s/%s matches no record; check the spelling.",
          o$source_a, o$alias_a, o$source_b, o$alias_b
        ))
        next
      }
      if (cluster_id[ia] != cluster_id[ib]) {
        moved <- which(cluster_id == cluster_id[ib])
        cluster_id[moved] <- cluster_id[ia]
        forced <- c(forced, moved, ia)
      }
    }
  }
  sites$.cluster <- cluster_id

  groups <- split(seq_len(nrow(sites)), sites$.cluster)
  groups <- groups[order(vapply(groups, min, integer(1)))] # stable, priority order
  canonical <- character(length(groups))
  for (g in seq_along(groups)) canonical[g] <- canonical_site_id(sites$label[groups[[g]]])

  # Two clusters can still produce the same canonical id (a name collision far
  # apart); keep them distinct rather than silently merging.
  dup <- duplicated(canonical) | duplicated(canonical, fromLast = TRUE)
  for (nm in unique(canonical[dup])) {
    which_g <- which(canonical == nm)
    canonical[which_g] <- paste0(nm, "~", seq_along(which_g))
  }

  # Flag every name key spread over several far-apart positions, whether or
  # not the canonical ids ended up identical (`GeoB1-23` / `GeoB12-3` do not,
  # yet they are one key).
  cluster_key <- vapply(groups, function(m) sites$.key[m[1]], character(1))
  collisions <- tibble::tibble()
  for (k in unique(cluster_key[duplicated(cluster_key)])) {
    which_g <- which(cluster_key == k)
    ids <- canonical[which_g]
    # each member's position and depth, so a person can tell at a glance
    # whether this is a sign error, a coarse coordinate or a different core
    where <- vapply(groups[which_g], function(m) {
      i <- m[1]
      sprintf("%s: %s (%.3f, %.3f%s)", sites$source[i], sites$site_id[i], sites$latitude[i], sites$longitude[i],
        if (is.na(sites$water_depth_m[i])) "" else sprintf(", %g m", sites$water_depth_m[i])
      )
    }, character(1))
    collisions <- dplyr::bind_rows(collisions, tibble::tibble(
      type = "same_name_far_apart", site_id = ids,
      detail = sprintf("same name at %d positions >%g km apart -- %s", length(ids), max_km, paste(where, collapse = " | "))
    ))
  }

  aliases <- list()
  rows <- list()
  for (g in seq_along(groups)) {
    m <- groups[[g]]
    rep_i <- m[1]
    first_non_na <- function(col) {
      v <- sites[[col]][m]
      v <- v[!is.na(v)]
      if (length(v)) v[1] else sites[[col]][rep_i][NA]
    }
    id <- canonical[g]
    rows[[g]] <- tibble::tibble(
      site_id = id, label = id,
      latitude = sites$latitude[rep_i], longitude = sites$longitude[rep_i],
      water_depth_m = first_non_na("water_depth_m"), basin = first_non_na("basin"),
      program = classify_core_source(gsub("_", "-", id, fixed = TRUE)),
      sources = paste(unique(sites$source[m]), collapse = ","),
      n_aliases = length(m)
    )
    aliases[[g]] <- tibble::tibble(
      source = sites$source[m], alias = sites$site_id[m], site_id = id,
      match_method = ifelse(seq_along(m) == 1, "representative",
        ifelse(m %in% forced, "override", ifelse(!is.na(sites$.iodp[m]), "drilling_site_key", "name_key"))
      ),
      distance_km = haversine_km(sites$latitude[rep_i], sites$longitude[rep_i], sites$latitude[m], sites$longitude[m])
    )
  }
  out_sites <- dplyr::bind_rows(rows)

  # Different names, same place, different sources: probably the same core
  # under another name, or a genuine neighbour -- a person decides. Pairs
  # within 2 km only; O(n^2) is fine at a few thousand sites.
  near <- tibble::tibble()
  if (nrow(out_sites) > 1) {
    n <- nrow(out_sites)
    ii <- rep(seq_len(n - 1), times = (n - 1):1)
    jj <- unlist(lapply(seq_len(n - 1), function(i) (i + 1):n))
    d <- haversine_km(out_sites$latitude[ii], out_sites$longitude[ii], out_sites$latitude[jj], out_sites$longitude[jj])
    close <- which(!is.na(d) & d <= 2)
    if (length(close) > 0) {
      near <- tibble::tibble(
        type = "near_different_name", site_id = out_sites$site_id[ii[close]], other_site_id = out_sites$site_id[jj[close]],
        detail = sprintf("%.2f km apart, different names", d[close])
      )
    }
  }

  list(
    sites = out_sites,
    site_aliases = dplyr::bind_rows(aliases),
    review = dplyr::bind_rows(collisions, near)
  )
}

#' Re-key a source's tables onto canonical sites
#'
#' @param tables A named list of bundle tables from one source (`series`,
#'   `series_catalog`, `age_control`, `site_links`); `sites` is ignored.
#' @param aliases The `site_aliases` from [build_crosswalk()].
#' @param source The source these tables came from.
#' @return `tables` with `site_id` replaced by the canonical id.
#' @keywords internal
#' @export
apply_crosswalk <- function(tables, aliases, source) {
  map <- aliases[aliases$source == source, ]
  lookup <- stats::setNames(map$site_id, map$alias)
  out <- tables[setdiff(names(tables), "sites")]
  for (nm in names(out)) {
    old_id <- out[[nm]]$site_id
    new <- unname(lookup[old_id])
    # An NA site id is legitimate (an index entry not yet matched to a site)
    # and stays NA; any *named* site must have a crosswalk entry.
    unmapped <- is.na(new) & !is.na(old_id)
    if (any(unmapped)) {
      rlang::abort(sprintf(
        "Table `%s` has site ids with no crosswalk entry for source \"%s\": %s",
        nm, source, paste(utils::head(unique(old_id[unmapped]), 5), collapse = ", ")
      ))
    }
    out[[nm]]$site_id <- new
  }
  out
}

#' Find the bundle site for each of your own core names
#'
#' Matches by IODP site key, else by separator-insensitive name, ignoring
#' position (your cores' coordinates are for the *neighbour* search, not for
#' deciding identity), and checks the aliases too so a source's own spelling
#' resolves. A name that matches more than one site is `NA` rather than a
#' guess.
#'
#' @param labels Your core names.
#' @param sites,aliases The bundle's `sites` and `site_aliases`.
#' @return A tibble: `label`, `site_id` (`NA` if none/ambiguous), `match`
#'   (`"drilling_site"`, `"name"`, `"ambiguous"` or `"none"`).
#' @keywords internal
#' @export
resolve_site_labels <- function(labels, sites, aliases = NULL) {
  cand <- tibble::tibble(site_id = sites$site_id, label = sites$label)
  if (!is.null(aliases)) cand <- dplyr::bind_rows(cand, tibble::tibble(site_id = aliases$site_id, label = aliases$alias))
  cand$name_key <- site_name_key(cand$label)
  cand$iodp <- iodp_site_key(cand$label)

  purrr::map_dfr(labels, function(l) {
    ik <- iodp_site_key(l)
    hit <- if (!is.na(ik)) cand$site_id[!is.na(cand$iodp) & cand$iodp == ik] else character()
    how <- "drilling_site"
    if (length(hit) == 0) {
      hit <- cand$site_id[cand$name_key == site_name_key(l)]
      how <- "name"
    }
    hit <- unique(hit)
    tibble::tibble(
      label = l,
      site_id = if (length(hit) == 1) hit else NA_character_,
      match = if (length(hit) == 1) how else if (length(hit) > 1) "ambiguous" else "none"
    )
  })
}

#' Load hand-made crosswalk decisions
#'
#' A CSV of records to treat as one site whatever their positions say --
#' the way to act on `same_name_far_apart` items in the review once a person
#' has checked which coordinates are right. Columns: `source_a`, `alias_a`,
#' `source_b`, `alias_b`, `note`. The shipped file is empty: deciding which
#' source's coordinates are correct is a scientific call, made per case.
#'
#' @param path Path to `site_crosswalk_overrides.csv`.
#' @return A tibble with those columns.
#' @keywords internal
#' @export
load_crosswalk_overrides <- function(
    path = system.file("extdata", "site_crosswalk_overrides.csv", package = "sedacore")) {
  out <- read_csv_chr(path)
  needed <- c("source_a", "alias_a", "source_b", "alias_b")
  if (!all(needed %in% names(out))) {
    rlang::abort(sprintf("Crosswalk overrides must have columns: %s", paste(c(needed, "note"), collapse = ", ")))
  }
  out
}
