cat_row <- function(site, series, source, param = "d18O_planktonic", species = "G. bulloides",
                    role = "central", recalc = FALSE, n = 50L) {
  tibble::tibble(
    site_id = site, series_id = series, source = source, parameter_name = param,
    species = species, role = role, is_recalc = recalc, n_obs_aged = n
  )
}

test_that("normalize_species reconciles abbreviations, synonyms and qualifiers", {
  expect_equal(
    normalize_species(c("Cibicidoides wuellerstorfi", "C. wuellerstorfi", "Planulina wuellerstorfi", "Cibicides wuellerstorfi")),
    rep("c. wuellerstorfi", 4)
  )
  expect_equal(normalize_species("Globigerina bulloides"), "g. bulloides")
  expect_equal(normalize_species("G. ruber white"), "g. ruber white")
  expect_equal(normalize_species("Globigerinoides ruber (white)"), "g. ruber")
  expect_equal(normalize_species(c(NA, "", "  ")), rep(NA_character_, 3))
})

test_that("the same measurement from two sources yields one preferred series, by priority", {
  cat <- dplyr::bind_rows(
    cat_row("S1", "p1", "palmod", species = "G. bulloides"),
    cat_row("S1", "w1", "wpa", species = "Globigerina bulloides")
  )
  out <- mark_preferred_series(cat, c("palmod", "wpa"))
  expect_equal(out$preferred, c(TRUE, FALSE))
  expect_equal(mark_preferred_series(cat, c("wpa", "palmod"))$preferred, c(FALSE, TRUE))
})

test_that("different species, proxies and sites are not duplicates of each other", {
  cat <- dplyr::bind_rows(
    cat_row("S1", "a", "palmod", species = "G. bulloides"),
    cat_row("S1", "b", "palmod", species = "N. pachyderma"),
    cat_row("S1", "c", "palmod", param = "d13C_planktonic", species = "G. bulloides"),
    cat_row("S2", "d", "palmod", species = "G. bulloides")
  )
  expect_true(all(mark_preferred_series(cat, "palmod")$preferred))
})

test_that("planktonic and benthic isotopes are never each other's duplicates", {
  cat <- dplyr::bind_rows(
    cat_row("S1", "a", "palmod", param = "d18O_planktonic", species = NA),
    cat_row("S1", "b", "palmod", param = "d18O_benthic", species = NA)
  )
  expect_true(all(mark_preferred_series(cat, "palmod")$preferred))
})

test_that("within a source, a recalculated series beats the original, then more observations win", {
  cat <- dplyr::bind_rows(
    cat_row("S1", "orig", "palmod", param = "SST_UK37", species = NA, recalc = FALSE, n = 400L),
    cat_row("S1", "recalc", "palmod", param = "SST_UK37", species = NA, recalc = TRUE, n = 300L)
  )
  expect_equal(mark_preferred_series(cat, "palmod")$preferred, c(FALSE, TRUE))

  cat2 <- dplyr::bind_rows(
    cat_row("S1", "few", "palmod", n = 10L), cat_row("S1", "many", "palmod", n = 90L)
  )
  expect_equal(mark_preferred_series(cat2, "palmod")$preferred, c(FALSE, TRUE))
})

test_that("priority outranks recalculation and size", {
  cat <- dplyr::bind_rows(
    cat_row("S1", "p", "palmod", n = 5L, recalc = FALSE),
    cat_row("S1", "w", "wpa", n = 500L, recalc = TRUE)
  )
  expect_equal(mark_preferred_series(cat, c("palmod", "wpa"))$preferred, c(TRUE, FALSE))
})

test_that("non-central series (bounds, errors) are never preferred, and cannot block the central one", {
  cat <- dplyr::bind_rows(
    cat_row("S1", "lo", "palmod", role = "lower", n = 999L),
    cat_row("S1", "mid", "palmod", role = "central", n = 5L)
  )
  out <- mark_preferred_series(cat, "palmod")
  expect_equal(out$preferred, c(FALSE, TRUE))
})

test_that("the choice is deterministic regardless of row order", {
  cat <- dplyr::bind_rows(
    cat_row("S1", "b", "palmod", n = 50L), cat_row("S1", "a", "palmod", n = 50L)
  )
  expect_equal(mark_preferred_series(cat, "palmod")$series_id[mark_preferred_series(cat, "palmod")$preferred], "a")
  expect_equal(
    mark_preferred_series(cat[2:1, ], "palmod")$series_id[mark_preferred_series(cat[2:1, ], "palmod")$preferred], "a"
  )
})

test_that("duplicate series ids are refused (they would make the choice ambiguous)", {
  cat <- dplyr::bind_rows(cat_row("S1", "x", "palmod"), cat_row("S2", "x", "wpa"))
  expect_error(mark_preferred_series(cat, c("palmod", "wpa")), "must be unique")
})

test_that("add_preferred_to_series carries the flag onto observations, and refuses orphans", {
  cat <- mark_preferred_series(dplyr::bind_rows(
    cat_row("S1", "p1", "palmod"), cat_row("S1", "w1", "wpa")
  ), c("palmod", "wpa"))
  series <- tibble::tibble(series_id = c("p1", "p1", "w1"), value = 1:3)
  expect_equal(add_preferred_to_series(series, cat)$preferred, c(TRUE, TRUE, FALSE))
  expect_error(add_preferred_to_series(tibble::tibble(series_id = "ghost"), cat), "not in the catalog")
})
