raw <- tibble::tibble(
  `Depth sed [m]` = c(0.1, 0.2, 0.3),
  `G. inflata δ18O [‰ PDB]` = c(1.1, 1.2, 1.3),
  `G. inflata δ13C [‰ PDB]` = c(0.5, 0.6, 0.7),
  `TOC [%]` = c(2.1, 2.2, 2.3),
  `Some Unrelated Column [xyz]` = c(9, 9, 9)
)

test_that("harmonize_columns matches species-prefixed isotope columns", {
  h <- harmonize_columns(raw, core_id = "core1", doi = "doi1", unmatched = "silent")

  expect_true(is_param_long(h))
  expect_setequal(unique(h$parameter_name), c("d18O", "d13C", "TOC", "Some Unrelated Column"))
  expect_equal(unique(h$depth_cm), c(10, 20, 30))
  expect_equal(unique(h$species[h$parameter_group == "isotopes"]), "G. inflata")
  expect_true(all(is.na(h$species[h$parameter_group == "toc"])))
})

test_that("harmonize_columns rounds depth_cm to clear floating-point noise from unit conversion", {
  # 0.28 * 100 == 28.000000000000004 in IEEE 754 double arithmetic --
  # exposed by readr::write_tsv()'s full-precision serialization and by
  # exact-depth matching in interpolate_depth_series() against clean
  # sample depths, even though R's default console print hides it.
  noisy <- tibble::tibble(
    `Depth sed [m]` = c(0.28, 0.58, 1.12),
    `TOC [%]` = c(2.1, 2.2, 2.3)
  )
  h <- harmonize_columns(noisy, core_id = "core1", doi = "doi1", unmatched = "silent")
  expect_identical(sort(unique(h$depth_cm)), c(28, 58, 112))
})

test_that("harmonize_columns never drops a genuinely numeric column, even if unmatched", {
  h <- harmonize_columns(raw, core_id = "core1", doi = "doi1", unmatched = "silent")
  fallback <- h[h$parameter_name == "Some Unrelated Column", ]

  expect_equal(nrow(fallback), 3)
  expect_equal(fallback$parameter_group, rep("uncategorized", 3))
  expect_equal(fallback$unit, rep("xyz", 3))
  expect_equal(fallback$value, c(9, 9, 9))
  expect_equal(attr(h, "fallback_columns"), "Some Unrelated Column [xyz]")
  expect_length(attr(h, "dropped_columns"), 0)
})

test_that("unmatched = 'warn' warns about fallback columns but still includes them", {
  expect_warning(
    h <- harmonize_columns(raw, core_id = "core1", doi = "doi1", unmatched = "warn"),
    class = "sedacore_unmatched_columns"
  )
  expect_true(any(h$parameter_name == "Some Unrelated Column"))
})

test_that("unmatched = 'error' aborts before including fallback columns", {
  expect_error(
    harmonize_columns(raw, core_id = "core1", doi = "doi1", unmatched = "error"),
    class = "sedacore_unmatched_columns"
  )
})

test_that("harmonize_columns errors when depth_col is absent", {
  bad <- dplyr::select(raw, -`Depth sed [m]`)
  expect_error(harmonize_columns(bad, core_id = "core1", unmatched = "silent"))
})

test_that("harmonize_columns excludes a genuinely non-numeric column and reports it as dropped", {
  raw2 <- tibble::tibble(
    `Depth sed [m]` = c(0.1, 0.2),
    `Comment` = c("tie point", NA)
  )
  h <- harmonize_columns(raw2, core_id = "core1", unmatched = "silent")
  expect_equal(nrow(h), 0)
  expect_true(is_param_long(h))
  expect_equal(attr(h, "dropped_columns"), "Comment")
})

test_that("harmonize_columns matches real bulk-property column names (PANGAEA.96475)", {
  bulk <- tibble::tibble(
    `Depth sed [m]` = c(0.1, 0.2),
    `Water wm [%]` = c(45.2, 44.8),
    `Poros [% vol]` = c(78.1, 77.9),
    `DBD [g/cm**3]` = c(0.6, 0.61),
    `WBD [g/cm**3]` = c(1.4, 1.41)
  )
  h <- harmonize_columns(bulk, core_id = "core1", unmatched = "silent")
  expect_setequal(
    unique(h$parameter_name),
    c("water_content", "porosity", "bulk_density_dry", "bulk_density_wet")
  )
  expect_length(attr(h, "fallback_columns"), 0)
})

test_that("harmonize_columns matches real TOC companion columns (PANGAEA.54815)", {
  toc <- tibble::tibble(
    `Depth sed [m]` = c(0.1, 0.2),
    `CaCO3 [%]` = c(56.2, 58.5),
    `TOC [%]` = c(4.11, 3.94),
    `TN [%]` = c(0.538, 0.5),
    `C/N` = c(7.64, 7.87),
    `δ13C Corg [‰ PDB]` = c(-20.3, -20.2)
  )
  h <- harmonize_columns(toc, core_id = "core1", unmatched = "silent")
  expect_setequal(
    unique(h$parameter_name),
    c("CaCO3", "TOC", "TN", "C_N_ratio", "d13C_org")
  )
  expect_length(attr(h, "fallback_columns"), 0)
})

test_that("harmonize_columns matches real alkenone/SST biomarker columns (PANGAEA.54904)", {
  alk <- tibble::tibble(
    `Depth sed [m]` = c(0.03, 0.08),
    `C37:3Me [ng/g]` = c(1233, 935),
    `C37:2Me [ng/g]` = c(3169, 2294),
    `K37 [ng/g]` = c(4402, 3229),
    `UK'37` = c(0.72, 0.71),
    `SST (1-12) [°C]` = c(24.1, 23.9)
  )
  h <- harmonize_columns(alk, core_id = "core1", unmatched = "silent")
  expect_setequal(
    unique(h$parameter_name),
    c("C37_3Me", "C37_2Me", "K37", "UK37", "SST")
  )
  expect_true(all(h$parameter_group[h$parameter_name %in% c("C37_3Me", "C37_2Me", "K37", "UK37", "SST")] == "biomarkers"))
  expect_length(attr(h, "fallback_columns"), 0)
})

test_that("harmonize_columns matches real pore-water chemistry columns (PANGAEA.104672/87950/52909)", {
  porewater <- tibble::tibble(
    `Depth sed [m]` = c(0.1, 0.2),
    `Fe [µmol/l]` = c(1.2, 3.4),
    `[PO4]3- [µmol/l]` = c(0.5, 0.9),
    `Si(OH)4 [µmol/l]` = c(12.1, 15.3)
  )
  h <- harmonize_columns(porewater, core_id = "core1", unmatched = "silent")
  expect_setequal(unique(h$parameter_name), c("Fe", "PO4", "Si_OH_4"))
  expect_true(all(h$parameter_group == "porewater_chemistry"))
  expect_length(attr(h, "fallback_columns"), 0)
})

test_that("harmonize_columns matches age-model columns; numeric-looking fallback columns are kept, genuine text is dropped (PANGAEA.734358)", {
  age <- tibble::tibble(
    `Depth sed [m]` = c(0.1, 0.2),
    `Isotopic event` = c("2.2", "3.1"),
    `Age model [ka]` = c(12.5, 25.1),
    `Comm` = c("tie point", NA)
  )
  h <- suppressWarnings(harmonize_columns(age, core_id = "core1", unmatched = "silent"))
  expect_setequal(unique(h$parameter_name), c("age_model", "Isotopic event"))
  expect_equal(attr(h, "fallback_columns"), "Isotopic event")
  expect_equal(attr(h, "dropped_columns"), "Comm")
})

test_that("harmonize_columns matches a Bacon-style 'Age [ka BP] (Median Age)' header, plus d18O quantile uncertainty bounds", {
  # Reported by a user: an age-depth model dataset reporting the median
  # age point-estimate and quantile uncertainty bounds as separate
  # columns -- the plain 'Age [ka BP]' synonym doesn't match this exact
  # header, so age_model silently disappeared from the harmonized table.
  # Same story for the species-prefixed d18O quantile columns: the plain
  # d18O regex is anchored at the end of the string, so it doesn't match
  # these differently-suffixed uncertainty-bound headers either.
  age <- tibble::tibble(
    `Depth sed [m]` = c(0.1, 0.2),
    `Age [ka BP] (Median Age)` = c(12.5, 25.1),
    `C. wuellerstorfi δ18O [‰ PDB] (Quantile (Q0.025))` = c(2.1, 2.3),
    `C. wuellerstorfi δ18O [‰ PDB] (Quantile (Q0.975))` = c(2.5, 2.7),
    `C. wuellerstorfi δ18O [‰ PDB]` = c(2.3, 2.5)
  )
  h <- suppressWarnings(harmonize_columns(age, core_id = "core1", unmatched = "silent"))

  expect_true("age_model" %in% h$parameter_name)
  expect_equal(h$value[h$parameter_name == "age_model"], c(12.5, 25.1))
  expect_equal(h$parameter_group[h$parameter_name == "age_model"], c("age", "age"))

  expect_equal(h$value[h$parameter_name == "d18O"], c(2.3, 2.5))
  expect_equal(h$value[h$parameter_name == "d18O_q025"], c(2.1, 2.3))
  expect_equal(h$value[h$parameter_name == "d18O_q975"], c(2.5, 2.7))
  expect_true(all(h$species[h$parameter_name %in% c("d18O", "d18O_q025", "d18O_q975")] == "C. wuellerstorfi"))
  expect_true(all(h$parameter_group[h$parameter_name %in% c("d18O_q025", "d18O_q975")] == "isotopes"))
})

test_that("apply_dataset_override is a no-op when no override is registered", {
  out <- apply_dataset_override("10.1594/PANGAEA.doesnotexist", raw)
  expect_identical(out, raw)
})

test_that("apply_dataset_override applies a registered override", {
  overrides <- list(
    "10.1594/PANGAEA.test" = function(raw_data) {
      dplyr::mutate(raw_data, extra = 1)
    }
  )
  out <- apply_dataset_override("10.1594/PANGAEA.test", raw, overrides = overrides)
  expect_true("extra" %in% names(out))
})

test_that("harmonize_columns recognizes real sea-ice proxy headers (IRD, IP25, PIP25 indices, sea-ice diatoms)", {
  # Headers copied from real PANGAEA datasets: 929721 (IRD, method-
  # qualified), 877365 (IRD per g), 923852 (IRD flux), 905513 (IP25),
  # 887833 (PBIP25/PDIP25), 891102 (P(III)IP25, ng/g), 924721 (diatoms).
  sea_ice <- tibble::tibble(
    `Depth sed [m]` = c(1, 2),
    `IRD [#] (Counting >150 µm fraction)` = c(10, 20),
    `IRD [#/g]` = c(1.5, 2.5),
    `IRD flux [#/cm**2/ka]` = c(3, 4),
    `IP25/sed [µg/g]` = c(0.1, 0.2),
    `IP25/TOC [µg/g]` = c(0.3, 0.4),
    `PBIP25` = c(0.5, 0.6),
    `PDIP25` = c(0.7, 0.8),
    `P(III)IP25 (c=0.63)` = c(0.2, 0.3),
    `F. curta gr [%]` = c(5, 6)
  )
  h <- suppressWarnings(harmonize_columns(sea_ice, core_id = "core1", unmatched = "silent"))

  expect_setequal(
    unique(h$parameter_name),
    c("IRD_count", "IRD_per_g", "IRD_flux_count", "IP25_sed", "IP25_toc",
      "PBIP25", "PDIP25", "P_III_IP25", "F_curta_gr")
  )
  expect_true(all(h$parameter_group == "sea_ice"))
  expect_length(attr(h, "fallback_columns"), 0)
})

test_that("sea_ice_conc only matches PIP25-derived concentrations, not the satellite-based column", {
  conc <- tibble::tibble(
    `Depth sed [m]` = c(1, 2),
    `Sea ice conc [%] (PI(III)IP25-based)` = c(40, 50),
    `Sea ice conc [%] (satelite based)` = c(70, 80)
  )
  h <- suppressWarnings(harmonize_columns(conc, core_id = "core1", unmatched = "silent"))

  expect_equal(h$value[h$parameter_name == "sea_ice_conc"], c(40, 50))
  expect_equal(h$raw_column_name[h$parameter_name == "sea_ice_conc"], rep("Sea ice conc [%] (PI(III)IP25-based)", 2))
})

test_that("IP25 reported in ng/g converts to the canonical ug/g via unit_conversions.csv", {
  x <- new_param_long(tibble::tibble(
    core_id = "core1", pangaea_doi = "d", depth_cm = 1, parameter_group = "sea_ice",
    parameter_name = "IP25_sed", species = NA_character_, value = 500, unit = "ng/g",
    raw_column_name = "IP25/sed [ng/g]", is_interpolated = FALSE,
    interp_gap_cm = NA_real_, retrieved_at = Sys.time()
  ))
  out <- standardize_units(x)
  expect_equal(out$value, 0.5)
  expect_equal(out$unit, "µg/g")
})
