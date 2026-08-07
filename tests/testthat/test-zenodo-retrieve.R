test_that("load_zenodo_synonym_map reads the bundled table", {
  synonyms <- load_zenodo_synonym_map()
  expect_true(all(c("raw_pattern", "match_type", "parameter_group", "parameter_name") %in% names(synonyms)))
  expect_true("Organic carbon (wt%)" %in% synonyms$raw_pattern)
})

test_that("zenodo_search_iodp errors when neither query nor expedition is given", {
  expect_error(zenodo_search_iodp())
})

test_that("zenodo_search_iodp finds the real confirmed Expedition 362 CHNS record", {
  vcr::use_cassette("zenodo_search_iodp_362", {
    res <- zenodo_search_iodp(expedition = "362", query = "CHNS", count = 10)
  })

  expect_true(all(c("doi", "title", "description", "n_files") %in% names(res)))
  expect_true(nrow(res) > 0)
  expect_true(any(res$doi == "10.5281/zenodo.19867114"))
})

test_that("zenodo_fetch_dataset retrieves and parses the real CHNS record's data", {
  vcr::use_cassette("zenodo_fetch_dataset_362_chns", {
    fetched <- zenodo_fetch_dataset("10.5281/zenodo.19867114")
  })

  expect_true(nrow(fetched) > 0)
  expect_true(all(c("doi", "title", "source_file", "data") %in% names(fetched)))
  expect_true(all(vapply(fetched$data, is.data.frame, logical(1))))
  # Real confirmed columns from this record.
  first <- fetched$data[[1]]
  expect_true("Organic carbon (wt%)" %in% names(first) || "Nitrogen (wt%)" %in% names(first))
})

test_that("zenodo_fetch_dataset output harmonizes real CHNS data via the Zenodo synonym map", {
  vcr::use_cassette("zenodo_fetch_dataset_362_chns", {
    fetched <- zenodo_fetch_dataset("10.5281/zenodo.19867114")
  })

  # Pick the first file with real (non-NA) Carbon/Nitrogen values --
  # confirmed live that this specific record's "Organic carbon (wt%)"
  # (-> TOC) is always blank, even though the column exists structurally
  # in every file (a legitimate real-world gap: this expedition's CHNS
  # export doesn't report organic-vs-total carbon separately, only total
  # elemental composition). "Carbon (wt%)" (-> total_carbon) and
  # "Nitrogen (wt%)" (-> TN) are the ones with real data here.
  has_real_carbon <- function(d) {
    all(c("Carbon (wt%)", "Nitrogen (wt%)") %in% names(d)) &&
      any(!is.na(suppressWarnings(as.numeric(d[["Carbon (wt%)"]]))))
  }
  target <- Filter(has_real_carbon, fetched$data)
  skip_if(length(target) == 0, "no fetched file had real Carbon/Nitrogen values")

  h <- harmonize_columns(
    target[[1]], core_id = "362-U1480", doi = fetched$doi[1],
    depth_col = "Top depth CSF-A (m)", synonym_map = load_zenodo_synonym_map(),
    unmatched = "silent"
  )

  expect_true(is_param_long(h))
  expect_true("total_carbon" %in% h$parameter_name)
  expect_true("TN" %in% h$parameter_name)
  expect_true(all(h$unit[h$parameter_name == "total_carbon"] == "wt%"))
  expect_true(all(h$parameter_group[h$parameter_name == "total_carbon"] == "chns"))
})

test_that("zenodo_record_id_from_doi extracts the numeric record ID", {
  expect_equal(zenodo_record_id_from_doi("10.5281/zenodo.19867114"), "19867114")
  expect_error(zenodo_record_id_from_doi("not-a-zenodo-doi"))
})

test_that("looks_like_tabular_file recognizes csv/tsv/txt only", {
  expect_true(looks_like_tabular_file("data.csv"))
  expect_true(looks_like_tabular_file("data.TSV"))
  expect_true(looks_like_tabular_file("data.txt"))
  expect_false(looks_like_tabular_file("readme.html"))
  expect_false(looks_like_tabular_file("checksums.md5"))
})
