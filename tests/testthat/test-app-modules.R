# Offline tests of the Shiny modules via shiny::testServer(). The network
# calls (PANGAEA/Zenodo/neighbour search) are injected, so nothing here
# touches the internet.

mk_long <- function(core, name, depths, vals, group = "toc") {
  new_param_long(tibble::tibble(
    core_id = core, pangaea_doi = "d", depth_cm = depths, parameter_group = group,
    parameter_name = name, species = NA_character_, value = vals, unit = "u",
    raw_column_name = "r", is_interpolated = FALSE, interp_gap_cm = NA_real_,
    retrieved_at = Sys.time()
  ))
}

mk_samples <- function(core, ages) {
  new_sample_registry(tibble::tibble(
    core_id = core, depth_cm = seq(50, by = 10, length.out = length(ages)), latitude = 1, longitude = 1,
    water_depth_m = 1, sample_id = paste0("s", seq_along(ages)), age_ka = ages
  ))
}

mk_neighbours <- function(core, ids) {
  new_neighbour_registry(tibble::tibble(
    core_id = core, neighbour_doi = "x", neighbour_core_id = ids, latitude = 1, longitude = 1,
    distance_km = 1, source_program = "other", date_time = Sys.time(), citation = "c",
    available_parameters = "TOC", included_in_export = TRUE
  ))
}

test_that("new_app_state starts empty and log_message prepends timestamped lines", {
  shiny::isolate({
    st <- new_app_state()
    expect_equal(nrow(st$param_long), 0)
    expect_equal(nrow(st$depth_matched_long), 0)
    log_message(st, "first")
    log_message(st, "second")
    expect_match(st$log[1], "^\\[\\d\\d:\\d\\d:\\d\\d\\] second$")
    expect_match(st$log[2], "first$")
  })
})

# -- cores ------------------------------------------------------------------------------------

test_that("mod_cores_server auto-loads the bundled example and exposes the active core", {
  st <- new_app_state()
  shiny::testServer(mod_cores_server, args = list(state = st), {
    session$flushReact()
    expect_gt(nrow(st$core_registry), 0)
    expect_true(inherits(st$sample_registry, "sample_registry"))
    expect_true(any(grepl("^\\[.*\\] Loaded \\d+ cores from example_core_metadata.tsv\\.$", shiny::isolate(st$log))))
    session$setInputs(active_core = st$core_registry$core_id[1])
    expect_equal(session$returned(), st$core_registry$core_id[1])
  })
})

test_that("mod_cores_server logs, rather than crashes on, an unreadable core file", {
  st <- new_app_state()
  bad <- tempfile(fileext = ".tsv")
  writeLines("not\ta\tcore\tfile", bad)
  shiny::testServer(mod_cores_server, args = list(state = st), {
    session$setInputs(core_file = data.frame(datapath = bad, name = "bad.tsv"), load_btn = 1)
    expect_true(any(grepl("ERROR loading cores:", shiny::isolate(st$log))))
  })
})

# -- filter ---------------------------------------------------------------------------------------

test_that("mod_filter_server returns the selection", {
  st <- new_app_state()
  shiny::testServer(mod_filter_server, args = list(state = st), {
    session$setInputs(measurement_filter = c("TOC", "TN"))
    expect_equal(session$returned(), c("TOC", "TN"))
  })
})

test_that("known_parameter_names covers both synonym tables, including sea-ice proxies", {
  nm <- known_parameter_names()
  expect_true(all(c("TOC", "age_model", "total_carbon", "IRD_count", "IP25_sed") %in% nm))
  expect_false(is.unsorted(nm))
})

# -- fetchers (shared by retrieve + neighbours) -------------------------------------------------------

fake_fetch <- function(doi) {
  tibble::tibble(
    doi = doi,
    data = list(tibble::tibble(`Depth sed [m]` = c(0.1, 0.2), `TOC [%]` = c(1, 2), `TN [%]` = c(0.1, 0.2)))
  )
}

test_that("make_fetchers appends harmonized rows and logs (PANGAEA route)", {
  shiny::isolate({
    st <- new_app_state()
    f <- make_fetchers(st, function() NULL, fetch_pangaea = fake_fetch)
    suppressWarnings(f$pangaea("10.1594/PANGAEA.1", "core1"))
    expect_setequal(unique(st$param_long$parameter_name), c("TOC", "TN"))
    expect_true(all(st$param_long$core_id == "core1"))
    expect_true(any(grepl("Fetching 10.1594/PANGAEA.1 for core1", st$log)))
    expect_true(any(grepl("Harmonized 4 rows from 10.1594/PANGAEA.1", st$log)))
  })
})

test_that("make_fetchers applies the measurement filter as it is *at fetch time*", {
  shiny::isolate({
    st <- new_app_state()
    selected <- NULL
    f <- make_fetchers(st, function() selected, fetch_pangaea = fake_fetch)
    selected <- "TOC" # changed after the fetchers were built
    suppressWarnings(f$pangaea("d1", "core1"))
    expect_equal(unique(st$param_long$parameter_name), "TOC")
  })
})

test_that("make_fetchers logs a fetch error instead of throwing", {
  shiny::isolate({
    st <- new_app_state()
    f <- make_fetchers(st, function() NULL, fetch_pangaea = function(doi) stop("network down"))
    expect_no_error(f$pangaea("d1", "core1"))
    expect_true(any(grepl("ERROR fetching: network down", st$log)))
    expect_equal(nrow(st$param_long), 0)
  })
})

test_that("make_fetchers' Zenodo route uses Zenodo's depth column and says so when a record has no tables", {
  shiny::isolate({
    st <- new_app_state()
    zen <- function(doi) {
      tibble::tibble(
        doi = doi,
        data = list(tibble::tibble(`Top depth CSF-A (m)` = c(1, 2), `Carbon (wt%)` = c(0.5, 0.6)))
      )
    }
    f <- make_fetchers(st, function() NULL, fetch_zenodo = zen)
    suppressWarnings(f$zenodo("10.5281/zenodo.1", "362-U1480A"))
    expect_equal(unique(st$param_long$parameter_name), "total_carbon")

    f2 <- make_fetchers(st, function() NULL, fetch_zenodo = function(doi) tibble::tibble())
    f2$zenodo("10.5281/zenodo.2", "362-U1480A")
    expect_true(any(grepl("No tabular data file found in Zenodo record 10.5281/zenodo.2", st$log)))
  })
})

# -- retrieve ------------------------------------------------------------------------------------------------

new_retrieve_args <- function(st, filter = function() NULL, core = function() "T1", ...) {
  list(
    state = st, active_core = core, measurement_filter = filter,
    fetchers = make_fetchers(st, filter, fetch_pangaea = fake_fetch), ...
  )
}

test_that("mod_retrieve_server: searching PANGAEA stores results and passes the filter through", {
  st <- new_app_state()
  shiny::isolate({
    cr <- load_core_registries(
      system.file("extdata", "example_core_metadata.tsv", package = "sedacore")
    )$core_registry
    cr$core_id[1] <- "T1"
    st$core_registry <- cr
  })
  seen <- NULL
  fake_search <- function(row, count, required_parameters, progress) {
    seen <<- required_parameters
    tibble::tibble(
      core_id = "T1", score = 1, doi = "10.1594/PANGAEA.9", size = 1, size_measure = "datasets",
      citation = "c", available_parameters = "TOC"
    )
  }
  shiny::testServer(mod_retrieve_server,
    args = new_retrieve_args(st, filter = function() "TOC", search_pangaea = fake_search),
    {
      session$setInputs(search_btn = 1)
      expect_equal(seen, "TOC")
      expect_equal(st$search_results$doi, "10.1594/PANGAEA.9")
      expect_true(any(grepl("requiring one of: TOC", shiny::isolate(st$log))))
    }
  )
})

test_that("mod_retrieve_server: fetching with nothing selected tells the user how", {
  st <- new_app_state()
  shiny::testServer(mod_retrieve_server, args = new_retrieve_args(st), {
    session$setInputs(fetch_selected_btn = 1, zenodo_fetch_selected_btn = 1)
    log <- shiny::isolate(st$log)
    expect_true(any(grepl("No rows selected in the search results table", log)))
    expect_true(any(grepl("No rows selected in the Zenodo search results table", log)))
  })
})

test_that("mod_retrieve_server: Zenodo search needs an expedition or terms, and auto-fills the expedition from an IODP core id", {
  st <- new_app_state()
  shiny::testServer(mod_retrieve_server, args = new_retrieve_args(st), {
    session$setInputs(zenodo_expedition = "", zenodo_query = "", zenodo_search_btn = 1)
    expect_true(any(grepl("Enter an IODP expedition number", shiny::isolate(st$log))))
  })

  st2 <- new_app_state()
  core <- shiny::reactiveVal("362-U1480A")
  # (updateTextInput only messages the browser, so testServer can't read the
  # filled value back; the parsing itself is covered by
  # iodp_expedition_from_core_id()'s own test. Here: the observer runs clean.)
  shiny::testServer(mod_retrieve_server, args = new_retrieve_args(st2, core = core), {
    expect_no_error(session$flushReact())
  })
})

test_that("mod_retrieve_server: Zenodo search stores results", {
  st <- new_app_state()
  fake_zen <- function(query, expedition, count) {
    tibble::tibble(doi = "10.5281/zenodo.1", title = "t", description = "d", n_files = 2L)
  }
  shiny::testServer(mod_retrieve_server, args = new_retrieve_args(st, search_zenodo = fake_zen), {
    session$setInputs(zenodo_expedition = "362", zenodo_query = "CHNS", zenodo_search_btn = 1)
    expect_equal(st$zenodo_search_results$doi, "10.5281/zenodo.1")
  })
})

test_that("mod_retrieve_server: importing your own CSV appends rows with provenance, and needs a file first", {
  st <- new_app_state()
  csv <- system.file("extdata", "external_import_template.csv", package = "sedacore")
  shiny::testServer(mod_retrieve_server, args = new_retrieve_args(st), {
    session$setInputs(external_import_btn = 1)
    expect_true(any(grepl("Choose a CSV file", shiny::isolate(st$log))))

    session$setInputs(
      external_import_file = data.frame(datapath = csv, name = "mine.csv"),
      external_import_source = "Smith 2020", external_import_doi = "",
      external_import_btn = 2
    )
    expect_equal(nrow(st$param_long), 2)
    expect_true(all(st$param_long$source == "Smith 2020"))
    expect_true(all(is.na(st$param_long$pangaea_doi)))
    expect_true(any(grepl("Imported 2 rows for T1 from mine.csv", shiny::isolate(st$log))))
  })
})

# -- neighbours ------------------------------------------------------------------------------------------------

test_that("mod_neighbours_server: search stores results; refuses with no program selected", {
  st <- new_app_state()
  shiny::isolate(st$core_registry <- load_core_registries(
    system.file("extdata", "example_core_metadata.tsv", package = "sedacore")
  )$core_registry)
  seen <- list()
  fake_find <- function(registry, core, radius_km, count, source_programs, required_parameters, progress) {
    seen <<- list(radius = radius_km, programs = source_programs, req = required_parameters)
    mk_neighbours(core, "N1")
  }
  shiny::testServer(mod_neighbours_server,
    args = list(
      state = st, active_core = function() "T1", measurement_filter = function() "TOC",
      fetchers = make_fetchers(st, function() NULL, fetch_pangaea = fake_fetch), find = fake_find
    ),
    {
      session$setInputs(
        radius_km = 80, neighbour_count = 10, neighbour_source_programs = c("GeoB", "IODP"),
        neighbour_btn = 1
      )
      expect_equal(seen, list(radius = 80, programs = c("GeoB", "IODP"), req = "TOC"))
      expect_equal(st$neighbour_results$neighbour_core_id, "N1")

      session$setInputs(neighbour_source_programs = character(), neighbour_btn = 2)
      expect_true(any(grepl("Select at least one source program", shiny::isolate(st$log))))
    }
  )
})

test_that("mod_neighbours_server: fetching nothing selected is explained, and a missing core id falls back to the DOI", {
  st <- new_app_state()
  fetched <- list()
  fetchers <- list(pangaea = function(doi, core_id) fetched[[length(fetched) + 1]] <<- c(doi, core_id))
  nb <- mk_neighbours("T1", NA_character_)
  nb$neighbour_doi <- "10.1594/PANGAEA.5"
  shiny::isolate(st$neighbour_results <- nb)
  shiny::testServer(mod_neighbours_server,
    args = list(state = st, active_core = function() "T1", measurement_filter = function() NULL, fetchers = fetchers),
    {
      session$setInputs(fetch_neighbours_btn = 1)
      expect_true(any(grepl("No rows selected in the neighbours table", shiny::isolate(st$log))))

      session$setInputs(neighbour_results_rows_selected = 1, fetch_neighbours_btn = 2)
      expect_equal(fetched[[1]], c("10.1594/PANGAEA.5", "10.1594/PANGAEA.5"))
    }
  )
})

# -- integrate --------------------------------------------------------------------------------------------------------

test_that("mod_integrate_server builds the depth/age-matched table and logs why a neighbour was skipped", {
  st <- new_app_state()
  shiny::isolate({
    st$sample_registry <- mk_samples("T1", c(500, 600))
    st$param_long <- dplyr::bind_rows(
      mk_long("A", "TOC", c(0, 20), c(1, 2)), mk_long("A", "age_model", c(0, 20), c(0, 20), "age"),
      mk_long("B", "TOC", c(0, 20), c(1, 2))
    )
    st$neighbour_results <- mk_neighbours("T1", c("A", "B"))
  })
  shiny::testServer(mod_integrate_server,
    args = list(state = st, active_core = function() "T1", measurement_filter = function() NULL),
    {
      # The scenario that once crashed the whole app (NA match_basis rows
      # alongside a refused neighbour).
      expect_no_error(session$setInputs(interp_btn = 1))
      expect_gt(nrow(st$depth_matched_long), 0)
      expect_true(all(st$depth_matched_long$target_core_id == "T1"))
      expect_true(any(grepl("Skipped B: no age model available", shiny::isolate(st$log))))
    }
  )
})

test_that("mod_integrate_server re-running replaces only that core's previous rows", {
  st <- new_app_state()
  shiny::isolate({
    st$sample_registry <- dplyr::bind_rows(mk_samples("T1", c(5, 6)), mk_samples("T2", c(5, 6)))
    st$param_long <- dplyr::bind_rows(
      mk_long("T1", "TOC", c(50, 60), c(1, 2)), mk_long("T2", "TOC", c(50, 60), c(3, 4))
    )
  })
  core <- shiny::reactiveVal("T1")
  shiny::testServer(mod_integrate_server,
    args = list(state = st, active_core = core, measurement_filter = function() NULL),
    {
      session$setInputs(interp_btn = 1)
      n_t1 <- sum(st$depth_matched_long$target_core_id == "T1")
      core("T2")
      session$setInputs(interp_btn = 2)
      core("T1")
      session$setInputs(interp_btn = 3) # re-run T1
      expect_equal(sum(st$depth_matched_long$target_core_id == "T1"), n_t1) # not doubled
      expect_gt(sum(st$depth_matched_long$target_core_id == "T2"), 0) # T2 untouched
    }
  )
})

test_that("mod_integrate_server explains missing sample depths and missing data", {
  st <- new_app_state()
  shiny::isolate(st$sample_registry <- mk_samples("T1", c(1, 2)))
  shiny::testServer(mod_integrate_server,
    args = list(state = st, active_core = function() "T1", measurement_filter = function() NULL),
    {
      session$setInputs(interp_btn = 1)
      expect_true(any(grepl("No harmonized data yet for the active core; fetch some first", shiny::isolate(st$log))))
    }
  )

  st2 <- new_app_state()
  shiny::isolate(st2$sample_registry <- mk_samples("OTHER", c(1, 2)))
  shiny::testServer(mod_integrate_server,
    args = list(state = st2, active_core = function() "T1", measurement_filter = function() NULL),
    {
      session$setInputs(interp_btn = 1)
      expect_true(any(grepl("No sample depths on file for the active core", shiny::isolate(st2$log))))
    }
  )
})

test_that("the depth/age-matched table renders before any interpolation has run", {
  st <- new_app_state()
  shiny::testServer(mod_integrate_server,
    args = list(state = st, active_core = function() "T1", measurement_filter = function() NULL),
    {
      expect_no_error(output$depth_matched_table)
    }
  )
})

# -- outputs ----------------------------------------------------------------------------------------------------------------

test_that("mod_outputs_server renders every view, on empty and populated state", {
  st <- new_app_state()
  shiny::testServer(mod_outputs_server, args = list(state = st, active_core = function() "T1"), {
    for (o in c("measurement_summary", "param_table", "target_checklist_table", "target_checklist_summary", "log")) {
      expect_no_error(output[[o]])
    }
    st$param_long <- dplyr::bind_rows(
      mk_long("T1", "TOC", c(0, 10), c(1, 2)), mk_long("T1", "age_model", c(0, 10), c(0, 5), "age")
    )
    for (o in c("measurement_summary", "param_table", "target_checklist_table", "target_checklist_summary")) {
      expect_no_error(output[[o]])
    }
    expect_match(output$target_checklist_summary, "^For T1 -- ")
  })
})

test_that("mod_outputs_server's clear button empties both tables and logs it", {
  st <- new_app_state()
  shiny::isolate({
    st$param_long <- mk_long("T1", "TOC", 1, 1)
    st$depth_matched_long <- dplyr::mutate(mk_long("T1", "TOC", 1, 1), target_core_id = "T1")
  })
  shiny::testServer(mod_outputs_server, args = list(state = st, active_core = function() "T1"), {
    session$setInputs(clear_data_btn = 1)
    expect_equal(nrow(st$param_long), 0)
    expect_equal(nrow(st$depth_matched_long), 0)
    expect_true(any(grepl("Cleared harmonized and depth-matched data", shiny::isolate(st$log))))
  })
})

test_that("mod_outputs_server's exports write the harmonized table", {
  st <- new_app_state()
  shiny::isolate(st$param_long <- mk_long("T1", "TOC", c(1, 2), c(5, 6)))
  shiny::testServer(mod_outputs_server, args = list(state = st, active_core = function() "GeoB 1/2"), {
    path <- output$export_btn
    expect_true(file.exists(path))
    expect_equal(nrow(readr::read_tsv(path, show_col_types = FALSE)), 2)
  })
})

# -- the whole app ----------------------------------------------------------------------------------------------------------

test_that("the assembled app builds, and its server starts without error", {
  app <- sedacore_app()
  expect_s3_class(app, "shiny.appobj")
  # Starting the full server runs every module's startup logic (core
  # auto-load, filter seeding) together; it must not error.
  shiny::testServer(app_server, {
    expect_no_error(session$flushReact())
  })
})
