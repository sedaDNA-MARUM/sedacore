# Module: getting data for the active core.
#
# Three routes into state$param_long: PANGAEA search + fetch, Zenodo (IODP
# community) search + fetch, and importing your own CSV. All of them end in
# the same harmonize-filter-append step (see make_fetchers()), so everything
# downstream treats every source identically.

#' @noRd
mod_retrieve_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$h4("Retrieve from PANGAEA"),
    actionButton(ns("search_btn"), "Search PANGAEA for active core"),
    DT::dataTableOutput(ns("search_results")),
    actionButton(ns("fetch_selected_btn"), "Fetch + harmonize selected", class = "btn-primary"),
    tags$hr(),
    tags$h4("Search Zenodo (IODP)"),
    helpText(
      "Recent IODP expeditions aren't necessarily mirrored to PANGAEA yet",
      "(confirmed: Expeditions 390-395 return zero PANGAEA results) --",
      "Zenodo's IODP community often has real core-level data for them",
      "already, e.g. elemental (CHNS) analyses. Only CHNS-type columns",
      "are recognized so far; other Zenodo IODP measurement types will",
      "come in as uncategorized fallback columns until added to",
      "zenodo_column_synonyms.csv. An expedition alone can match",
      "thousands of records (every instrument/data type on that",
      "expedition) -- add search terms (e.g. \"CHNS\", \"geochemistry\")",
      "to actually narrow it down to the data type you want."
    ),
    textInput(ns("zenodo_expedition"), "IODP expedition number (optional)", placeholder = "e.g. 362"),
    textInput(ns("zenodo_query"), "Search terms (optional)", placeholder = "e.g. CHNS, geochemistry"),
    actionButton(ns("zenodo_search_btn"), "Search Zenodo for active core"),
    DT::dataTableOutput(ns("zenodo_search_results")),
    actionButton(ns("zenodo_fetch_selected_btn"), "Fetch + harmonize selected", class = "btn-primary"),
    tags$hr(),
    tags$h4("Import external data"),
    helpText(
      "For data not reachable through PANGAEA or Zenodo (e.g. a paper's",
      "supplementary table, another database, a manual LORE export) --",
      "reshape it into the template below and import it here; it flows",
      "through the same harmonized/depth-matched/checklist/export",
      "pipeline as everything else, with no special-casing."
    ),
    downloadButton(ns("external_import_template_btn"), "Download template CSV"),
    fileInput(ns("external_import_file"), "External data CSV", accept = ".csv"),
    textInput(ns("external_import_source"), "Source (optional)", placeholder = "e.g. Smith et al. 2020, Table 2"),
    textInput(ns("external_import_doi"), "DOI (optional)", placeholder = "e.g. 10.1594/PANGAEA.xxxxx"),
    actionButton(ns("external_import_btn"), "Import for active core", class = "btn-primary")
  )
}

#' @noRd
mod_retrieve_server <- function(id, state, active_core, measurement_filter, fetchers,
                                search_pangaea = pg_search_cores,
                                search_zenodo = zenodo_search_iodp) {
  moduleServer(id, function(input, output, session) {
    log_msg <- function(msg) log_message(state, msg)

    # -- PANGAEA -----------------------------------------------------------
    observeEvent(input$search_btn, {
      req(state$core_registry, active_core())
      row <- state$core_registry[state$core_registry$core_id == active_core(), ]
      # Same measurement_filter as harmonize/interpolate and neighbour
      # search. When set this is real additive cost (each result has to be
      # fetched to check its actual parameters, unlike find_neighbours()'s
      # equivalent filter, which reuses fetches it makes anyway) -- hence
      # the progress bar.
      required_parameters <- if (length(measurement_filter()) > 0) measurement_filter() else NULL
      filter_note <- if (!is.null(required_parameters)) {
        sprintf(", requiring one of: %s", paste(required_parameters, collapse = ", "))
      } else {
        ""
      }
      log_msg(sprintf("Searching PANGAEA for %s%s...", active_core(), filter_note))

      res <- withProgress(message = "Searching PANGAEA", value = 0, {
        tryCatch(
          search_pangaea(
            row,
            count = 20, required_parameters = required_parameters,
            progress = function(i, n, doi) {
              setProgress(value = i / n, detail = sprintf("Checking result %d of %d...", i, n))
            }
          ),
          error = function(e) {
            log_msg(paste("ERROR:", conditionMessage(e)))
            NULL
          }
        )
      })

      if (!is.null(res)) {
        state$search_results <- res
        log_msg(sprintf("Found %d candidate datasets.", nrow(res)))
      }
    })

    output$search_results <- DT::renderDataTable({
      DT::datatable(
        state$search_results[, c("doi", "score", "size", "size_measure", "available_parameters", "citation")],
        selection = "multiple", options = list(pageLength = 5, scrollX = TRUE)
      )
    })

    observeEvent(input$fetch_selected_btn, {
      sel <- input$search_results_rows_selected
      if (length(sel) == 0) {
        log_msg("No rows selected in the search results table -- click a row to select it, then try again.")
        return(invisible())
      }
      for (doi in state$search_results$doi[sel]) {
        fetchers$pangaea(doi, active_core())
      }
    })

    # -- Zenodo ------------------------------------------------------------
    # Zenodo IODP records are per expedition (e.g. "IODP Expedition 362
    # Elemental analysis (CHNS)"), not per core, so a core_id text search
    # wouldn't match anything meaningful -- auto-fill the expedition number
    # when the active core's ID looks like an IODP hole (e.g. "362-U1480A").
    observeEvent(active_core(), {
      req(active_core())
      expedition <- iodp_expedition_from_core_id(active_core())
      if (!is.na(expedition)) {
        updateTextInput(session, "zenodo_expedition", value = expedition)
      }
    })

    observeEvent(input$zenodo_search_btn, {
      req(active_core())
      expedition <- trimws(input$zenodo_expedition %||% "")
      if (!nzchar(expedition)) expedition <- NULL
      query <- trimws(input$zenodo_query %||% "")
      if (!nzchar(query)) query <- NULL
      if (is.null(expedition) && is.null(query)) {
        log_msg(paste(
          "Enter an IODP expedition number and/or search terms to search",
          "Zenodo (expedition is auto-filled when the active core's ID",
          "looks like an IODP hole, e.g. \"362-U1480A\")."
        ))
        return(invisible())
      }
      log_msg(sprintf(
        "Searching Zenodo (expedition: %s; terms: %s)...",
        expedition %||% "any", query %||% "none"
      ))
      res <- tryCatch(
        search_zenodo(query = query, expedition = expedition, count = 20),
        error = function(e) {
          log_msg(paste("ERROR:", conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(res)) {
        state$zenodo_search_results <- res
        log_msg(sprintf("Found %d candidate Zenodo record(s).", nrow(res)))
      }
    })

    output$zenodo_search_results <- DT::renderDataTable({
      display <- state$zenodo_search_results[, c("doi", "title", "n_files", "description")]
      display$doi <- doi_link(display$doi, base_url = "https://doi.org/")
      # rownames add a display column before the data ones, so the doi
      # (1st data column) is rendered column 2; escape = -2 leaves it as raw
      # HTML (the <a> link) while every other column is escaped as normal.
      DT::datatable(
        display,
        selection = "multiple", options = list(pageLength = 5, scrollX = TRUE),
        escape = -2
      )
    })

    observeEvent(input$zenodo_fetch_selected_btn, {
      sel <- input$zenodo_search_results_rows_selected
      if (length(sel) == 0) {
        log_msg("No rows selected in the Zenodo search results table -- click a row to select it, then try again.")
        return(invisible())
      }
      for (doi in state$zenodo_search_results$doi[sel]) {
        fetchers$zenodo(doi, active_core())
      }
    })

    # -- your own data -----------------------------------------------------
    output$external_import_template_btn <- downloadHandler(
      filename = function() "external_import_template.csv",
      content = function(file) {
        file.copy(
          system.file("extdata", "external_import_template.csv", package = "sedacore"),
          file
        )
      }
    )

    observeEvent(input$external_import_btn, {
      req(active_core())
      if (is.null(input$external_import_file)) {
        log_msg("Choose a CSV file (matching the template) to import first.")
        return(invisible())
      }
      source_label <- trimws(input$external_import_source %||% "")
      if (!nzchar(source_label)) source_label <- NA_character_
      doi_label <- trimws(input$external_import_doi %||% "")
      if (!nzchar(doi_label)) doi_label <- NA_character_

      h <- tryCatch(
        import_external_data(
          input$external_import_file$datapath,
          core_id = active_core(), source = source_label, doi = doi_label
        ),
        error = function(e) {
          log_msg(paste("ERROR importing:", conditionMessage(e)))
          NULL
        }
      )
      if (is.null(h) || nrow(h) == 0) {
        return(invisible())
      }

      # Same optional measurement filter as the PANGAEA/Zenodo fetches.
      n_before_filter <- nrow(h)
      h <- keep_selected_measurements(h, measurement_filter())
      if (nrow(h) > 0) {
        state$param_long <- dplyr::bind_rows(state$param_long, h)
        log_msg(sprintf(
          "Imported %d rows for %s from %s.",
          nrow(h), active_core(), input$external_import_file$name
        ))
      } else {
        log_msg(sprintf(
          "0 rows from %s after restricting to selected measurements (%d rows available before filtering).",
          input$external_import_file$name, n_before_filter
        ))
      }
    })
  })
}
