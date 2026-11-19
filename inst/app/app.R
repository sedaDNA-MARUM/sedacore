library(shiny)
library(leaflet)
library(sedacore)

ui <- fluidPage(
  tags$head(
    tags$link(rel = "preconnect", href = "https://fonts.googleapis.com"),
    tags$link(
      rel = "stylesheet",
      href = "https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,600&display=swap"
    ),
    # The sidebar (core selection, PANGAEA search, neighbour search,
    # interpolate, export) is much taller than the main panel -- without
    # this, scrolling down to reach the lower sidebar controls scrolls
    # the map/tables out of view entirely. Giving the sidebar its own
    # internal scroll (pinned to the viewport) keeps the page itself
    # only as tall as the main panel, so the map/tables never scroll
    # out of view.
    tags$style(HTML(
      ".col-sm-4 > .well {
        position: sticky;
        top: 15px;
        max-height: calc(100vh - 30px);
        overflow-y: auto;
      }"
    ))
  ),
  titlePanel(
    tags$div(
      tags$span(
        "sedacore",
        style = paste(
          "font-family: 'Fraunces', Georgia, serif; font-weight: 600;",
          "font-size: 3.6rem; color: #e6550d; letter-spacing: 0.3px;"
        )
      ),
      tags$span(
        " — core contextual data explorer",
        style = "font-weight: 300; color: #666; font-size: 1.6rem; margin-left: 6px;"
      )
    ),
    windowTitle = "sedacore"
  ),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      fileInput("core_file", "Core metadata TSV",
        accept = c(".tsv", ".txt"),
        placeholder = "(bundled example used if empty)"
      ),
      actionButton("load_btn", "Load / reload cores", icon = icon("upload")),
      tags$hr(),
      selectInput("active_core", "Active core", choices = NULL),
      tags$hr(),
      tags$h4("Measurement filter"),
      selectizeInput("measurement_filter", "Restrict to selected measurements",
        choices = NULL, multiple = TRUE,
        options = list(placeholder = "(leave empty = harmonize/interpolate everything)")
      ),
      actionLink("use_target_checklist_btn", "Use target checklist"),
      helpText(
        "PANGAEA data is wildly inconsistent across datasets -- when set,",
        "this restricts harmonizing and interpolating to only the selected",
        "parameter_name(s), instead of pulling in everything a dataset",
        "happens to contain, and also requires PANGAEA search and",
        "neighbour search results to have at least one of these",
        "parameters (e.g. only show results with TOC data) -- checking",
        "this for a plain search costs extra requests (each result has",
        "to be fetched), so it's slower than leaving it empty. Leave",
        "empty for the default (keep/show everything)."
      ),
      tags$hr(),
      tags$h4("Retrieve from PANGAEA"),
      actionButton("search_btn", "Search PANGAEA for active core"),
      DT::dataTableOutput("search_results"),
      actionButton("fetch_selected_btn", "Fetch + harmonize selected", class = "btn-primary"),
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
      textInput("zenodo_expedition", "IODP expedition number (optional)", placeholder = "e.g. 362"),
      textInput("zenodo_query", "Search terms (optional)", placeholder = "e.g. CHNS, geochemistry"),
      actionButton("zenodo_search_btn", "Search Zenodo for active core"),
      DT::dataTableOutput("zenodo_search_results"),
      actionButton("zenodo_fetch_selected_btn", "Fetch + harmonize selected", class = "btn-primary"),
      tags$hr(),
      tags$h4("Import external data"),
      helpText(
        "For data not reachable through PANGAEA or Zenodo (e.g. a paper's",
        "supplementary table, another database, a manual LORE export) --",
        "reshape it into the template below and import it here; it flows",
        "through the same harmonized/depth-matched/checklist/export",
        "pipeline as everything else, with no special-casing."
      ),
      downloadButton("external_import_template_btn", "Download template CSV"),
      fileInput("external_import_file", "External data CSV", accept = ".csv"),
      textInput("external_import_source", "Source (optional)", placeholder = "e.g. Smith et al. 2020, Table 2"),
      textInput("external_import_doi", "DOI (optional)", placeholder = "e.g. 10.1594/PANGAEA.xxxxx"),
      actionButton("external_import_btn", "Import for active core", class = "btn-primary"),
      tags$hr(),
      tags$h4("Neighbouring cores"),
      numericInput("radius_km", "Radius (km)", value = 50, min = 1, max = 2000),
      numericInput("neighbour_count", "Max candidates to check", value = 50, min = 5, max = 300, step = 10),
      helpText(
        "PANGAEA's search relevance isn't distance-based, so in densely-",
        "sampled regions (e.g. well-studied MARUM survey areas) real nearby",
        "cores can rank surprisingly low -- raise this if you suspect",
        "known-close cores are missing from the results. Higher values",
        "take longer (each candidate is a separate PANGAEA request)."
      ),
      checkboxGroupInput("neighbour_source_programs", "Source program(s)",
        choices = c("GeoB", "IODP", "ODP", "DSDP", "other"),
        selected = c("GeoB", "IODP", "ODP", "DSDP", "other"),
        inline = TRUE
      ),
      actionButton("neighbour_btn", "Find neighbours"),
      DT::dataTableOutput("neighbour_results"),
      actionButton("fetch_neighbours_btn", "Fetch + harmonize selected neighbours", class = "btn-primary"),
      tags$hr(),
      tags$h4("Interpolate"),
      actionButton("interp_btn", "Interpolate to active core's sample depths"),
      tags$hr(),
      downloadButton("export_btn", "Export harmonized TSV"),
      downloadButton("export_depth_matched_btn", "Export depth-matched TSV"),
      actionButton("clear_data_btn", "Clear harmonized + depth-matched data", class = "btn-danger")
    ),
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("Map", leaflet::leafletOutput("map", height = 600)),
        tabPanel("Measurements", DT::dataTableOutput("measurement_summary")),
        tabPanel(
          "Target checklist",
          tags$p(
            "Fixed list of measurements of interest for sedaDNA work (edit ",
            tags$code("inst/extdata/target_measurements.csv"),
            " to change it) -- checked against the active core's own retrieved data first",
            " (never interpolated estimates), falling back to already-fetched neighbour",
            " cores' data (flagged \"found in neighbour\", not conflated with the core's",
            " own real data) so it's clear what still needs sourcing from PANGAEA or",
            " measuring at MARUM.",
            style = "color: #666;"
          ),
          textOutput("target_checklist_summary"),
          tags$br(),
          DT::dataTableOutput("target_checklist_table")
        ),
        tabPanel("Harmonized data", DT::dataTableOutput("param_table")),
        tabPanel("Depth/age-matched data", DT::dataTableOutput("depth_matched_table")),
        tabPanel("Log", verbatimTextOutput("log"))
      )
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    core_registry = NULL,
    sample_registry = NULL,
    param_long = empty_param_long(),
    depth_matched_long = empty_depth_matched_long(),
    search_results = empty_search_results(),
    zenodo_search_results = empty_zenodo_search_results(),
    neighbour_results = empty_neighbour_registry(),
    log = character()
  )

  append_log <- function(msg) {
    rv$log <- c(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg), rv$log)
  }

  # Known parameter_name vocabulary from *every* synonym table (PANGAEA +
  # Zenodo/IODP), not just PANGAEA's -- so Zenodo-only names (total_carbon
  # etc.) are selectable in the measurement filter from the start too.
  known_parameter_names <- function() {
    sort(unique(c(load_synonym_map()$parameter_name, load_zenodo_synonym_map()$parameter_name)))
  }

  # Seed the measurement-filter choices up front, then grow it with
  # anything new that shows up via harmonize's fallback path (a raw
  # column not yet in either table) so those stay selectable too,
  # without ever dropping the user's current selection. "age_model"
  # stays a selectable choice here too -- selecting it still means
  # something (requiring it when searching/filtering neighbours by
  # available_parameters), even though it's never *excluded* from
  # storage or the Harmonized data table by this same control (see
  # keep_selected_measurements() and the age_ka column below).
  updateSelectizeInput(session, "measurement_filter", choices = known_parameter_names())
  observeEvent(rv$param_long, {
    known <- isolate(input$measurement_filter)
    all_choices <- sort(unique(c(known_parameter_names(), rv$param_long$parameter_name)))
    updateSelectizeInput(session, "measurement_filter",
      choices = all_choices, selected = known
    )
  })

  observeEvent(input$use_target_checklist_btn, {
    updateSelectizeInput(session, "measurement_filter",
      choices = sort(unique(c(known_parameter_names(), rv$param_long$parameter_name, target_measurement_names()))),
      selected = target_measurement_names()
    )
  })

  load_cores <- function(path) {
    registries <- load_core_registries(path)
    sr <- registries$sample_registry
    cr <- registries$core_registry
    rv$sample_registry <- sr
    rv$core_registry <- cr
    updateSelectInput(session, "active_core", choices = cr$core_id, selected = cr$core_id[1])
    append_log(sprintf("Loaded %d cores from %s.", nrow(cr), basename(path)))
  }

  observeEvent(input$load_btn,
    {
      path <- if (!is.null(input$core_file)) {
        input$core_file$datapath
      } else {
        system.file("extdata", "example_core_metadata.tsv", package = "sedacore")
      }
      tryCatch(load_cores(path), error = function(e) append_log(paste("ERROR loading cores:", conditionMessage(e))))
    },
    ignoreNULL = FALSE
  )

  output$map <- leaflet::renderLeaflet({
    req(rv$core_registry)
    leaflet_core_map(rv$core_registry, neighbours = rv$neighbour_results)
  })

  observeEvent(input$map_marker_click, {
    click <- input$map_marker_click
    if (!is.null(click$id) && !startsWith(click$id, "neighbour-")) {
      updateSelectInput(session, "active_core", selected = click$id)
    }
  })

  # Shared by fetch_and_harmonize() (PANGAEA) and fetch_and_harmonize_zenodo()
  # (Zenodo). The harmonize + measurement-filter + message logic lives in
  # harmonize_for_store() (tested); this just appends the rows and logs.
  harmonize_and_store <- function(raw_data, core_id, doi, depth_col, synonym_map, synonym_csv_name) {
    res <- harmonize_for_store(
      raw_data, core_id, doi, depth_col, synonym_map,
      measurement_filter = input$measurement_filter, synonym_csv_name = synonym_csv_name
    )
    if (!is.null(res$data)) {
      rv$param_long <- dplyr::bind_rows(rv$param_long, res$data)
    }
    for (msg in res$messages) append_log(msg)
  }

  fetch_and_harmonize <- function(doi, core_id) {
    append_log(sprintf("Fetching %s for %s...", doi, core_id))
    fetched <- tryCatch(pg_fetch_dataset(doi), error = function(e) {
      append_log(paste("ERROR fetching:", conditionMessage(e)))
      NULL
    })
    if (is.null(fetched)) {
      return(invisible())
    }
    for (i in seq_len(nrow(fetched))) {
      harmonize_and_store(
        fetched$data[[i]], core_id, fetched$doi[i],
        depth_col = "Depth sed [m]", synonym_map = load_synonym_map(),
        synonym_csv_name = "column_synonyms.csv"
      )
    }
  }

  # Zenodo/IODP records can bundle several files per DOI (e.g.
  # data_by_hole.zip and data_by_splice.zip on the same record) --
  # zenodo_fetch_dataset() already expands those into one row per file
  # found, same shape as pg_fetch_dataset(), so this loop is identical in
  # structure to fetch_and_harmonize() above, just pointed at Zenodo's
  # own depth-column convention (CSF-A depth scale, not PANGAEA's
  # "Depth sed") and its own synonym table.
  fetch_and_harmonize_zenodo <- function(doi, core_id) {
    append_log(sprintf("Fetching Zenodo record %s for %s...", doi, core_id))
    fetched <- tryCatch(zenodo_fetch_dataset(doi), error = function(e) {
      append_log(paste("ERROR fetching from Zenodo:", conditionMessage(e)))
      NULL
    })
    if (is.null(fetched) || nrow(fetched) == 0) {
      append_log(sprintf("No tabular data file found in Zenodo record %s.", doi))
      return(invisible())
    }
    for (i in seq_len(nrow(fetched))) {
      harmonize_and_store(
        fetched$data[[i]], core_id, fetched$doi[i],
        depth_col = "Top depth CSF-A (m)", synonym_map = load_zenodo_synonym_map(),
        synonym_csv_name = "zenodo_column_synonyms.csv"
      )
    }
  }

  observeEvent(input$search_btn, {
    req(rv$core_registry, input$active_core)
    row <- rv$core_registry[rv$core_registry$core_id == input$active_core, ]
    # Reuses the same measurement_filter control as harmonize/interpolate
    # and neighbour search -- when set, this is real additive cost (each
    # result has to be fetched to check its actual parameters, unlike
    # find_neighbours()'s equivalent filter, which reuses fetches it
    # makes anyway), hence the progress bar.
    required_parameters <- if (length(input$measurement_filter) > 0) input$measurement_filter else NULL
    filter_note <- if (!is.null(required_parameters)) {
      sprintf(", requiring one of: %s", paste(required_parameters, collapse = ", "))
    } else {
      ""
    }
    append_log(sprintf("Searching PANGAEA for %s%s...", input$active_core, filter_note))

    res <- withProgress(message = "Searching PANGAEA", value = 0, {
      tryCatch(
        pg_search_cores(
          row, count = 20, required_parameters = required_parameters,
          progress = function(i, n, doi) {
            setProgress(value = i / n, detail = sprintf("Checking result %d of %d...", i, n))
          }
        ),
        error = function(e) {
          append_log(paste("ERROR:", conditionMessage(e)))
          NULL
        }
      )
    })

    if (!is.null(res)) {
      rv$search_results <- res
      append_log(sprintf("Found %d candidate datasets.", nrow(res)))
    }
  })

  output$search_results <- DT::renderDataTable({
    DT::datatable(
      rv$search_results[, c("doi", "score", "size", "size_measure", "available_parameters", "citation")],
      selection = "multiple", options = list(pageLength = 5, scrollX = TRUE)
    )
  })

  observeEvent(input$fetch_selected_btn, {
    sel <- input$search_results_rows_selected
    if (length(sel) == 0) {
      append_log("No rows selected in the search results table -- click a row to select it, then try again.")
      return(invisible())
    }
    for (doi in rv$search_results$doi[sel]) {
      fetch_and_harmonize(doi, input$active_core)
    }
  })

  # Zenodo IODP records are per-expedition (e.g. "IODP Expedition 362
  # Elemental analysis (CHNS)"), not per-core like PANGAEA datasets, so a
  # core_id text search wouldn't match anything meaningful -- auto-fill
  # the expedition number when the active core's ID looks like an IODP
  # hole (e.g. "362-U1480A"), matching classify_core_source()'s own
  # IODP-ID convention, so the common case needs no manual typing.
  observeEvent(input$active_core, {
    req(input$active_core)
    expedition <- iodp_expedition_from_core_id(input$active_core)
    if (!is.na(expedition)) {
      updateTextInput(session, "zenodo_expedition", value = expedition)
    }
  })

  observeEvent(input$zenodo_search_btn, {
    req(input$active_core)
    expedition <- trimws(input$zenodo_expedition %||% "")
    if (!nzchar(expedition)) expedition <- NULL
    query <- trimws(input$zenodo_query %||% "")
    if (!nzchar(query)) query <- NULL
    if (is.null(expedition) && is.null(query)) {
      append_log(paste(
        "Enter an IODP expedition number and/or search terms to search",
        "Zenodo (expedition is auto-filled when the active core's ID",
        "looks like an IODP hole, e.g. \"362-U1480A\")."
      ))
      return(invisible())
    }
    append_log(sprintf(
      "Searching Zenodo (expedition: %s; terms: %s)...",
      expedition %||% "any", query %||% "none"
    ))
    res <- tryCatch(
      zenodo_search_iodp(query = query, expedition = expedition, count = 20),
      error = function(e) {
        append_log(paste("ERROR:", conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(res)) {
      rv$zenodo_search_results <- res
      append_log(sprintf("Found %d candidate Zenodo record(s).", nrow(res)))
    }
  })

  output$zenodo_search_results <- DT::renderDataTable({
    display <- rv$zenodo_search_results[, c("doi", "title", "n_files", "description")]
    display$doi <- doi_link(display$doi, base_url = "https://doi.org/")
    DT::datatable(
      display,
      selection = "multiple", options = list(pageLength = 5, scrollX = TRUE),
      escape = -2
    )
  })

  observeEvent(input$zenodo_fetch_selected_btn, {
    sel <- input$zenodo_search_results_rows_selected
    if (length(sel) == 0) {
      append_log("No rows selected in the Zenodo search results table -- click a row to select it, then try again.")
      return(invisible())
    }
    for (doi in rv$zenodo_search_results$doi[sel]) {
      fetch_and_harmonize_zenodo(doi, input$active_core)
    }
  })

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
    req(input$active_core)
    if (is.null(input$external_import_file)) {
      append_log("Choose a CSV file (matching the template) to import first.")
      return(invisible())
    }
    source_label <- trimws(input$external_import_source %||% "")
    if (!nzchar(source_label)) source_label <- NA_character_
    doi_label <- trimws(input$external_import_doi %||% "")
    if (!nzchar(doi_label)) doi_label <- NA_character_

    h <- tryCatch(
      import_external_data(
        input$external_import_file$datapath, core_id = input$active_core,
        source = source_label, doi = doi_label
      ),
      error = function(e) {
        append_log(paste("ERROR importing:", conditionMessage(e)))
        NULL
      }
    )
    if (is.null(h) || nrow(h) == 0) {
      return(invisible())
    }

    # Same optional measurement filter as harmonize_and_store() applies
    # to PANGAEA/Zenodo fetches, for consistency.
    n_before_filter <- nrow(h)
    h <- keep_selected_measurements(h, input$measurement_filter)
    if (nrow(h) > 0) {
      rv$param_long <- dplyr::bind_rows(rv$param_long, h)
      append_log(sprintf(
        "Imported %d rows for %s from %s.",
        nrow(h), input$active_core, input$external_import_file$name
      ))
    } else {
      append_log(sprintf(
        "0 rows from %s after restricting to selected measurements (%d rows available before filtering).",
        input$external_import_file$name, n_before_filter
      ))
    }
  })

  # find_neighbours() makes many sequential PANGAEA fetches (one per bbox
  # candidate) and can take a while. withProgress()/setProgress() push a
  # real progress bar to the browser *during* the call -- Shiny's official
  # mechanism for this, using a direct websocket write that bypasses the
  # normal reactive flush cycle -- so there's live feedback on exactly
  # which candidate is being checked, and periodic traffic on the
  # connection during the search rather than one long silence.
  #
  # (An ExtendedTask + future/promises background-process version of this
  # was tried first, to make the search fully non-blocking, but produced
  # an unexplained hang of the whole session partway through even with
  # globals passed explicitly -- reverted in favour of this simpler,
  # well-established synchronous-with-progress approach.)
  observeEvent(input$neighbour_btn, {
    req(rv$core_registry, input$active_core)
    if (length(input$neighbour_source_programs) == 0) {
      append_log("Select at least one source program before searching for neighbours.")
      return(invisible())
    }
    # Reuses the same measurement_filter control used to restrict
    # harmonize/interpolate -- one "measurements I care about" setting,
    # applied everywhere. Costs no extra PANGAEA requests (see
    # find_neighbours()'s required_parameters docs): every candidate's
    # data is already fetched to read its coordinates.
    required_parameters <- if (length(input$measurement_filter) > 0) input$measurement_filter else NULL
    filter_note <- if (!is.null(required_parameters)) {
      sprintf(", requiring one of: %s", paste(required_parameters, collapse = ", "))
    } else {
      ""
    }
    append_log(sprintf(
      "Searching neighbours within %s km of %s (programs: %s, up to %s candidates%s)...",
      input$radius_km, input$active_core,
      paste(input$neighbour_source_programs, collapse = ", "), input$neighbour_count, filter_note
    ))

    res <- withProgress(message = "Searching PANGAEA for neighbours", value = 0, {
      tryCatch(
        find_neighbours(
          rv$core_registry, input$active_core,
          radius_km = input$radius_km, count = input$neighbour_count,
          source_programs = input$neighbour_source_programs,
          required_parameters = required_parameters,
          progress = function(i, n, doi) {
            setProgress(value = i / n, detail = sprintf("Checking candidate %d of %d...", i, n))
          }
        ),
        error = function(e) {
          append_log(paste("ERROR:", conditionMessage(e)))
          NULL
        }
      )
    })

    if (!is.null(res)) {
      rv$neighbour_results <- res
      append_log(sprintf(
        "Found %d neighbouring datasets (sorted: GeoB/IODP first, most recent first).",
        nrow(res)
      ))
    }
  })

  output$neighbour_results <- DT::renderDataTable({
    display <- rv$neighbour_results[, c(
      "neighbour_core_id", "neighbour_doi", "distance_km", "source_program",
      "date_time", "available_parameters", "citation"
    )]
    display$neighbour_doi <- doi_link(display$neighbour_doi)
    DT::datatable(
      display,
      selection = "multiple", options = list(pageLength = 5, scrollX = TRUE, order = list()),
      # rownames add a display column before the data ones, so
      # neighbour_doi (2nd data column) is the 3rd rendered column;
      # escape = -3 leaves it as raw HTML (the <a> link) while every
      # other column still gets escaped as normal.
      escape = -3
    )
  })

  observeEvent(input$fetch_neighbours_btn, {
    sel <- input$neighbour_results_rows_selected
    if (length(sel) == 0) {
      append_log("No rows selected in the neighbours table -- click a row to select it, then try again.")
      return(invisible())
    }
    for (i in sel) {
      doi <- rv$neighbour_results$neighbour_doi[i]
      nb_core_id <- rv$neighbour_results$neighbour_core_id[i]
      if (is.na(nb_core_id)) nb_core_id <- doi
      fetch_and_harmonize(doi, nb_core_id)
    }
  })

  # Builds a SEPARATE depth-matched table (rv$depth_matched_long), rather
  # than appending into rv$param_long -- the raw/harmonized table stays a
  # complete, untouched record of everything retrieved, and this is a
  # purpose-built view aligned to the active core's own sample depths.
  #
  # Candidate measurements come from the active core's own harmonized data
  # PLUS any already-fetched neighbour cores found for it -- some cores
  # (e.g. GeoB18549-1) have no PANGAEA data of their own at all, so their
  # neighbours' measurements are the only depth series available to match
  # against their sample depths. report_depth_coverage() always returns a
  # row for every one of the active core's sample depths x every
  # parameter present in the candidate data, even when nothing could be
  # matched there (flagged via match_status), so missing coverage is
  # always visible rather than the target depth silently vanishing from
  # the output. The output's `core_id` still reflects the true
  # measurement source; a `target_core_id` column records which active
  # core it was matched against, so borrowed-from-a-neighbour rows are
  # never silently relabelled as native data.
  observeEvent(input$interp_btn, {
    req(rv$sample_registry, input$active_core)
    sample_depths <- rv$sample_registry$depth_cm[rv$sample_registry$core_id == input$active_core]
    if (length(sample_depths) == 0 || all(is.na(sample_depths))) {
      append_log("No sample depths on file for the active core -- interpolation needs real depths to match against.")
      return(invisible())
    }

    # Candidate selection and every message below live in tested helpers
    # (interp_candidate_data(), depth_coverage_messages()); the observer
    # only wires inputs to them. age_model rows always survive the filter.
    candidates <- interp_candidate_data(
      rv$param_long, rv$neighbour_results, input$active_core, input$measurement_filter
    )
    if (!is.null(candidates$message)) {
      append_log(candidates$message)
      return(invisible())
    }

    report <- tryCatch(
      report_depth_coverage(
        candidates$data, sample_depths,
        target_core_id = input$active_core, sample_registry = rv$sample_registry
      ),
      error = function(e) {
        append_log(paste("ERROR interpolating:", conditionMessage(e)))
        NULL
      }
    )
    if (is.null(report)) {
      return(invisible())
    }

    report$target_core_id <- input$active_core
    rv$depth_matched_long <- dplyr::bind_rows(
      rv$depth_matched_long[rv$depth_matched_long$target_core_id != input$active_core, ],
      report
    )

    for (msg in depth_coverage_messages(report, input$active_core, sample_depths)) append_log(msg)
  })

  # An age model isn't shown as its own measurement row here -- it's
  # promoted to an age_ka column applied to every other row, so it's
  # always visible regardless of which parameter_name(s) the user has
  # selected in measurement_filter (unlike TOC/TN/etc, it's core context,
  # not a measurement you'd opt in/out of). Uses the same
  # core_age_model()/depth_to_age() machinery as cross-core neighbour
  # matching, just applied within each row's own core.
  output$param_table <- DT::renderDataTable({
    dt <- DT::datatable(attach_age_ka(rv$param_long), options = list(pageLength = 15, scrollX = TRUE))
    DT::formatStyle(dt, "is_interpolated",
      backgroundColor = DT::styleEqual(c(TRUE, FALSE), c("#fff3cd", "white"))
    )
  })

  # Same rationale as output$param_table above: age_model doesn't appear
  # as its own parameter row here either -- it's the age_ka column every
  # other row already carries (matched_age_ka, renamed for display), so a
  # separate "age_model" row would just be a redundant, less useful
  # version of the same information.
  output$depth_matched_table <- DT::renderDataTable({
    display <- rv$depth_matched_long[rv$depth_matched_long$parameter_name != "age_model", ]
    display <- dplyr::rename(display, age_ka = "matched_age_ka")
    display <- dplyr::relocate(display, "age_ka", .after = "depth_cm")

    dt <- DT::datatable(display, options = list(pageLength = 15, scrollX = TRUE))
    DT::formatStyle(dt, "match_status",
      backgroundColor = DT::styleEqual(
        c(
          "measured", "interpolated", "no data in that depth range", "gap too large to interpolate",
          "no age model for cross-core match", "target core has no age_ka for this sample depth"
        ),
        c("white", "#fff3cd", "#f8d7da", "#f8d7da", "#f8d7da", "#f8d7da")
      )
    )
  })

  # Quick-reference list of which measurement types have been retrieved so
  # far, for visualization/orientation before diving into the full table.
  output$measurement_summary <- DT::renderDataTable({
    DT::datatable(summarise_measurements(rv$param_long), options = list(pageLength = 20, scrollX = TRUE))
  })

  # Reactive (no button needed -- cheap to compute), always reflecting the
  # active core's current checklist status; see check_target_measurements().
  # Falls back to already-fetched neighbour cores' data (found via a prior
  # "Find neighbours" search for this active core) for parameters the core
  # has none of its own -- flagged distinctly ("found in neighbour"), not
  # conflated with the core's own real data.
  target_checklist <- reactive({
    req(input$active_core)
    neighbour_ids <- unique(rv$neighbour_results$neighbour_core_id[
      rv$neighbour_results$core_id == input$active_core
    ])
    check_target_measurements(rv$param_long, input$active_core, neighbour_core_ids = neighbour_ids)
  })

  output$target_checklist_summary <- renderText({
    checklist <- target_checklist()
    missing <- checklist$parameter_name[checklist$status == "missing"]
    in_neighbour <- checklist$parameter_name[checklist$status == "found in neighbour"]
    parts <- character()
    if (length(in_neighbour) > 0) {
      parts <- c(parts, sprintf("found only in a neighbour core: %s", paste(in_neighbour, collapse = ", ")))
    }
    if (length(missing) > 0) {
      parts <- c(parts, sprintf("missing entirely: %s", paste(missing, collapse = ", ")))
    }
    if (length(parts) == 0) {
      sprintf("All %d target measurement(s) found for %s.", nrow(checklist), input$active_core)
    } else {
      sprintf("For %s -- %s.", input$active_core, paste(parts, collapse = "; "))
    }
  })

  output$target_checklist_table <- DT::renderDataTable({
    dt <- DT::datatable(target_checklist(), options = list(pageLength = 20, scrollX = TRUE))
    DT::formatStyle(dt, "status",
      backgroundColor = DT::styleEqual(
        c("found", "found in neighbour", "missing"), c("white", "#fff3cd", "#f8d7da")
      )
    )
  })

  observeEvent(input$clear_data_btn, {
    rv$param_long <- empty_param_long()
    # depth_matched_long is derived entirely from param_long (via
    # report_depth_coverage()) -- clearing one without the other would
    # leave stale depth-matched rows pointing at measurements that no
    # longer exist.
    rv$depth_matched_long <- empty_depth_matched_long()
    append_log("Cleared harmonized and depth-matched data (fetched raw data on PANGAEA is untouched; re-fetch to start again).")
  })

  output$log <- renderText({
    paste(rv$log, collapse = "\n")
  })

  output$export_btn <- downloadHandler(
    filename = function() export_filename(input$active_core, "harmonized"),
    content = function(file) write_harmonized_tsv(rv$param_long, file)
  )

  output$export_depth_matched_btn <- downloadHandler(
    filename = function() export_filename(input$active_core, "depth_matched"),
    content = function(file) write_harmonized_tsv(rv$depth_matched_long, file)
  )

  # Auto-load the bundled example on startup.
  observe({
    if (is.null(rv$core_registry)) {
      default_path <- system.file("extdata", "example_core_metadata.tsv", package = "sedacore")
      isolate(tryCatch(
        load_cores(default_path),
        error = function(e) append_log(paste("ERROR:", conditionMessage(e)))
      ))
    }
  })
}

shinyApp(ui, server)
