library(shiny)
library(leaflet)
library(sedacore)

empty_param_long <- function() {
  new_param_long(tibble::tibble(
    core_id = character(), pangaea_doi = character(), depth_cm = numeric(),
    parameter_group = character(), parameter_name = character(),
    species = character(), value = numeric(), unit = character(),
    raw_column_name = character(), is_interpolated = logical(),
    interp_gap_cm = numeric(), retrieved_at = as.POSIXct(character())
  ))
}

# Like empty_param_long(), plus target_core_id: which active core's sample
# depths a row was matched against. Usually equal to the row's own
# core_id, but differs when the row's measurements were *borrowed* from a
# neighbour core (see the interp_btn handler) -- kept as an explicit,
# separate column so that's always visible rather than implied.
empty_depth_matched_long <- function() {
  x <- empty_param_long()
  x$target_core_id <- character()
  x
}

# Turns a raw PANGAEA DOI (e.g. "10.1594/PANGAEA.228750") into a
# clickable link to its landing page, for DT columns rendered with
# escape = FALSE for that column. NA-safe (NA stays NA, not "NA" text).
doi_link <- function(doi) {
  ifelse(is.na(doi), NA_character_,
    sprintf(
      '<a href="https://doi.pangaea.de/%s" target="_blank" rel="noopener noreferrer">%s</a>',
      doi, doi
    )
  )
}

empty_search_results <- function() {
  tibble::tibble(
    core_id = character(), score = numeric(), doi = character(),
    size = numeric(), size_measure = character(), citation = character()
  )
}

empty_neighbour_results <- function() {
  new_neighbour_registry(tibble::tibble(
    core_id = character(), neighbour_doi = character(),
    neighbour_core_id = character(), latitude = numeric(), longitude = numeric(),
    distance_km = numeric(), source_program = character(),
    date_time = as.POSIXct(character()), citation = character(),
    available_parameters = character(), included_in_export = logical()
  ))
}

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
        "happens to contain, and also requires neighbour search results to",
        "have at least one of these parameters (e.g. only show neighbours",
        "with TOC data). Leave empty for the default (keep/show everything)."
      ),
      tags$hr(),
      tags$h4("Retrieve from PANGAEA"),
      actionButton("search_btn", "Search PANGAEA for active core"),
      DT::dataTableOutput("search_results"),
      actionButton("fetch_selected_btn", "Fetch + harmonize selected", class = "btn-primary"),
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
        tabPanel("Depth-matched data", DT::dataTableOutput("depth_matched_table")),
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
    neighbour_results = empty_neighbour_results(),
    log = character()
  )

  append_log <- function(msg) {
    rv$log <- c(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg), rv$log)
  }

  # Seed the measurement-filter choices from the full known parameter_name
  # vocabulary (column_synonyms.csv) up front, then grow it with anything
  # new that shows up via harmonize's fallback path (a raw column not yet
  # in that table) so those stay selectable too, without ever dropping the
  # user's current selection.
  updateSelectizeInput(session, "measurement_filter",
    choices = sort(unique(load_synonym_map()$parameter_name))
  )
  observeEvent(rv$param_long, {
    known <- isolate(input$measurement_filter)
    all_choices <- sort(unique(c(load_synonym_map()$parameter_name, rv$param_long$parameter_name)))
    updateSelectizeInput(session, "measurement_filter",
      choices = all_choices, selected = known
    )
  })

  observeEvent(input$use_target_checklist_btn, {
    updateSelectizeInput(session, "measurement_filter",
      choices = sort(unique(c(load_synonym_map()$parameter_name, rv$param_long$parameter_name, load_target_measurements()$parameter_name))),
      selected = load_target_measurements()$parameter_name
    )
  })

  load_cores <- function(path) {
    raw <- readr::read_tsv(path, show_col_types = FALSE)
    sr <- as_sample_registry(raw)
    cr <- collapse_to_core_registry(sr)
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
      h <- tryCatch(
        harmonize_columns(fetched$data[[i]], core_id = core_id, doi = fetched$doi[i], unmatched = "warn"),
        error = function(e) {
          append_log(paste("ERROR harmonizing:", conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(h) && nrow(h) > 0) {
        fallback <- attr(h, "fallback_columns")
        dropped <- attr(h, "dropped_columns")
        n_before_filter <- nrow(h)
        # Optional user-set filter (see "Measurement filter" in the
        # sidebar) -- applied *after* harmonize_columns() so nothing is
        # ever lost during parsing itself, only at the point of deciding
        # what to actually keep.
        if (length(input$measurement_filter) > 0) {
          h <- h[h$parameter_name %in% input$measurement_filter, ]
        }
        if (nrow(h) > 0) {
          rv$param_long <- dplyr::bind_rows(rv$param_long, h)
          append_log(sprintf("Harmonized %d rows from %s.", nrow(h), doi))
        } else {
          append_log(sprintf(
            "0 rows from %s after restricting to selected measurements (%d rows available before filtering).",
            doi, n_before_filter
          ))
        }
        if (length(fallback) > 0) {
          append_log(sprintf(
            "Not in column_synonyms.csv, imported as-is: %s",
            paste(fallback, collapse = ", ")
          ))
        }
        if (length(dropped) > 0) {
          append_log(sprintf(
            "Excluded (non-numeric, e.g. free-text comment fields): %s",
            paste(dropped, collapse = ", ")
          ))
        }
      }
    }
  }

  observeEvent(input$search_btn, {
    req(rv$core_registry, input$active_core)
    row <- rv$core_registry[rv$core_registry$core_id == input$active_core, ]
    append_log(sprintf("Searching PANGAEA for %s...", input$active_core))
    res <- tryCatch(pg_search_cores(row, count = 20), error = function(e) {
      append_log(paste("ERROR:", conditionMessage(e)))
      NULL
    })
    if (!is.null(res)) {
      rv$search_results <- res
      append_log(sprintf("Found %d candidate datasets.", nrow(res)))
    }
  })

  output$search_results <- DT::renderDataTable({
    DT::datatable(
      rv$search_results[, c("doi", "score", "size", "size_measure", "citation")],
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

    neighbour_ids <- unique(rv$neighbour_results$neighbour_core_id[
      rv$neighbour_results$core_id == input$active_core
    ])
    neighbour_ids <- neighbour_ids[!is.na(neighbour_ids)]
    candidate_core_ids <- union(input$active_core, neighbour_ids)
    current <- rv$param_long[rv$param_long$core_id %in% candidate_core_ids, ]
    # Independent of any filtering already applied at harmonize time --
    # covers data fetched before the filter was set, and lets the filter
    # be changed freely without re-fetching.
    if (length(input$measurement_filter) > 0) {
      current <- current[current$parameter_name %in% input$measurement_filter, ]
    }
    if (nrow(current) == 0) {
      msg <- "No harmonized data yet for the active core"
      if (length(neighbour_ids) > 0) msg <- paste0(msg, " or its fetched neighbours")
      if (length(input$measurement_filter) > 0) msg <- paste0(msg, " matching the selected measurement filter")
      append_log(paste0(msg, "; fetch some first."))
      return(invisible())
    }

    report <- tryCatch(report_depth_coverage(current, sample_depths), error = function(e) {
      append_log(paste("ERROR interpolating:", conditionMessage(e)))
      NULL
    })
    if (is.null(report)) {
      return(invisible())
    }

    report$target_core_id <- input$active_core
    rv$depth_matched_long <- dplyr::bind_rows(
      rv$depth_matched_long[rv$depth_matched_long$target_core_id != input$active_core, ],
      report
    )

    n_measured <- sum(report$match_status == "measured")
    n_interp <- sum(report$match_status == "interpolated")
    n_no_data <- sum(report$match_status == "no data in that depth range")
    n_gap <- sum(report$match_status == "gap too large to interpolate")
    borrowed_ids <- unique(report$core_id[!is.na(report$core_id) & report$core_id != input$active_core])
    source_note <- if (length(borrowed_ids) > 0) {
      sprintf(" (some borrowed from neighbour core(s): %s)", paste(borrowed_ids, collapse = ", "))
    } else {
      ""
    }
    append_log(sprintf(
      "Depth coverage for %s's %d sample depth(s) x %d parameter(s): %d measured, %d interpolated, %d out of retrieved range, %d gaps too large%s.",
      input$active_core, length(unique(sample_depths)), length(unique(report$parameter_name)),
      n_measured, n_interp, n_no_data, n_gap, source_note
    ))

    if (n_no_data > 0) {
      no_data_summary <- report[report$match_status == "no data in that depth range", ]
      no_data_summary <- unique(no_data_summary[, c("parameter_group", "parameter_name", "measured_depth_min", "measured_depth_max")])
      details <- sprintf(
        "%s (%s): retrieved %.0f-%.0f cm",
        no_data_summary$parameter_name, no_data_summary$parameter_group,
        no_data_summary$measured_depth_min, no_data_summary$measured_depth_max
      )
      append_log(sprintf(
        "Depth mismatch: your samples are at %.0f-%.0f cm, but retrieved data doesn't cover that range for: %s.",
        min(sample_depths, na.rm = TRUE), max(sample_depths, na.rm = TRUE), paste(details, collapse = "; ")
      ))
    }
  })

  output$param_table <- DT::renderDataTable({
    dt <- DT::datatable(rv$param_long, options = list(pageLength = 15, scrollX = TRUE))
    DT::formatStyle(dt, "is_interpolated",
      backgroundColor = DT::styleEqual(c(TRUE, FALSE), c("#fff3cd", "white"))
    )
  })

  output$depth_matched_table <- DT::renderDataTable({
    dt <- DT::datatable(rv$depth_matched_long, options = list(pageLength = 15, scrollX = TRUE))
    DT::formatStyle(dt, "match_status",
      backgroundColor = DT::styleEqual(
        c("measured", "interpolated", "no data in that depth range", "gap too large to interpolate"),
        c("white", "#fff3cd", "#f8d7da", "#f8d7da")
      )
    )
  })

  # Quick-reference list of which measurement types have been retrieved so
  # far, for visualization/orientation before diving into the full table.
  output$measurement_summary <- DT::renderDataTable({
    summary_tbl <- rv$param_long |>
      dplyr::group_by(parameter_group, parameter_name, unit) |>
      dplyr::summarise(
        n_cores = dplyr::n_distinct(core_id),
        cores = paste(sort(unique(core_id)), collapse = ", "),
        n_measured = sum(!is_interpolated),
        n_interpolated = sum(is_interpolated),
        .groups = "drop"
      ) |>
      dplyr::arrange(parameter_group, parameter_name)
    DT::datatable(summary_tbl, options = list(pageLength = 20, scrollX = TRUE))
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
    filename = function() {
      core_part <- if (!is.null(input$active_core) && nzchar(input$active_core)) {
        gsub("[^A-Za-z0-9_-]+", "_", input$active_core)
      } else {
        "export"
      }
      sprintf("sedacore_%s_%s.tsv", core_part, format(Sys.Date()))
    },
    content = function(file) write_harmonized_tsv(rv$param_long, file)
  )

  output$export_depth_matched_btn <- downloadHandler(
    filename = function() {
      core_part <- if (!is.null(input$active_core) && nzchar(input$active_core)) {
        gsub("[^A-Za-z0-9_-]+", "_", input$active_core)
      } else {
        "export"
      }
      sprintf("sedacore_depth_matched_%s_%s.tsv", core_part, format(Sys.Date()))
    },
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
