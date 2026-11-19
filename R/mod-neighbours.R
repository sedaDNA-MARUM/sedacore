# Module: neighbouring cores of the active core.
#
# find_neighbours() makes many sequential PANGAEA fetches (one per bbox
# candidate) and can take a while. withProgress()/setProgress() push a real
# progress bar to the browser *during* the call -- Shiny's official
# mechanism for this -- so there is live feedback on which candidate is
# being checked, and periodic traffic on the connection rather than one
# long silence. (An ExtendedTask + future/promises version, to make the
# search fully non-blocking, was tried first and hung the whole session
# partway through even with globals passed explicitly -- reverted in favour
# of this simpler synchronous-with-progress approach.)

#' @noRd
mod_neighbours_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$h4("Neighbouring cores"),
    numericInput(ns("radius_km"), "Radius (km)", value = 50, min = 1, max = 2000),
    numericInput(ns("neighbour_count"), "Max candidates to check", value = 50, min = 5, max = 300, step = 10),
    helpText(
      "PANGAEA's search relevance isn't distance-based, so in densely-",
      "sampled regions (e.g. well-studied MARUM survey areas) real nearby",
      "cores can rank surprisingly low -- raise this if you suspect",
      "known-close cores are missing from the results. Higher values",
      "take longer (each candidate is a separate PANGAEA request)."
    ),
    checkboxGroupInput(ns("neighbour_source_programs"), "Source program(s)",
      choices = c("GeoB", "IODP", "ODP", "DSDP", "other"),
      selected = c("GeoB", "IODP", "ODP", "DSDP", "other"),
      inline = TRUE
    ),
    actionButton(ns("neighbour_btn"), "Find neighbours"),
    DT::dataTableOutput(ns("neighbour_results")),
    actionButton(ns("fetch_neighbours_btn"), "Fetch + harmonize selected neighbours", class = "btn-primary")
  )
}

#' @noRd
mod_neighbours_server <- function(id, state, active_core, measurement_filter, fetchers,
                                  find = find_neighbours) {
  moduleServer(id, function(input, output, session) {
    log_msg <- function(msg) log_message(state, msg)

    observeEvent(input$neighbour_btn, {
      req(state$core_registry, active_core())
      if (length(input$neighbour_source_programs) == 0) {
        log_msg("Select at least one source program before searching for neighbours.")
        return(invisible())
      }
      # The same measurement_filter control used everywhere: one
      # "measurements I care about" setting. Costs no extra PANGAEA
      # requests here (see find_neighbours()'s required_parameters docs):
      # every candidate's data is already fetched to read its coordinates.
      required_parameters <- if (length(measurement_filter()) > 0) measurement_filter() else NULL
      filter_note <- if (!is.null(required_parameters)) {
        sprintf(", requiring one of: %s", paste(required_parameters, collapse = ", "))
      } else {
        ""
      }
      log_msg(sprintf(
        "Searching neighbours within %s km of %s (programs: %s, up to %s candidates%s)...",
        input$radius_km, active_core(),
        paste(input$neighbour_source_programs, collapse = ", "), input$neighbour_count, filter_note
      ))

      res <- withProgress(message = "Searching PANGAEA for neighbours", value = 0, {
        tryCatch(
          find(
            state$core_registry, active_core(),
            radius_km = input$radius_km, count = input$neighbour_count,
            source_programs = input$neighbour_source_programs,
            required_parameters = required_parameters,
            progress = function(i, n, doi) {
              setProgress(value = i / n, detail = sprintf("Checking candidate %d of %d...", i, n))
            }
          ),
          error = function(e) {
            log_msg(paste("ERROR:", conditionMessage(e)))
            NULL
          }
        )
      })

      if (!is.null(res)) {
        state$neighbour_results <- res
        log_msg(sprintf(
          "Found %d neighbouring datasets (sorted: GeoB/IODP first, most recent first).",
          nrow(res)
        ))
      }
    })

    output$neighbour_results <- DT::renderDataTable({
      display <- state$neighbour_results[, c(
        "neighbour_core_id", "neighbour_doi", "distance_km", "source_program",
        "date_time", "available_parameters", "citation"
      )]
      display$neighbour_doi <- doi_link(display$neighbour_doi)
      # rownames add a display column before the data ones, so
      # neighbour_doi (2nd data column) is the 3rd rendered column;
      # escape = -3 leaves it as raw HTML (the <a> link) while every other
      # column still gets escaped as normal.
      DT::datatable(
        display,
        selection = "multiple", options = list(pageLength = 5, scrollX = TRUE, order = list()),
        escape = -3
      )
    })

    observeEvent(input$fetch_neighbours_btn, {
      sel <- input$neighbour_results_rows_selected
      if (length(sel) == 0) {
        log_msg("No rows selected in the neighbours table -- click a row to select it, then try again.")
        return(invisible())
      }
      for (i in sel) {
        doi <- state$neighbour_results$neighbour_doi[i]
        nb_core_id <- state$neighbour_results$neighbour_core_id[i]
        if (is.na(nb_core_id)) nb_core_id <- doi
        fetchers$pangaea(doi, nb_core_id)
      }
    })
  })
}
