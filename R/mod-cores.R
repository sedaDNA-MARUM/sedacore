# Module: core list, active core, and the map.
#
# Owns the sample/core registries (loaded from an uploaded TSV, or the
# bundled example on startup), the "active core" selector, and the map;
# clicking a map marker selects that core. Returns the active core as a
# reactive for every other module to use.

#' @noRd
mod_cores_sidebar_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fileInput(ns("core_file"), "Core metadata TSV",
      accept = c(".tsv", ".txt"),
      placeholder = "(bundled example used if empty)"
    ),
    actionButton(ns("load_btn"), "Load / reload cores", icon = icon("upload")),
    tags$hr(),
    selectInput(ns("active_core"), "Active core", choices = NULL)
  )
}

#' @noRd
mod_cores_map_ui <- function(id) {
  leaflet::leafletOutput(NS(id, "map"), height = 600)
}

#' @noRd
mod_cores_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    load_cores <- function(path) {
      registries <- load_core_registries(path)
      cr <- registries$core_registry
      state$sample_registry <- registries$sample_registry
      state$core_registry <- cr
      updateSelectInput(session, "active_core", choices = cr$core_id, selected = cr$core_id[1])
      log_message(state, sprintf("Loaded %d cores from %s.", nrow(cr), basename(path)))
    }

    observeEvent(input$load_btn,
      {
        path <- if (!is.null(input$core_file)) {
          input$core_file$datapath
        } else {
          system.file("extdata", "example_core_metadata.tsv", package = "sedacore")
        }
        tryCatch(
          load_cores(path),
          error = function(e) log_message(state, paste("ERROR loading cores:", conditionMessage(e)))
        )
      },
      ignoreNULL = FALSE
    )

    output$map <- leaflet::renderLeaflet({
      req(state$core_registry)
      leaflet_core_map(state$core_registry, neighbours = state$neighbour_results)
    })

    observeEvent(input$map_marker_click, {
      click <- input$map_marker_click
      if (!is.null(click$id) && !startsWith(click$id, "neighbour-")) {
        updateSelectInput(session, "active_core", selected = click$id)
      }
    })

    # Auto-load the bundled example on startup.
    observe({
      if (is.null(state$core_registry)) {
        default_path <- system.file("extdata", "example_core_metadata.tsv", package = "sedacore")
        isolate(tryCatch(
          load_cores(default_path),
          error = function(e) log_message(state, paste("ERROR:", conditionMessage(e)))
        ))
      }
    })

    reactive(input$active_core)
  })
}
