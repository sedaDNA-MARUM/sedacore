#' @noRd
app_ui <- function() {
  fluidPage(
    tags$head(
      tags$link(rel = "preconnect", href = "https://fonts.googleapis.com"),
      tags$link(
        rel = "stylesheet",
        href = "https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,600&display=swap"
      ),
      # The sidebar is much taller than the main panel -- without this,
      # scrolling down to reach the lower sidebar controls scrolls the
      # map/tables out of view entirely. Giving the sidebar its own internal
      # scroll (pinned to the viewport) keeps the page itself only as tall
      # as the main panel, so the map/tables never scroll out of view.
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
          " \u2014 core contextual data explorer",
          style = "font-weight: 300; color: #666; font-size: 1.6rem; margin-left: 6px;"
        )
      ),
      windowTitle = "sedacore"
    ),
    sidebarLayout(
      sidebarPanel(
        width = 4,
        mod_cores_sidebar_ui("cores"), tags$hr(),
        mod_filter_ui("filter"), tags$hr(),
        mod_retrieve_ui("retrieve"), tags$hr(),
        mod_neighbours_ui("neighbours"), tags$hr(),
        mod_integrate_ui("integrate"), tags$hr(),
        mod_outputs_sidebar_ui("outputs")
      ),
      mainPanel(
        width = 8,
        tabsetPanel(
          tabPanel("Map", mod_cores_map_ui("cores")),
          tabPanel("Measurements", mod_outputs_measurements_ui("outputs")),
          tabPanel("Target checklist", mod_outputs_checklist_ui("outputs")),
          tabPanel("Harmonized data", mod_outputs_harmonized_ui("outputs")),
          tabPanel("Depth/age-matched data", mod_integrate_table_ui("integrate")),
          tabPanel("Log", mod_outputs_log_ui("outputs"))
        )
      )
    )
  )
}
