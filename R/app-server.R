#' @noRd
app_server <- function(input, output, session) {
  state <- new_app_state()

  # Two modules own the controls everything else depends on; the rest only
  # ever see them as reactives, never as raw inputs.
  active_core <- mod_cores_server("cores", state)
  measurement_filter <- mod_filter_server("filter", state)
  fetchers <- make_fetchers(state, measurement_filter)

  mod_retrieve_server("retrieve", state, active_core, measurement_filter, fetchers)
  mod_neighbours_server("neighbours", state, active_core, measurement_filter, fetchers)
  mod_integrate_server("integrate", state, active_core, measurement_filter)
  mod_outputs_server("outputs", state, active_core)
}

#' The sedacore Shiny app
#'
#' Returns the app object without running it (so it can be run with
#' [shiny::runApp()], deployed, or driven by `shiny::testServer()`).
#' [run_sedacore_app()] is the usual way to launch it.
#'
#' @return A `shiny.appobj`.
#' @export
sedacore_app <- function() {
  shiny::shinyApp(app_ui(), app_server)
}
