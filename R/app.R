#' Launch the sedacore Shiny app
#'
#' A thin wrapper around [shiny::runApp()]; all business logic lives in
#' the already-tested `R/*.R` functions (map/search/retrieve/harmonize/
#' neighbour-search/interpolate/export) -- `inst/app/app.R` only wires
#' reactives to those functions.
#'
#' @param ... Passed to [shiny::runApp()].
#' @return Does not return; launches the app.
#' @export
run_sedacore_app <- function(...) {
  app_dir <- system.file("app", package = "sedacore")
  if (app_dir == "") {
    rlang::abort("Could not find the sedacore Shiny app; is the package installed correctly?")
  }
  shiny::runApp(app_dir, ...)
}
