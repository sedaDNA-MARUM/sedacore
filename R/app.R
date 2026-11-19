#' Launch the sedacore Shiny app
#'
#' A thin wrapper around [shiny::runApp()] on [sedacore_app()]. All the
#' logic lives in tested package functions and modules (`R/app-helpers.R`,
#' `R/mod-*.R`); nothing is hidden in a script.
#'
#' @param ... Passed to [shiny::runApp()] (e.g. `port`, `host`,
#'   `launch.browser`).
#' @return Does not return; launches the app.
#' @export
run_sedacore_app <- function(...) {
  shiny::runApp(sedacore_app(), ...)
}
