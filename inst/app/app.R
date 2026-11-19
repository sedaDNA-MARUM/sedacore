# The app itself lives in the package (R/app-ui.R, R/app-server.R and the
# R/mod-*.R modules) so that it can be tested; this file only exists so
# shiny::runApp(system.file("app", package = "sedacore")) keeps working.
sedacore::sedacore_app()
