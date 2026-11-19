# Module: the measurement filter ("the parameters I care about").
#
# One control, applied everywhere -- harmonizing, interpolating, and (as
# required_parameters) PANGAEA and neighbour search. Returns the selection
# as a reactive.

#' @noRd
mod_filter_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$h4("Measurement filter"),
    selectizeInput(ns("measurement_filter"), "Restrict to selected measurements",
      choices = NULL, multiple = TRUE,
      options = list(placeholder = "(leave empty = harmonize/interpolate everything)")
    ),
    actionLink(ns("use_target_checklist_btn"), "Use target checklist"),
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
    )
  )
}

#' @noRd
mod_filter_server <- function(id, state) {
  moduleServer(id, function(input, output, session) {
    # Seed the choices up front, then grow them with anything new that
    # shows up via harmonize's fallback path (a raw column not yet in
    # either synonym table), never dropping the user's current selection.
    # "age_model" stays selectable: selecting it still means something
    # (requiring it when filtering neighbours by available_parameters),
    # even though it is never *excluded* from storage or the harmonized
    # table by this control (see keep_selected_measurements()).
    updateSelectizeInput(session, "measurement_filter", choices = known_parameter_names())

    observeEvent(state$param_long, {
      selected <- isolate(input$measurement_filter)
      updateSelectizeInput(session, "measurement_filter",
        choices = sort(unique(c(known_parameter_names(), state$param_long$parameter_name))),
        selected = selected
      )
    })

    observeEvent(input$use_target_checklist_btn, {
      updateSelectizeInput(session, "measurement_filter",
        choices = sort(unique(c(known_parameter_names(), state$param_long$parameter_name, target_measurement_names()))),
        selected = target_measurement_names()
      )
    })

    reactive(input$measurement_filter)
  })
}
