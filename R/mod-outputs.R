# Module: read-only views of what has been retrieved, plus export and clear.
#
# Measurements summary, target checklist, harmonized table, log; the two
# TSV exports; and the button that clears harmonized + depth-matched data.

#' @noRd
mod_outputs_sidebar_ui <- function(id) {
  ns <- NS(id)
  tagList(
    downloadButton(ns("export_btn"), "Export harmonized TSV"),
    downloadButton(ns("export_depth_matched_btn"), "Export depth-matched TSV"),
    actionButton(ns("clear_data_btn"), "Clear harmonized + depth-matched data", class = "btn-danger")
  )
}

#' @noRd
mod_outputs_measurements_ui <- function(id) DT::dataTableOutput(NS(id, "measurement_summary"))

#' @noRd
mod_outputs_harmonized_ui <- function(id) DT::dataTableOutput(NS(id, "param_table"))

#' @noRd
mod_outputs_log_ui <- function(id) verbatimTextOutput(NS(id, "log"))

#' @noRd
mod_outputs_checklist_ui <- function(id) {
  ns <- NS(id)
  tagList(
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
    textOutput(ns("target_checklist_summary")),
    tags$br(),
    DT::dataTableOutput(ns("target_checklist_table"))
  )
}

#' @noRd
mod_outputs_server <- function(id, state, active_core) {
  moduleServer(id, function(input, output, session) {
    output$measurement_summary <- DT::renderDataTable({
      DT::datatable(summarise_measurements(state$param_long), options = list(pageLength = 20, scrollX = TRUE))
    })

    # An age model isn't shown as its own measurement row: it is promoted to
    # an age_ka column on every other row, so it's always visible whichever
    # parameter_name(s) are selected in the measurement filter (see
    # attach_age_ka()).
    output$param_table <- DT::renderDataTable({
      dt <- DT::datatable(attach_age_ka(state$param_long), options = list(pageLength = 15, scrollX = TRUE))
      DT::formatStyle(dt, "is_interpolated",
        backgroundColor = DT::styleEqual(c(TRUE, FALSE), c("#fff3cd", "white"))
      )
    })

    # Reactive (no button needed -- cheap to compute), always reflecting the
    # active core's current checklist status; see check_target_measurements().
    # Falls back to already-fetched neighbours' data (found via a prior
    # "Find neighbours" search) for parameters the core has none of its own
    # -- flagged distinctly ("found in neighbour"), never conflated with the
    # core's own real data.
    target_checklist <- reactive({
      req(active_core())
      neighbour_ids <- unique(state$neighbour_results$neighbour_core_id[
        state$neighbour_results$core_id == active_core()
      ])
      check_target_measurements(state$param_long, active_core(), neighbour_core_ids = neighbour_ids)
    })

    output$target_checklist_summary <- renderText({
      checklist_summary(target_checklist(), active_core())
    })

    output$target_checklist_table <- DT::renderDataTable({
      dt <- DT::datatable(target_checklist(), options = list(pageLength = 20, scrollX = TRUE))
      DT::formatStyle(dt, "status",
        backgroundColor = DT::styleEqual(
          c("found", "found in neighbour", "missing"), c("white", "#fff3cd", "#f8d7da")
        )
      )
    })

    output$log <- renderText(paste(state$log, collapse = "\n"))

    observeEvent(input$clear_data_btn, {
      state$param_long <- empty_param_long()
      # depth_matched_long is derived entirely from param_long (via
      # report_depth_coverage()) -- clearing one without the other would
      # leave stale depth-matched rows pointing at measurements that no
      # longer exist.
      state$depth_matched_long <- empty_depth_matched_long()
      log_message(state, "Cleared harmonized and depth-matched data (fetched raw data on PANGAEA is untouched; re-fetch to start again).")
    })

    output$export_btn <- downloadHandler(
      filename = function() export_filename(active_core(), "harmonized"),
      content = function(file) write_harmonized_tsv(state$param_long, file)
    )

    output$export_depth_matched_btn <- downloadHandler(
      filename = function() export_filename(active_core(), "depth_matched"),
      content = function(file) write_harmonized_tsv(state$depth_matched_long, file)
    )
  })
}
