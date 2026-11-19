# Module: depth/age matching of the active core's sample depths.
#
# Builds a SEPARATE table (state$depth_matched_long) rather than appending
# into state$param_long: the harmonized table stays a complete, untouched
# record of everything retrieved, and this is a purpose-built view aligned
# to the active core's own sample depths. Candidates are the core's own data
# PLUS any already-fetched neighbours (some cores have no data of their own,
# so their neighbours' series are all there is to match against). The
# output's `core_id` still reflects the true measurement source; a
# `target_core_id` column records which active core a row was matched
# against, so borrowed rows are never silently relabelled as native data.
# All selection and message logic lives in tested helpers
# (interp_candidate_data(), depth_coverage_messages()).

#' @noRd
mod_integrate_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tags$h4("Interpolate"),
    actionButton(ns("interp_btn"), "Interpolate to active core's sample depths")
  )
}

#' @noRd
mod_integrate_table_ui <- function(id) {
  DT::dataTableOutput(NS(id, "depth_matched_table"))
}

#' @noRd
mod_integrate_server <- function(id, state, active_core, measurement_filter) {
  moduleServer(id, function(input, output, session) {
    log_msg <- function(msg) log_message(state, msg)

    observeEvent(input$interp_btn, {
      req(state$sample_registry, active_core())
      sample_depths <- state$sample_registry$depth_cm[state$sample_registry$core_id == active_core()]
      if (length(sample_depths) == 0 || all(is.na(sample_depths))) {
        log_msg("No sample depths on file for the active core -- interpolation needs real depths to match against.")
        return(invisible())
      }

      candidates <- interp_candidate_data(
        state$param_long, state$neighbour_results, active_core(), measurement_filter()
      )
      if (!is.null(candidates$message)) {
        log_msg(candidates$message)
        return(invisible())
      }

      report <- tryCatch(
        report_depth_coverage(
          candidates$data, sample_depths,
          target_core_id = active_core(), sample_registry = state$sample_registry
        ),
        error = function(e) {
          log_msg(paste("ERROR interpolating:", conditionMessage(e)))
          NULL
        }
      )
      if (is.null(report)) {
        return(invisible())
      }

      report$target_core_id <- active_core()
      state$depth_matched_long <- dplyr::bind_rows(
        state$depth_matched_long[state$depth_matched_long$target_core_id != active_core(), ],
        report
      )
      for (msg in depth_coverage_messages(report, active_core(), sample_depths)) log_msg(msg)
    })

    # age_model doesn't appear as its own parameter row here either -- it's
    # the age_ka column every other row already carries (see
    # depth_matched_display()).
    output$depth_matched_table <- DT::renderDataTable({
      display <- depth_matched_display(state$depth_matched_long)
      dt <- DT::datatable(display, options = list(pageLength = 15, scrollX = TRUE))
      # The empty table (before the first interpolation) has no
      # match_status column yet; DT errors if asked to style a column that
      # isn't there -- which broke this tab until something was run.
      if ("match_status" %in% names(display)) {
        dt <- DT::formatStyle(dt, "match_status",
          backgroundColor = DT::styleEqual(
            c(
              "measured", "interpolated", "no data in that depth range", "gap too large to interpolate",
              "no age model for cross-core match", "target core has no age_ka for this sample depth"
            ),
            c("white", "#fff3cd", "#f8d7da", "#f8d7da", "#f8d7da", "#f8d7da")
          )
        )
      }
      dt
    })
  })
}
