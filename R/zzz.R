#' @noRd
.onLoad <- function(libname, pkgname) {
  # sedacore deals in scientific text (delta-18-O, permil, degree signs)
  # throughout, not just PANGAEA fetches -- so this needs to run at load
  # time, not just before network calls. See check_utf8_locale().
  check_utf8_locale(startup = TRUE)
}
