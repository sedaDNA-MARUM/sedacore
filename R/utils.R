#' Guard against locale-dependent corruption of PANGAEA text
#'
#' PANGAEA column names, units, and species labels routinely contain
#' non-ASCII characters (delta, permil, degree, mu signs). Empirically
#' confirmed: when the R session's `LC_CTYPE` is not UTF-8-capable (e.g.
#' the bare `"C"` locale), those characters are silently and irreversibly
#' replaced with literal `"<U+03B4>"`-style escape text somewhere in the
#' download/parse path -- and `pangaear::pg_data()` caches that
#' already-corrupted result to disk, so even fixing the locale afterwards
#' won't help until the cache entry is cleared. Rather than let this
#' silently corrupt species names and units downstream, fail loudly before
#' doing anything.
#'
#' @param startup Set `TRUE` when called from `.onLoad()`: uses
#'   `packageStartupMessage()` (suppressible via
#'   `suppressPackageStartupMessages()`, required for `R CMD check`
#'   portability) instead of `rlang::inform()`.
#' @return Invisibly `TRUE`. Aborts only if no UTF-8 locale is available on
#'   the system to switch to.
#' @noRd
check_utf8_locale <- function(startup = FALSE) {
  if (isTRUE(l10n_info()[["UTF-8"]])) {
    return(invisible(TRUE))
  }

  original <- Sys.getlocale("LC_CTYPE")
  candidates <- c("en_US.UTF-8", "C.UTF-8", "en_GB.UTF-8", "UTF-8")
  for (loc in candidates) {
    switched <- suppressWarnings(tryCatch(Sys.setlocale("LC_CTYPE", loc), error = function(e) ""))
    if (!identical(switched, "") && isTRUE(l10n_info()[["UTF-8"]])) {
      msg <- sprintf(
        "LC_CTYPE was \"%s\" (not UTF-8); switched to \"%s\" for this session to avoid corrupting PANGAEA text.",
        original, loc
      )
      if (startup) {
        packageStartupMessage(msg)
      } else {
        rlang::inform(msg, class = "sedacore_locale_autoswitch")
      }
      return(invisible(TRUE))
    }
  }

  rlang::abort(
    c(
      "R session locale is not UTF-8-capable; refusing to fetch PANGAEA data.",
      "x" = sprintf("Current LC_CTYPE: \"%s\"; no UTF-8 locale found to switch to.", original),
      "i" = paste(
        "PANGAEA column names/units use non-ASCII characters",
        "(delta, permil, degree signs) that get silently corrupted into",
        "literal \"<U+xxxx>\" text under a non-UTF-8 locale, and that",
        "corruption gets cached to disk by pangaear."
      ),
      "i" = "Fix: install a UTF-8 locale on this system, e.g. en_US.UTF-8."
    ),
    class = "sedacore_locale_error"
  )
}

#' Detect PANGAEA text corrupted by a non-UTF-8 locale
#'
#' Defense-in-depth backstop for `check_utf8_locale()`: even with a UTF-8
#' locale confirmed up front, abort loudly if any already-corrupted text
#' shows up (e.g. from a stale `pangaear` disk cache written under a
#' previous non-UTF-8 session) instead of silently returning wrong data.
#'
#' @param x A character vector.
#' @return `TRUE` if any element looks like `"<U+03B4>"`-style corruption.
#' @noRd
has_mojibake <- function(x) {
  any(grepl("<U\\+[0-9A-Fa-f]{4,6}>", x))
}

#' @noRd
abort_if_mojibake <- function(x, doi) {
  if (has_mojibake(x)) {
    rlang::abort(
      c(
        sprintf("Corrupted (mojibake) text detected in PANGAEA dataset %s.", doi),
        "i" = paste(
          "This usually means the dataset was fetched and cached under a",
          "non-UTF-8 locale in a previous session."
        ),
        "i" = "Fix: clear the pangaear cache (pangaear::pg_cache$delete_all()) and re-fetch."
      ),
      class = "sedacore_mojibake_error"
    )
  }
  invisible(TRUE)
}
