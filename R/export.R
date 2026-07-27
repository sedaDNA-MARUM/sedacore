#' Write a harmonized `param_long` table to TSV
#'
#' @param x A `param_long` tibble.
#' @param path Output file path.
#' @return `path`, invisibly.
#' @export
write_harmonized_tsv <- function(x, path) {
  validate_param_long(x)
  readr::write_tsv(x, path)
  invisible(path)
}
