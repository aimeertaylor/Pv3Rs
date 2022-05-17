#' Test if x contains integer numbers.
#'
#' Internal function copied from \code{\link[base]{integer}}. It is differs from
#' \code{\link[base]{is.integer}} because tests if x is an integer type, not an
#' integer number.
#'
#' @noRd
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol}
