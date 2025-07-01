#' Test if x is an integer number.
#'
#' Internal function copied from [base::integer()]. It is differs from
#' [base::is.integer()] because [base::is.integer()] tests if x is an integer
#' type (for compatibility with C and Fortran code).
#'
#' @noRd
is_wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}
