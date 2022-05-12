#' Function to test if x contains integer numbers.
#'
#' This function is copied from \code{\link[base]{integer}}. It is differs from \code{\link[base]{is.integer}}
#' because tests if x is an integer type, not an integer number.
#'
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol}
