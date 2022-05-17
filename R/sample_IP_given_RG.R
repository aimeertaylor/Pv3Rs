#' Sample an IBD partition given a relationship graph
#'
#' \code{sample_IP_given_RG} is a wrapper function that converts lineages
#' sampled by \code{\link{sample_lineages_from_infinite_pop}} (default) or by
#' \code{\link{sample_lineages_from_finite_pop}} (if \code{outbred = FALSE})
#' into an IBD partition; see enumerate_IP.
#'
#' @param RG A relationship graph; see \code{\link{enumerate_RG}}.
#' @param outbred Logical specifying whether or not the parasite genotypes belong to
#'   an outbred population or not.
#' @param ... Additional arguments to be passed to
#'   \code{\link{sample_lineages_from_finite_pop}} if \code{outbred = FALSE}.
#'
#' @examples
#' RG <- sample_RG(3)
#' sample_IP_given_RG(RG)
#'
#' @export
sample_IP_given_RG <- function(RG, outbred = TRUE, ...) {

  if (outbred) {
    lineages <- sample_lineages_from_infinite_pop(RG)
  } else {
    lineages <- sample_lineages_from_finite_pop(RG, ...)
  }

  IP <- lineages_to_IP(lineages)
  return(IP)
}


