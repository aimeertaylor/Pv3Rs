#' Sample IBD partition/s given a relationship graph
#'
#' \code{sample_IPs_given_RG} is a wrapper function that converts per-marker lineages
#' sampled by \code{\link{sample_lineages}} into per-marker IBD partitions; see enumerate_IP.
#'
#' @param RG A relationship graph; see \code{\link{enumerate_RGs_prev}}.
#' @param n_m Number of markers.
#' @param outbred Logical specifying whether or not the parasite genotypes belong to
#'   an outbred population or not.
#' @param ... Additional arguments to be passed to
#'   \code{\link{sample_lineages}} when \code{outbred = FALSE}.
#'
#' @examples
#' set.seed(1)
#' RG <- sample_RG_prev(3)
#' sample_IPs_given_RG(RG, n_m = 3)
#'
#' @export
sample_IPs_given_RG <- function(RG, n_m, outbred = TRUE, ...) {
  lineages <- sample_lineages(RG, n_m, outbred, ...)
  IPs <- apply(lineages, 2, convert_lineages_to_IP)
  return(IPs)
}
