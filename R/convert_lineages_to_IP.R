#' Converts a vector of lineages to an IBD partition.
#'
#' \code{convert_lineages_to_IP} is copied from \code{\link[partitions]{vec_to_eq}} and
#' modified so that the output is comparable with the output of
#' \code{\link[partitions]{listParts}}. Otherwise stated, the output is ordered so that large
#' clusters precede small clusters and so that low genotype indices precede
#' high genotype indices.
#'
#' @param lineages_m For a single marker, a character vector lineages whose names are parasite genotype
#'   indices.
#'
#' @examples
#' set.seed(1)
#' RG <- sample_RG(3)
#' plot_RG(RG)
#' lineages <- sample_lineages(RG, n_m = 3)
#' apply(lineages, 2, convert_lineages_to_IP)
#'
#' @export
convert_lineages_to_IP <- function(lineages_m) {
  out <- split(names(lineages_m), lineages_m, lex.order = T)
  out <- out[names(sort(sapply(out, min)))] # Order s.t. low genotype indices precede high
  len <- sapply(out, length)
  out <- out[order(len, decreasing = T)] # Order s.t. large clusters precede small clusters
  class(out) <- c(class(out), "equivalence")
  return(out)
}
