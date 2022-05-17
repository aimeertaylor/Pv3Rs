#' Converts a vector of per-genotype lineages to an IBD partition.
#'
#' \code{lineage_to_IP} is copied from \code{\link[partitions]{vec_to_eq}} and
#' modified so that the output is comparable with the output of
#' \code{\link[partitions]{listParts}}. Otherwise stated, the output is ordered so that large
#' clusters precede small clusters and so that low genotype indices precede
#' high genotype indices.
#'
#' @param lin_vec A character vector of parasites lineages whose names are genotype
#'   indices.
#'
#' @examples
#' RG <- sample_RG(3)
#' lineages <- sample_lineages_from_infinite_pop(RG)
#' lineages_to_IP(lineages)
#'
#' @export
lineages_to_IP <- function(lin_vec) {
  out <- split(names(lin_vec), lin_vec, lex.order = T)
  out <- out[names(sort(sapply(out, min)))] # Order s.t. low genotype indices precede high
  len <- sapply(out, length)
  out <- out[order(len, decreasing = T)] # Order s.t. large clusters precede small clusters
  class(out) <- c(class(out), "equivalence")
  return(out)
}



