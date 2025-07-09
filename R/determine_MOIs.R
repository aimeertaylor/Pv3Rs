#' Determine multiplicities of infection (MOIs)
#'
#' Returns a MOI estimate for each episode based on allelic diversity across
#' markers.
#'
#' A true MOI is a number of genetically distinct groups of clonal parasites
#' within an infection. Give or take *de novo* mutations, all parasites within
#' a clonal group share the same DNA sequence, which we call a genotype. As
#' such, MOIs are distinct parasite genotype counts. Under the Pv3Rs model
#' assumption that there are no genotyping errors, the true MOI of an episode is
#' greater than or equal to the maximum distinct allele count for any marker in
#' the data on that episode. In other words, under the assumption of no
#' genotyping errors, maximum distinct allelic counts are the most
#' parsimonious MOI estimates compatible with the data. By default, these MOI
#' estimates are used by \code{\link{compute_posterior}}.
#'
#' @param y List of lists encoding allelic data; see
#'   \code{\link{compute_posterior}} for more details. The outer list contains
#'   episodes in chronological order. The inner list contains named markers per
#'   episode. For each marker, one must specify an allelic vector: a set of
#'   distinct alleles detected at that marker; or `NA` if marker data are
#'   missing.
#' @param return.names Logical; if TRUE and `y` has named episodes, episode
#'   names are returned.
#'
#' @return Numeric vector containing one MOI estimate per episode, each estimate
#'   representing the maximum number of distinct alleles observed at any marker per
#'   episode.
#'
#' @examples
#' y <- list(enrol = list(m1 = c("A", "B"), m2 = c("A"), m3 = c("C")),
#'           recur = list(m1 = c("B"), m2 = c("B", "C"), m3 = c("A", "B", "C")))
#' determine_MOIs(y) # returns c(2, 3)
#'
#' @export
determine_MOIs <- function(y, return.names = FALSE) {
  if(return.names) {
    if(is.null(names(y))) warning("unnamed episodes")
    sapply(prep_data(y), function(x) max(sapply(x, length)))
  } else {
    unname(sapply(prep_data(y), function(x) max(sapply(x, length))))
  }
}
