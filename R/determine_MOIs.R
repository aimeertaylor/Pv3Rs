#' Determine MOIs from unphased data
#'
#' The MOIs returned correspond to the parsimonious explanation of the data,
#' i.e., the minimal MOIs required to support the observed allelic diversity.
#'
#' For a given episode, the minimal MOI required to support the observed allelic
#' diversity is equal to the maximum number of per-marker alleles observed
#' across markers. These MOIs are the default used by
#' \code{\link{compute_posterior}}.
#'
#' At present, \code{Pv3Rs} only supports prevalence data (categorical data that
#' signal the detection of alleles), not quantitative (proportional abundance)
#' data. Allele repeats at markers with observed data, and repeat \code{NA}s at
#' markers with missing data, are removed. \code{NA}s in allelic vectors that
#' also contain non-\code{NA} values are ignored.
#'
#' @param y Observed data in the form of a list of lists. The outer list is a
#'   list of episodes in chronological order. The inner list is a list of named
#'   markers per episode. For each marker, one must specify an allelic vector: a
#'   set of distinct alleles detected at that marker.
#' @param return.names Logical; if TRUE and episodes are named, episode names
#'   are returned.
#'
#' @return Returns a vector of the (minimum) multiplicity of infection (MOI) for
#'   each infection.
#'
#' @examples
#' # two infections, three markers
#' y <- list(
#'   list(m1 = c("A", "B"), m2 = c("A"), m3 = c("C")),
#'   list(m1 = c("B"), m2 = c("B", "C"), m3 = c("A", "B", "C"))
#' )
#' determine_MOIs(y) # should be c(2, 3)
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
