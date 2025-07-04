#' Determine multiplicities of infection (MOIs)
#'
#' Calculates a MOI for each episode based on allelic diversity
#' across markers.
#'
#' For a given episode, the MOI is the maximum number distinct alleles observed
#' at any marker within that episode. These MOIs are the default used by
#' \code{\link{compute_posterior}}.
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
#' @return Numeric vector containing one MOI per episode, each MOI representing
#'   the maximum number of alleles observed at any marker per episode.
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
