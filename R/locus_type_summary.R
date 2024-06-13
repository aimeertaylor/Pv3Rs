#' Summarise locus-wise data types
#'
#' @details A function to summarise the data at each locus as one of four types:
#'
#' * All match (all genotypes have the same allele).
#' * All diff. (all genotypes have a different allele).
#' * Intra-match (some intra-episode genotypes have the same allele)
#' * Inter-match (some inter-episode genotypes have the same allele).
#'
#'   The number of apparent genotypes is the sum of the per-episode MOIs. When
#'   the total number of apparent genotypes exceeds 3, "Intra-match" excludes
#'   any "Inter-match" (i.e., it is equivalent to Inter-matches only), whereas
#'   an "Inter-match does not exclude an "Intra-match
#'
#' @param y A list of lists for two episodes; see [compute_posterior()] for more
#'   details.
#' @param m A positive whole number indexing a marker.
#' @return A vector of strings summarising the data at each locus.
#'
#' @examples
#' # example code
#' y <- list(
#'   list(m1 = c("A", "B"), m2 = c("A"), m3 = c("C")),
#'   list(m1 = c("B"), m2 = c("B", "C"), m3 = c("A", "B", "C", "D"))
#' )
#'
#' locus_type_summary(y)
#' @export
locus_type_summary <- function (y) {

  total_MOI <- sum(determine_MOIs(y))

  if (length(y) != 2) stop("This function only works for data on paired episodes")

  if (total_MOI > 3) {
    warning("Intra-match excludes Inter-match; Inter-match does not exclude Intra-match")
  }

  markers <- 1:length(y[[1]])

  marker_types <- sapply(markers, function(m) {

    y0 <- y[[1]][[m]]
    y1 <- y[[2]][[m]]

    if (length(intersect(y0, y1)) == 0) {
      if (length(union(y0, y1)) == total_MOI) {
        return("All diff.")
      } else {
        return("Intra-match")
      }
    } else {
      if (setequal(y0, y1)) {
        return("All match")
      } else {
        return("Inter-match")
      }
    }
  })

  return(marker_types)
}
