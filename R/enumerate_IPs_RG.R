#' Enumerate IBD partitions consistent with a given relationship graph
#'
#' Enumerate IBD partitions consistent with a given relationship graph; see
#' vignette("enumerate", package = "Pv3Rs") for more details.
#'
#' @details
#' For a IBD partition to be consistent with the relationship graph given, it
#' must satisfy the following:
#' \itemize{
#'   \item{genotypes within each IBD cell have clonal and sibling relationships
#'   only,}
#'   \item{genotypes that are clones must be in the same IBD cell,}
#'   \item{each cluster of sibling units span at most 2 IBD cells (corresponds
#'   to two parents).}
#' }
#'
#' Note that all IBD partitions are equally likely.
#'
#' @param RG A relationship graph; see \code{\link{enumerate_RGs}} for details.
#' @param compat Logical, if true, a list of \code{partitions} equivalence
#' objects are returned. Otherwise, returns a list of IBD parititon vectors in
#' reference to the sibling clusters.
#'
#' @examples
#' set.seed(3)
#' RG <- sample_RG(c(2, 1, 2))
#' enumerate_IPs_RG(RG)
#'
#' @export
enumerate_IPs_RG <- function(RG, compat = TRUE) {
  # for each sibling cluster, `split_two` lists all ways for the cluster to
  # be split into 1 or 2 IBD cells
  ibd_per_sib <- lapply(RG$sib, split_two)
  # take cartesian product over sibling clusters
  ibd_list <- apply(
    expand.grid(ibd_per_sib, stringsAsFactors = F), 1,
    function(x) unname(purrr::flatten(x))
  )

  if (compat) {
    # return in terms of genotypes instead
    RG.clone <- RG$clone
    return(lapply(ibd_list, function(p) { # for each IBD cell
      # `RG.clone[x]`` lists the genotype names for sibling unit `x`
      out <- lapply(p, function(x) unname(unlist(RG.clone[x])))
      class(out) <- c("list", "equivalence")
      out
    }))
  } else {
    return(ibd_list)
  }
}
