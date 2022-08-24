#' Enumerate IBD partitions given RG
#'
#' The rules are (a) each IBD unit has clonal and sibling relationships only,
#' (b) clones are IBD, (c) each sibling cluster is part of at most 2 IBD units.
#'
#' Note that all IBD partitions are equally likely.
#'
#' @param RG A relationship graph; see \code{\link{enumerate_RGs}} for details.
#' @param compat Logical, if true, the data type of the output matches with the
#' input IP of \code{\link{compute_pr_IP_RG}}
#'
#' @export
enumerate_IPs_RG <- function(RG, compat=TRUE) {
  gs <- igraph::V(RG)$name
  gs_count <- length(gs)

  # IBD partition over each sibling cluster
  ibd_per_sib <- lapply(RG$sib, split_two)
  # take cartesian product over sibling clusters
  ibd_list <- apply(expand.grid(ibd_per_sib), 1, function(x) unname(flatten(x)))

  if(compat) {
    # return in terms of genotypes instead
    RG.clone <- RG$clone
    return(lapply(ibd_list, function(p) {
      out <- lapply(p, function(x) unname(unlist(RG.clone[x])))
      class(out) <- c("list", "equivalence")
      out
    }))
  } else return(ibd_list)
}

