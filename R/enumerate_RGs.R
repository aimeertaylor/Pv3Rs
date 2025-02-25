#' Enumerate relationship graphs
#'
#' A relationship graph is a complete graph on all genotypes (one per vertex),
#' where each edge is annotated as a clone, sibling, or stranger edge. The
#' enumerated relationship graphs satisfy the following constraints:
#' \itemize{
#'   \item{The subgraph induced by the clone edges is a cluster graph.}
#'   \item{The subgraph induced by the clone edges and sibling edges is a cluster
#'   graph.}
#'   \item{Clone edges are only allowed for two genotypes from different
#'   infections.}
#' }
#'
#' Relationship graphs are enumerated by generating nested set partitions that
#' meet certain constraints; see `vignette("enumerate", package = "Pv3Rs")` for a
#' detailed description. In summary, since the clone edges induce a cluster
#' graph, the information encoded by clonal relationships is equivalent to a
#' partition of the genotypes. Note that genotypes from the same infection
#' cannot belong to the same clonal partition cell. Subsequent information
#' encoded by sibling relationships is equivalent to further partitioning the
#' clonal partition. There are no constraints when enumerating the sibling
#' partitions. The data structure returned encodes each graph as a nested set
#' partition. Each partition is represented in the form of a list of vectors
#' (`clone` and `sib`) and as a membership vector (`clone.vec` and `sib.vec`),
#' where each entry identifies the partition cell that the corresponding index
#' belongs to.
#'
#' @param MOIs A numeric vector specifying, for each infection, the number of
#'   distinct parasite genotypes, a.k.a. the multiplicity of infection (MOI).
#' @param igraph Logical for whether to return \code{igraph} objects.
#'
#' @return A list of relationship graphs. If \code{igraph} is \code{FALSE},
#' each element is a list of four attributes:
#'   \describe{
#'     \item{clone}{A list of groups of genotypes that make up the clonal
#'     cells.}
#'     \item{clone.vec}{A numeric vector indicating the clonal membership of
#'     each genotype.}
#'     \item{sib}{A list of groups of clonal cells that make up the sibling
#'     cells.}
#'     \item{sib.vec}{A numeric vector indicating the sibling membership of
#'     each clonal cell.}
#'   }
#'   Otherwise, each element is an \code{igraph} object (see
#'   \code{\link{RG_to_igraph}}) along with these four attributes. Note that
#'   the weight matrix contains information equivalent to that of the four
#'   attributes.
#'
#' @examples
#' graphs <- enumerate_RGs(c(2, 1, 2), igraph=TRUE) # 250 graphs
#'
#' @export
enumerate_RGs <- function(MOIs, igraph = TRUE) {
  # Check MOIs are positive whole numbers
  if (!all(is.wholenumber(MOIs)) | any(MOIs < 1)) stop("MOIs should be positive integers")

  if(sum(MOIs) > 10) message("Total MOI > 10 may lead to high memory use")

  infection_count <- length(MOIs) # Number of time points
  gs_count <- sum(MOIs) # Number of genotypes

  # compute set partitions
  part.list <- list()
  for (i in 1:gs_count) {
    part.list[[i]] <- partitions::setparts(i)
  }

  gs <- paste0("g", 1:gs_count) # Genotype names
  ts_per_gs <- rep(1:infection_count, MOIs)

  # get clonal partitions, accounting for no intra-infection clones
  CP_list <- enumerate_CPs(MOIs)

  RG_i <- 0
  RGs <- list() # list of relationship graphs

  # Count number of valid graphs and create progress bar
  n.RG <- sum(sapply(CP_list, function(x) ncol(part.list[[max(x)]])))
  pbar <- msg_progress_bar(n.RG)
  message(paste("Number of valid relationship graphs (RGs) is", n.RG))

  for (CP in CP_list) { # for each clonal partition (membership vector)
    n.clones <- max(CP) # number of clonal cells
    clone.names <- paste0("c", 1:n.clones)
    # list of vectors of genotype names by clonal cell
    clones <- stats::setNames(split(gs, CP), clone.names)

    # given clonal relationships, generate all compatible sibling relationships
    # sibling partitions are all set partitions of the clonal cells
    sib.parts <- part.list[[n.clones]]
    n.parts <- ncol(sib.parts) # number of possible partitions
    for (j in 1:n.parts) { # for each sibling partition
      sib.vec <- sib.parts[, j] # membership vector
      n.sib.clones <- max(sib.vec) # number of sibling cells
      sib.clones <- stats::setNames(
        split(clone.names, sib.vec),
        paste0("s", 1:n.sib.clones)
      ) # list of vectors of clonal cell names by sibling cell
      RG <- list(
        clone = clones,
        clone.vec = CP,
        sib = sib.clones,
        sib.vec = sib.vec
      )

      if (igraph) RG <- RG_to_igraph(RG, gs, ts_per_gs)

      RG_i <- RG_i + 1
      RGs[[RG_i]] <- RG
      pbar$increment()
    }
  }

  RGs
}

