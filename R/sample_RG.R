#' Sample a relationship graph (RG)
#'
#' Uses the techniques in \code{\link{enumerate_RGs}} to sample a single RG
#' uniformly. All clonal partitions are generated, each weighted by its number
#' of consistent sibling partitions. A clonal partition is sampled proportional
#' to its weight, then a consistent sibling partition is drawn uniformly. The
#' resulting nested partition represents the RG; see \code{\link{enumerate_RGs}}
#' for details.
#'
#' @param MOIs Vector of per-episode multiplicities of infection (MOIs),
#' i.e., numbers of vertices per episode.
#' @param igraph Logical; if `TRUE` (default), returns the RG as an
#'   \code{igraph} object.
#'
#' @return An RG encoded either as an `igraph` object (default), or as a list;
#'   see \code{\link{enumerate_RGs}} for details.
#'
#' @examples
#' set.seed(1)
#' RG <- sample_RG(c(3, 2))
#' plot_RG(RG)
#'
#' @export
sample_RG <- function(MOIs, igraph = TRUE) {
  # Check MOIs are positive whole numbers
  if (!all(is_wholenumber(MOIs)) | any(MOIs < 1)) stop("MOIs should be positive integers")

  if(sum(MOIs) > 10) warning(
    "Total MOI > 10 may lead to high memory use", immediate=T
  )

  infection_count <- length(MOIs) # Number of time points
  gs_count <- sum(MOIs) # Number of genotypes

  # compute set partitions
  part.list <- list()
  for (i in 1:gs_count) {
    part.list[[i]] <- partitions::setparts(i)
  }

  gs <- paste0("g", 1:gs_count) # Genotype names

  # get clonal partitions, accounting for no intra-infection clones
  CP_list <- enumerate_CPs(MOIs)
  # number of sibling partitions for each clonal partition
  sizes <- sapply(CP_list, function(x) ncol(part.list[[max(x)]]))
  # random selection of clonal partition with prob prop to # sibling partitions
  CP <- unlist(sample(CP_list, 1, prob = sizes)) # membership vector

  n.clones <- max(CP) # number of clonal cells
  clone.names <- paste0("c", 1:n.clones)
  # list of vectors of genotype names by clonal cell
  clones <- stats::setNames(split(gs, CP), clone.names)

  # given clonal relationships, generate all compatible sibling relationships
  # sibling partitions are all set partitions of the clonal cells
  sib.parts <- part.list[[n.clones]]
  j <- sample(1:ncol(sib.parts), 1) # uniform selection of sibling partition

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

  if (igraph) RG <- RG_to_igraph(RG, MOIs)

  return(RG)
}
