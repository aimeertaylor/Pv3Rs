#' Sample a transitive relationship graph, alternate version
#'
#' Uses the techniques in \code{\link{enumerate_RGs}} to uniformly sample
#' a relationship graph. All clonal partitions are generated, and the number
#' of sibling partitions consistent with each clonal partition is determined. A
#' clonal partition is randomly selected with probability proportional to the
#' corresponding number of sibling partitions, and a sibling partition is then
#' uniformly sampled. The nested partition is equivalent to a relationship
#' graph. See \code{\link{enumerate_RGs}} for details on the nested
#' partition representation of a relationship graph.
#'
#' @param MOIs A numeric vector specifying, for each infection, the number of
#'   distinct parasite genotypes, a.k.a. the multiplicity of infection (MOI).
#' @param igraph Logical for whether to return an \code{igraph} object.
#'
#' @return A relationship graph, i.e. one entry of the list returned by
#'   \code{\link{enumerate_RGs}}.
#'
#' @examples
#' set.seed(20)
#' RG <- sample_RG(c(2, 2))
#'
#' @export
sample_RG <- function(MOIs, igraph = T) {
  # Check MOIs are positive whole numbers
  if (!all(is.wholenumber(MOIs)) | any(MOIs < 1)) stop("MOIs should be positive integers")

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
  ts_per_gs <- rep(1:infection_count, MOIs)

  # get clonal partitions, accounting for no intra-infection clones
  CP_list <- enumerate_CPs(MOIs)
  # number of sibling partitions for each clonal partition
  sizes <- sapply(CP_list, function(x) ncol(part.list[[max(x)]]))
  # random selection of clonal partition with prob prop to # sibling partitions
  CP <- unlist(sample(CP_list, 1, prob = sizes)) # membership vector

  n.clones <- max(CP) # number of clonal cells
  clone.names <- paste0("c", 1:n.clones)
  # list of vectors of genotype names by clonal cell
  clones <- setNames(split(gs, CP), clone.names)

  # given clonal relationships, generate all compatible sibling relationships
  # sibling partitions are all set partitions of the clonal cells
  sib.parts <- part.list[[n.clones]]
  j <- sample(1:ncol(sib.parts), 1) # uniform selection of sibling partition

  sib.vec <- sib.parts[, j] # membership vector
  n.sib.clones <- max(sib.vec) # number of sibling cells
  sib.clones <- setNames(
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

  return(RG)
}
