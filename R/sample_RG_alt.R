#' Sample a transitive relationship graph, alternate version
#'
#' Uses the techniques in \code{\link{enumerate_RGs_alt}} to uniformly sample
#' a relationship graph.
#'
#' @param MOIs A numeric vector specifying, for each infection, the number of
#'   distinct parasite genotypes, a.k.a. the multiplicity of infection (MOI).
#' @param igraph Logical for whether to return an \code{igraph} object.
#'
#' @return See \code{\link{enumerate_RGs_alt}}.
#'
#' @export
sample_RG_alt <- function(MOIs, igraph=T){
  # Check MOIs are positive whole numbers
  if (!all(is.wholenumber(MOIs)) | any(MOIs < 1)) stop("MOIs should be positive integers")

  infection_count <- length(MOIs) # Number of time points
  gs_count <- sum(MOIs) # Number of genotypes

  # remove restriction for now as this method should scale better
  # if (gs_count > 6 | infection_count > 3) stop("Sorry, too many RGs")
  # if (gs_count <= 1) stop("Sorry, no RGs")

  gs <- paste0("g", 1:gs_count) # Genotype names
  ts_per_gs <- rep(1:infection_count, MOIs)

  CP_list <- enumerate_CPs(MOIs)
  sizes <- sapply(CP_list, function(x) ncol(part.list[[max(x)]]))
  CP <- unlist(sample(CP_list, 1, prob=sizes))

  n.clones <- max(CP)
  clone.names <- paste0("c", 1:n.clones)
  clones <- setNames(split(gs, CP), clone.names)

  # given clonal relationships, generate all compatible sibling relationships
  if (n.clones > 10) stop("Currently supports up to 10 genotypes")

  sib.parts <- part.list[[n.clones]]
  j <- sample(1:ncol(sib.parts), 1)

  sib.vec <- sib.parts[,j]
  n.sib.clones <- max(sib.vec)
  sib.clones <- setNames(split(clone.names, sib.vec),
                         paste0("s", 1:n.sib.clones))

  RG <- list(clone=clones,
             clone.vec=CP,
             sib=sib.clones,
             sib.vec=sib.vec)

  if(igraph) RG <- RG_to_igraph(RG, gs, ts_per_gs)

  return(RG)
}
