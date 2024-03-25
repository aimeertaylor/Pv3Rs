#' Enumerate transitive relationship graphs, alternate version
#'
#' A relationship graph is a complete graph on all genotypes, where each edge
#' is annotated as a clone, sibling, or stranger edge. The enumerated
#' relationship graphs satisfy the following constraints:
#' \itemize{
#'   \item{The subgraph induced by the clone edges is a cluster graph.}
#'   \item{The subgraph induced by the clone edges and sibling edges is a cluster
#'   graph.}
#'   \item{Clone edges are only allowed for two genotypes from different
#'   infections.}
#' }
#'
#' This alternate version of \code{\link{enumerate_RGs}} is based on generating
#' set partitions. Since the clone edges induce a cluster graph, the
#' information encoded by clonal relationships is equivalent to a partition of
#' the genotypes. Note that genotypes from the same infection cannot belong to
#' the same partition cell. Subsequent information encoded by sibling
#' relationships is equivalent to further partitioning the clonal partition.
#' There are no constraints when enumerating the sibling partitions. The data
#' structure returned encodes each graph as a nested set partition. Each
#' partition is represented in the form of a list of vectors (`clone` and
#' `sib`) and as a membership vector (`clone.vec` and `sib.vec`), where each
#' entry identifies the partition cell the corresponding index belongs to.
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
#'   \code{\link{enumerate_RGs}}) along with these four attributes. Note that
#'   the weight matrix contains information equivalent to that of the four
#'   attributes.
#'
#' @examples
#' graphs <- enumerate_RGs_alt(c(2, 1, 2), igraph=T) # 250 graphs
#'
#' @export
enumerate_RGs_alt <- function(MOIs, igraph = TRUE) {
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

  RG_i <- 0
  RGs <- list() # list of relationship graphs

  # Count number of valid graphs and create progress bar
  n.RG <- sum(sapply(CP_list, function(x) ncol(part.list[[max(x)]])))
  pbar <- txtProgressBar(min = 0, max = n.RG) # min=0 in case n.RG is 1
  writeLines(paste("number of valid graphs is", n.RG))

  for (CP in CP_list) { # for each clonal partition (membership vector)
    n.clones <- max(CP) # number of clonal cells
    clone.names <- paste0("c", 1:n.clones)
    # list of vectors of genotype names by clonal cell
    clones <- setNames(split(gs, CP), clone.names)

    # given clonal relationships, generate all compatible sibling relationships
    # sibling partitions are all set partitions of the clonal cells
    sib.parts <- part.list[[n.clones]]
    n.parts <- ncol(sib.parts) # number of possible partitions
    for (j in 1:n.parts) { # for each sibling partition
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

      RG_i <- RG_i + 1
      RGs[[RG_i]] <- RG
      setTxtProgressBar(pbar, RG_i)
    }
  }
  writeLines("")

  RGs
}

#' Enumerate partitions induced by clonal relationships
#'
#' A clonal partition is a partition of genotypes where a pair of genotypes of
#' the same partition cell have a clonal relationship. Genotypes from the same
#' infection cannot be clones. This code enumerates all clonal partitions,
#' accounting for this intra-infection restriction.
#'
#' @param MOIs A numeric vector specifying, for each infection, the number of
#'   distinct parasite genotypes, a.k.a. the multiplicity of infection (MOI).
#'
#' @return A list of all possible partitions, where each partition is encoded
#' as a membership vector, which indices (genotype names) with the same entry
#' corresponding to genotypes being int he same partition cell.
#'
#' @examples
#' enumerate_CPs(c(2, 2))
#'
#' @export
enumerate_CPs <- function(MOIs) {
  if (!all(is.wholenumber(MOIs)) | any(MOIs < 1)) stop("MOIs should be positive integers")

  infection_count <- length(MOIs) # Number of time points
  gs_count <- sum(MOIs) # Number of genotypes

  ts_per_gs <- rep(1:infection_count, MOIs)
  gstarts <- head(c(1, cumsum(MOIs) + 1), -1) # first genotype of each infection

  # find all possible clonal relationships
  # initialise with lexicographically 'smallest' possible clone partition
  CP_init <- unlist(lapply(MOIs, function(x) 1:x))
  CP_part <- CP_init
  # prefix_max stores max element of each prefix of CP_part (excl current index)
  prefix_max <- rep(0, gs_count)
  # initialise prefix_max
  for (i in 1:gs_count) {
    if (i == 1) next
    prefix_max[i] <- max(prefix_max[i - 1], CP_part[i - 1])
  }

  CP_list <- list(CP_part)
  CP_i <- 1
  while (CP_part[gs_count] != gs_count) { # last partition is 1 2 .. gs_count
    # find genotype to increment membership entry
    for (i in gs_count:1) {
      if (CP_part[i] <= prefix_max[i]) {
        # some entry before i is already equal to CP_part[i], can increment
        break
      }
    }

    start_i <- i
    infection <- ts_per_gs[i]
    # clonal cells that are disallowed due to no intra-infection clones
    avoid <- head(CP_part[gstarts[infection]:i], -1)

    # try candidate values for new value of CP_part[i]
    candidate <- CP_part[i] + 1
    first <- TRUE
    # this loop also updates CP_part[i] for i > start_i until end of infection
    while (i <= gs_count & ts_per_gs[i] == infection) {
      # while in curr infection, account for the intra-infection restriction
      while (candidate %in% avoid) {
        candidate <- candidate + 1
      }
      CP_part[i] <- candidate
      if (first) {
        # new cell for genotype start_i also disallowed for remaining genotypes
        avoid <- c(avoid, candidate)
        candidate <- 1 # reset candidate for lexicographically closest partition
        first <- FALSE
      } else {
        candidate <- candidate + 1
      }
      i <- i + 1 # next genotype
    }

    # subsequent infections
    if (i <= gs_count) {
      # complete with lexicographically 'smallest' possible clone partition
      CP_part[i:gs_count] <- CP_init[i:gs_count]
    }

    # update prefix maximums for entries where CP_part was updated
    for (j in (start_i + 1):gs_count) {
      prefix_max[j] <- max(prefix_max[j - 1], CP_part[j - 1])
    }

    CP_i <- CP_i + 1
    CP_list[[CP_i]] <- CP_part
  }
  return(CP_list)
}
