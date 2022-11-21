#' Enumerate partitions induced by clonal relationships
#'
#' @noRd
enumerate_CPs <- function(MOIs) {
  if (!all(is.wholenumber(MOIs)) | any(MOIs < 1)) stop("MOIs should be positive integers")

  infection_count <- length(MOIs) # Number of time points
  gs_count <- sum(MOIs) # Number of genotypes

  # remove restriction for now as this method should scale better
  # if (gs_count > 6 | infection_count > 3) stop("Sorry, too many RGs")
  # if (gs_count <= 1) stop("Sorry, no RGs")

  ts_per_gs <- rep(1:infection_count, MOIs)
  gstarts <- head(c(1, cumsum(MOIs) + 1), -1) # first genotype of each infection

  # find all possible clonal relationships
  # initialise with lexicographically 'smallest' possible clone partition
  CP_init <-  unlist(lapply(MOIs, function(x) 1:x))
  CP_part <- CP_init
  # prefix_max stores largest element of each prefix of CP_part (excluding current index)
  prefix_max <- rep(0, gs_count)
  for(i in 1:gs_count) {
    if(i == 1) next
    prefix_max[i] <-  max(prefix_max[i-1], CP_part[i-1])
  }
  CP_list <- list(CP_part)
  CP_i <- 1
  while(CP_part[gs_count] != gs_count) {
    # find genotype to increment
    for(i in gs_count:1) {
      if(CP_part[i] <= prefix_max[i]) break
    }
    start_i <- i
    infection <- ts_per_gs[i]
    # clonal clones to avoid such that there are no clonal relationships within an infection
    avoid <- head(CP_part[gstarts[infection]:i], -1)
    candidate <- CP_part[i] + 1
    first <- TRUE
    while(i <= gs_count & ts_per_gs[i] == infection) { # while in current infection
      while(candidate %in% avoid) {
        candidate <- candidate + 1
      }
      CP_part[i] <- candidate
      if(first) {
        avoid <- c(avoid, candidate)
        candidate <- 1
        first <- FALSE
      } else candidate <- candidate + 1
      i <- i+1
    }
    if(i <= gs_count) {
      # complete with lexicographically 'smallest' possible clone partition
      CP_part[i:gs_count] <- CP_init[i:gs_count]
    }
    for(j in (start_i+1):gs_count) {
      prefix_max[j] <-  max(prefix_max[j-1], CP_part[j-1])
    }
    CP_i <- CP_i + 1
    CP_list[[CP_i]] <- CP_part
  }
  return(CP_list)
}

#' Enumerate transitive relationship graphs, alternate version
#'
#' This alternate version of \code{\link{enumerate_RGs}} is based on generating all
#' set partitions. This is done twice in a nested fashion: once for determining
#' clonal relationships, once for determining sibling relationships. Works for
#' up to 10 genotypes. A list of all partitions for sets of sizes 1 to 10 is
#' pre-computed.
#'
#' @param MOIs A numeric vector specifying, for each infection, the number of
#'   distinct parasite genotypes, a.k.a. the multiplicity of infection (MOI).
#'
#' @param igraph Logical for whether to return \code{igraph} objects.
#'
#' @return If \code{igraph} is \code{FALSE}, returns a list of four attributes:
#'   \describe{
#'     \item{clone}{A list of groups of genotypes that make up the clonal units.}
#'     \item{clone.vec}{A numeric vector indicating the clonal membership of each genotype.}
#'     \item{sib}{A list of groups of clonal units that make up the sibling clusters.}
#'     \item{sib.vec}{A numeric vector indicating the sibling membership of each clonal unit.}
#'   }
#'   Otherwise, returns an \code{igraph} object (see \code{\link{enumerate_RGs}})
#'   along with these four attributes. Note that the weight matrix contains the
#'   same information as these four attributes.
#'
#' @export
enumerate_RGs_alt <- function(MOIs, igraph=TRUE) {

  # Check MOIs are positive whole numbers
  if (!all(is.wholenumber(MOIs)) | any(MOIs < 1)) stop("MOIs should be positive integers")

  infection_count <- length(MOIs) # Number of time points
  gs_count <- sum(MOIs) # Number of genotypes

  # compute set partitions
  part.list <- list()
  for(i in 1:gs_count) {
    part.list[[i]] <- partitions::setparts(i)
  }

  # remove restriction for now as this method should scale better
  # if (gs_count > 6 | infection_count > 3) stop("Sorry, too many RGs")
  # if (gs_count <= 1) stop("Sorry, no RGs")

  gs <- paste0("g", 1:gs_count) # Genotype names
  ts_per_gs <- rep(1:infection_count, MOIs)

  CP_list <- enumerate_CPs(MOIs)

  RG_i <- 0
  RGs <- list()

  # Count number of valid graphs and create progress bar
  n.RG <- sum(sapply(CP_list, function(x) ncol(part.list[[max(x)]])))
  pbar <- txtProgressBar(min=0, max=n.RG) # min=0 in case n.RG is 1
  writeLines(paste("\nnumber of valid graphs is", n.RG))

  for(CP in CP_list) {
    n.clones <- max(CP)
    clone.names <- paste0("c", 1:n.clones)
    clones <- setNames(split(gs, CP), clone.names)

    # given clonal relationships, generate all compatible sibling relationships
    if (n.clones > 10) stop("Currently supports up to 10 genotypes")
    sib.parts <- part.list[[n.clones]]
    n.parts <- ncol(sib.parts)
    for(j in 1:n.parts) {
      sib.vec <- sib.parts[,j]
      n.sib.clones <- max(sib.vec)
      sib.clones <- setNames(split(clone.names, sib.vec),
                             paste0("s", 1:n.sib.clones))
      RG <- list(clone=clones,
                 clone.vec=CP,
                 sib=sib.clones,
                 sib.vec=sib.vec)

      if(igraph) RG <- RG_to_igraph(RG, gs, ts_per_gs)

      RG_i <- RG_i + 1
      RGs[[RG_i]] <- RG
      setTxtProgressBar(pbar, RG_i)
    }
  }
  writeLines("")

  RGs
}
