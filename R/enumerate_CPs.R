#' Enumerate clonal partitions
#'
#' A clonal partition is a partition of genotypes where a pair of genotypes of
#' the same partition cell have a clonal relationship. Genotypes from the same
#' infection cannot be clones. This code enumerates all clonal partitions,
#' accounting for this intra-infection restriction; see
#' [vignette("enumerate")](library/Pv3Rs/doc/enumerate.pdf) for more details.
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
  gstarts <- utils::head(c(1, cumsum(MOIs) + 1), -1) # first genotype of each infection

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
    avoid <- utils::head(CP_part[gstarts[infection]:i], -1)

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
