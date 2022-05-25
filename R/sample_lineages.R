#' Sample parasite lineages given a relationship graph
#'
#' The lineages of parasite genotypes are sampled from some population according to
#' their inter-genotype relationships specified in a relationship
#' graph.
#'
#' Functions \code{sample_lineages_from_infinite_pop} and
#' \code{sample_lineages_from_finite_pop} assume inter-genotype relationships
#' are either stranger, sibling or clone and sample lineages as follows.
#' \enumerate{
#'
#' \item For a stranger, sample a parental lineage.
#'
#' \item For all members of each clonal clique of the relationship graph, sample
#' one parental lineage and clone it.
#'
#' \item For a sibling clique of size k, sample two parental lineages and then k
#' filial lineages from the parental lineages.
#'
#' \item For a mixed sibling and clonal clique, sample two parental lineages and
#' then for all members of each clonal sub-clique, sample one filial lineage
#' from the two parental lineages and clone it.
#'
#' }
#'
#' @param RG A relationship graph; see \code{\link{enumerate_RGs}} for details.
#'
#' @section To-do:
#' To speed up, could group all rc stranger clusters of size one
#'
#' @name sample_lineages
NULL
#> NULL




#' @describeIn sample_lineages Sample from a infinite population
#'
#'   This function assumes parasite genotypes are drawn from an infinitely large
#'   outbred population. As such, each parental lineage is unique.
#'
#' @examples
#' RG <- sample_RG(3)
#' sample_lineages_from_infinite_pop(RG)
#'
#' @export
sample_lineages_from_infinite_pop <- function(RG) {

  # Extract number of genotypes
  genotype_count <- igraph::vcount(RG)

  # Generate lineages with one extra to compensate for dummy's lineage = "a"
  lineages = generate_lineages(genotype_count + 1)

  # Initiate lineage store
  RG_lineages <- c(dummy = "a", array(NA, dim = genotype_count, dimnames = list(igraph::V(RG)$name)))
  gs_per_rcs <- igraph::groups(igraph::components(RG)) # group extracts per-component genotypes into a list

  # For each relationship component
  for(i in 1:length(gs_per_rcs)) {

    rc <- igraph::induced_subgraph(RG, gs_per_rcs[[i]]) # extract component graph
    rc_gs <- igraph::V(rc)$name
    rc_size <- igraph::vcount(rc)
    if (length(rc_gs) != rc_size) stop("Not all vertices are named")
    rc_relationship_types <- unique(igraph::E(rc)$weight) # Get relationship types (could be empty - strangers are empty)

    if (length(rc_relationship_types) == 0) { # Stranger cluster
      # Check stranger cluster of size one
      if (rc_size != 1) stop("stranger cluster of size > 1")
      # Sample one lineage and store
      RG_lineages[rc_gs] <- setdiff(lineages, RG_lineages)[1] # Take a new lineage
    } else if (all(rc_relationship_types == 0.5)) { # Sibling cluster
      # Check sibling cluster of size two or more
      if (rc_size < 2) stop("sibling cluster of size < 2")
      # Sample parent lineages
      parents <- setdiff(lineages, RG_lineages)[1:2] # Take two new lineages
      RG_lineages[rc_gs] <- sample(parents, rc_size, replace = T) # Sample children lineages and store
    } else if (all(rc_relationship_types == 1)) { # Clonal cluster
      # Take one new lineage, clone it and store it
      RG_lineages[rc_gs] <- rep(setdiff(lineages, RG_lineages)[1], rc_size)
    } else if (setequal(rc_relationship_types, c(0.5, 1))) { # Mixed sibling and clonal cluster
      # Sample parent lineages
      parents <- setdiff(lineages, RG_lineages)[1:2] # Take two new lineages
      rcd <- rc # Duplicate relationship component
      rcd <- igraph::delete_edges(rcd, igraph::E(rc)[igraph::E(rc)$weight == 0.5]) # Delete sibling edges
      rcdC <- igraph::components(rcd) # Extract clonal components within mixed sibling-clone component
      lineageC <- sample(parents, size = rcdC$no, replace = T) # Sample one lineage per rcdC$no clonal components
      RG_lineages[names(rcdC$membership)] <- lineageC[rcdC$membership] # Clone lineages and store
    }
  }
  return(RG_lineages[-1])
}


#' @describeIn sample_lineages Sample from a finite population
#'
#'   This function assumes parasite genotypes are drawn from an finite
#'   population of 100 equifrequent lineages by default. As such, parental
#'   genotypes are not necessarily unique. Note that, by specifying a finite
#'   population representative of a brood of parasite genotypes versus parasite
#'   genotypes in the population at large, this function could be used to explore
#'   IBD structure following brood versus non-brood mating.
#'
#' @param lineage_probs List of population lineages and their frequencies
#' @param lineage_count Whole number specifying the number of lineages in the
#'   population to sample from.
#' @param equifrequent Logical specifying if the lineages in the population are
#'   equifrequent or not.
#' @param alpha The concentration parameter of a symmetric Dirichlet
#'   distribution over lineage frequencies when \code{equifrequent = FALSE}: set
#'   \code{0 < alpha < 1} to sample sparse frequencies, \code{alpha = 1} to sample
#'   frequncies uniformly, and \code{alpha >> 1} to sample frequencies close to
#'   \code{1/lineage_count}.
#' @param probs A probability vector specifying the lineage frequencies used
#'   only if \code{alpha = NULL}.
#'
#' @section To-do:
#' either return probs (in a way that doesn't mess up sample_IP_given_RG)
#'
#' @examples
#' RG <- sample_RG(3)
#' sample_lineages_from_finite_pop(RG, lineage_count = 5, equifrequent = FALSE, alpha = 0.1)
#'
#' @export
sample_lineages_from_finite_pop <- function(RG, lineage_probs = NULL, lineage_count = 100,
                                            equifrequent = TRUE, alpha = 1,
                                            probs = NULL) {

  if (is.null(lineage_probs)) { # Generate lineages and their frequencies
    if(equifrequent) {
      lineage_probs <- rep(1/lineage_count, lineage_count)
    } else {
      if (!is.NULL(alpha)) {
        if (alpha <= 0) stop("alpha must be positive")
        lineage_probs <- MCMCpack::rdirichlet(1, rep(alpha, lineage_count))
      } else {
        if (sum(probs) != 1 | length(probs) != lineage_count) {
          stop("probs needs to sum to one and have length equal to the lineage_count")
        }
      }
    }
     names(lineage_probs) <- generate_lineages(lineage_count)
  }

  lineages <- names(lineage_probs)
  genotype_count <- igraph::vcount(RG)
  RG_lineages <- array(NA, dim = genotype_count, dimnames = list(igraph::V(RG)$name)) # Initiate lineage store
  gs_per_rcs <- igraph::groups(igraph::components(RG)) # group extracts per-component genotypes into a list

  for(rc_i in 1:length(gs_per_rcs)) { # For each relationship component

    rc <- igraph::induced_subgraph(RG, gs_per_rcs[[rc_i]]) # extract component graph
    rc_gs <- igraph::V(rc)$name
    rc_size <- igraph::vcount(rc)
    if (length(rc_gs) != rc_size) stop("Not all vertices are named")
    rc_lineages <- array(NA, dim = rc_size, dimnames = list(rc_gs)) # Initiate lineage store
    rc_relationship_types <- unique(igraph::E(rc)$weight) # Get relationship types (could be empty - strangers don't have a "w" attribute)

    if (length(rc_relationship_types) == 0) { # Stranger cluster
      # Check stranger cluster of size one
      if (rc_size != 1) stop("stranger cluster of size > 1")
      # Sample one lineage and store
      RG_lineages[rc_gs] <- sample(lineages, size = 1, prob = lineage_probs)
    } else if (all(rc_relationship_types == 0.5)) { # Sibling cluster
      # Check sibling cluster of size two or more
      if (rc_size < 2) stop("sibling cluster of size < 2")
      # Sample parent lineages
      parents <- sample(lineages, size = 2, replace = T, prob = lineage_probs)
      # Sample children lineages and store
      RG_lineages[rc_gs] <- sample(parents, rc_size, replace = T)
    } else if (all(rc_relationship_types == 1)) { # Clonal cluster
      # Sample one lineage, clone and store it
      RG_lineages[rc_gs] <- rep(sample(lineages, 1, prob = lineage_probs), rc_size)
    } else if (setequal(rc_relationship_types, c(0.5, 1))) { # Mixed sibling and clonal cluster
      # Sample parent lineages
      parents <- sample(lineages, size = 2, replace = T, prob = lineage_probs)
      rcd <- rc # Duplicate relationship component
      rcd <- igraph::delete_edges(rcd, igraph::E(rc)[igraph::E(rc)$weight == 0.5]) # Delete sibling edges
      rcdC <- igraph::components(rcd) # Extract clonal components
      lineageC <- sample(parents, size = rcdC$no, replace = T) # Sample one lineage per rcdC$no clonal components
      RG_lineages[names(rcdC$membership)] <- lineageC[rcdC$membership] # Clone lineages and store
    }
  }
  return(RG_lineages)
}

#' @describeIn sample_lineages Generates n population lineages
#'
#'   Letters represent lineages with repeats if more than 26 lineages are
#'   required.
#'
#' @param n Number of lineages to generate.
#'
#' @examples
#' rev(generate_lineages(28))
#'
#' @export
generate_lineages <- function(n) {
  lineages <- letters
  while (length(lineages) < n) {
    reps <- ceiling(n/length(lineages))
    new_lineages <- apply(matrix(rep(letters, each = reps), nrow = reps), 2, paste0, collapse = "")
    lineages <- c(lineages, new_lineages)
  }
  return(lineages[1:n])
}


