#' Sample per-marker lineages for parasite genotypes in a relationship graph
#'
#' For each marker, the lineages of a set of first-generation offspring parasite
#' genotypes are sampled, according to their inter-genotype relationships
#' specified in a relationship graph, from a set of parents sampled from a
#' founder population, which can either be an infinitely large outbred
#' population or a finite inbred population. If its finite inbred, lineages can
#' either be equifreqent or not. Markers are treated as independent.
#'
#' \code{sample_lineages} assume inter-genotype relationships are either
#' stranger, sibling or clone and samples lineages as follows. Since clonal
#' lineages are drawn from two parental lineages, clones are mitotic copies of
#' one-another, not the products of selfing.
#' \enumerate{
#'
#' \item For a stranger, sample per-marker lineages from two parental lineages.
#'
#' \item For a clique of clones, sample per-marker lineages from two parental
#' clones for one member of the group and and clone it.
#'
#' \item For a clique of siblings, sample per-marker lineages from two
#' parental lineages.
#'
#' \item For a mixed clique of siblings and clones, for each marker, from two
#' parental lineages, sample a lineage for each clonal sub-clique and clone it.
#'
#' }
#'
#' When outbred = F (default), this function assumes parasite genotypes are drawn
#' from an inoutbredly large outbred population. As such, each parental lineage
#' is unique. When outbred = T, this function assumes parasite genotypes are
#' drawn from an inoutbredly large outbred population. As such, each parental
#' lineage is unique. More specifically, the function draws parasite genotypes
#' from an outbred population of 100 equifrequent lineages by default. As such,
#' parental genotypes are not necessarily unique. Note that, by specifying a
#' outbred population representative of a brood of parasite genotypes versus
#' parasite genotypes in the population at large, this function could be used to
#' explore IBD structure following brood versus non-brood mating.
#'
#' @param RG A relationship graph; see \code{\link{enumerate_RGs}} for details.
#' @param n_m Number of markers.
#' @param outbred Logical specifying if the parental genotypes should be sampled
#'   from a outbred population or not.
#' @param lineage_probs Vector of population lineage frequencies; if named,
#'   names will be taken as lineages; only used if \code{outbred = T}.
#' @param lineage_count Whole number specifying the number of lineages in the
#'   population to sample from; only used if \code{outbred = T} and
#'   \code{lineage_probs = NULL}.
#' @param equifrequent Logical specifying if the lineages in the population are
#'   equifrequent or not; only used if \code{outbred = T} and \code{lineage_probs
#'   = NULL}.
#' @param alpha The concentration parameter of a symmetric Dirichlet
#'   distribution over lineage frequencies when \code{equifrequent = FALSE}: set
#'   \code{0 < alpha < 1} to sample sparse frequencies, \code{alpha = 1}
#'   (default) to sample frequencies uniformly, and \code{alpha >> 1} to sample
#'   frequencies close to \code{1/lineage_count}; only used if \code{outbred = T},
#'   \code{lineage_probs = NULL} and \code{lequifrequent = F}.
#'
#' @examples
#' RG <- sample_RG(3)
#' sample_lineages(RG, n_m = 5)
#' sample_lineages(RG, n_m = 5, outbred = T)
#'
#' @export
sample_lineages <- function(RG, n_m,
                            outbred = F,
                            lineage_probs = NULL,
                            lineage_count = 10,
                            equifrequent = TRUE,
                            alpha = 1) {

  # Check n_m input
  if (!is.wholenumber(n_m) | n_m < 0) {
    stop("The marker count, n_m, must be a positive whole number")
  }

  # Extract number of genotypes
  genotype_count <- igraph::vcount(RG)

  # Extract per-component genotypes into a list
  gs_per_rcs <- igraph::groups(igraph::components(RG))

  # Generate lineages to sample
  if (outbred == F) {
    lineages = generate_lineages(length(gs_per_rcs) * 2) # Two lineages per component
  } else {
    # If lineage_probs are specified check and use them inc. names if named.
    if (!is.null(lineage_probs)) {
      if (sum(lineage_probs) != 1) stop("lineage frequencies must sum to one")
      if (any(lineage_probs < 0 | lineage_probs > 1)) stop("each lineage frequency must be between 0 and 1")
      lineage_count <- length(lineage_probs)
      if (is.null(names(lineage_probs))) names(lineage_probs) <- generate_lineages(lineage_count)
    } else { # Generate lineages and their frequencies
      if (equifrequent) {
        lineage_probs <- rep(1/lineage_count, lineage_count)
      } else {
        if (alpha <= 0) stop("alpha must be positive")
        lineage_probs <- MCMCpack::rdirichlet(1, rep(alpha, lineage_count))
        # Sort lineage probs s.t. lineages sampled are early in the alphabet (see generate_lineages)
        lineage_probs <- sort(lineage_probs, decreasing = T)
      }
      names(lineage_probs) <- generate_lineages(lineage_count)
    }
    lineages <- names(lineage_probs)
  }

  # Generate a store to keep track of lineages sampled
  lineages_sampled <- c()

  # Generate a matrix of lineages per marker (columns) for each genotype (row)
  lineages_per_marker <- array(NA, dim = c(genotype_count, n_m),
                               dimnames = list(igraph::V(RG)$name, paste0("m", 1:n_m)))

  # For each relationship component
  for(rc_i in 1:length(gs_per_rcs)) {

    # Extract component graph
    rc <- igraph::induced_subgraph(RG, gs_per_rcs[[rc_i]])
    rc_gs <- igraph::V(rc)$name
    rc_size <- igraph::vcount(rc)
    if (length(rc_gs) != rc_size) stop("Not all vertices are named")

    # Get relationship types (could be empty - strangers don't have weighted edges)
    rc_relationship_types <- unique(igraph::E(rc)$weight)

    # Sample parents
    if (outbred) {
      parents <- sample(x = lineages, size = 2, replace = T, prob = lineage_probs)
    } else {
      parents <- setdiff(lineages, lineages_sampled)[1:2]
    }

    # Update lineages sampled
    lineages_sampled <- unique(c(lineages_sampled, parents))


    if (length(rc_relationship_types) == 0) { # Stranger cluster
      # Check stranger cluster of size one
      if (rc_size != 1) stop("stranger cluster of size > 1")
      # Sample per-marker filial lineages
      lineages_per_marker[rc_gs,] <- sample(parents, n_m, replace = T)
    } else if (all(rc_relationship_types == 0.5)) { # Sibling cluster
      # Check sibling cluster of size two or more
      if (rc_size < 2) stop("sibling cluster of size < 2")
      # Sample per-marker filial lineages
      lineages_per_marker[rc_gs,] <- sample(parents, n_m * rc_size, replace = T)
    } else if (all(rc_relationship_types == 1)) { # Clonal cluster
      if (rc_size < 2) stop("clonal cluster of size < 2")
      # Sampler per-marker filial lineages (assumes clones are mitotic copies)
      filial_genotype <- sample(parents, n_m, replace = T)
      # Copy filial lineages
      lineages_per_marker[rc_gs,] <- matrix(filial_genotype, ncol = n_m, nrow = rc_size, byrow = T)
    } else if (setequal(rc_relationship_types, c(0.5, 1))) { # Mixed sibling and clonal cluster
      if (rc_size < 2) stop("mixed cluster of size < 2")
      # Duplicate relationship component
      rcd <- rc
      # Delete sibling edges
      rcd <- igraph::delete_edges(rcd, igraph::E(rc)[igraph::E(rc)$weight == 0.5])
      # Extract clonal components
      rcdC <- igraph::components(rcd)
      # For each marker
      for(m in 1:n_m) {
        # Sample one lineage per rcdC$no clonal components
        lineageC <- sample(parents, size = rcdC$no, replace = T)
        lineages_per_marker[names(rcdC$membership), m] <- lineageC[rcdC$membership] # Clone lineages and store
      }
    }
  }

  return(lineages_per_marker)
}






