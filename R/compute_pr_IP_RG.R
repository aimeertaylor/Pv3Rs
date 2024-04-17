#' Compute the probability of an IBD partition given relationship graph
#'
#' \code{compute_pr_IP_RG} assumes strangers, siblings, and clones are
#' represented numerically as 0, 0.5 and 1, respectively, and that the
#' relatedness (probability of IBD) of a given relationship is equal to its
#' numeric representation, implying that parasite genotypes are sampled from an
#' outbred population.
#'
#' Superseded by \code{\link{enumerate_IPs_RG}} as the probability distribution
#' of IBD partitions is in fact uniform.
#'
#' @param IP An IBD partition (IP), which is a \code{partitions} equivalence object; see \code{\link[partitions]{listParts}}.
#' @param RG A relationship graph (RG), which is an \code{igraph} graph; see \code{\link[igraph]{igraph}}.
#'
#' @return A probability from 0 to 1
#'
#' @examples
#'
#' #################################################
#' # Simple example of computation
#' #################################################
#'
#' n_genotypes <- 3
#' RG <- sample_RG_prev(n_genotypes)
#' IPs <- partitions::listParts(n_genotypes)
#' sapply(IPs, compute_pr_IP_RG, RG = RG)
#' plot_RG(RG)
#'
#'
#' #################################################
#' # Checking the computation against simulation
#' #################################################
#'
#' # Generate a relationship graph (RG) and enumerate IBD partitions (IPs)
#' set.seed(3)
#' RG <- sample_RG_prev(MOIs = c(2, 1, 1))
#' IPs <- enumerate_IPs(igraph::vcount(RG))
#'
#' # Compute the probability of the IPs given the RG
#' pr_IP_RG <- sapply(IPs, compute_pr_IP_RG, RG)
#' names(pr_IP_RG) <- sapply(IPs, convert_IP_to_string)
#'
#' # Make some stores for simulation results
#' sim_out_fr <- sim_in_fr <- pr_IP_RG
#' sim_out_fr[] <- sim_in_fr[] <- 0
#'
#' # Compute the frequency of simulated IPs given the RG graph and
#' # outbred population (well specified setting)
#' simulated_out <- sapply(1:1000, function(j) {
#'   convert_IP_to_string(sample_IPs_given_RG(RG, n_m = 1, outbred = T)[[1]])
#' })
#' simulated_out_fr <- table(simulated_out) / 1000
#' sim_out_fr[names(simulated_out_fr)] <- simulated_out_fr
#'
#' # Compute the frequency of simulated IPs given the RG graph and
#' # inbred population (miss specified setting))
#' simulated_in <- sapply(1:1000, function(j) {
#'   convert_IP_to_string(sample_IPs_given_RG(RG, n_m = 1, outbred = F, lineage_count = 10)[[1]])
#' })
#' simulated_in_fr <- table(simulated_in) / 1000
#' sim_in_fr[names(simulated_in_fr)] <- simulated_in_fr
#'
#' # Plot agreement and relationship graph
#' par(mfrow = c(1, 2), pty = "s")
#' plot_RG(RG)
#' max_xy <- max(c(pr_IP_RG, sim_out_fr, sim_in_fr))
#' plot(
#'   x = pr_IP_RG,
#'   y = as.vector(sim_out_fr[names(pr_IP_RG)]),
#'   xlab = "Probability", ylab = "Frequency",
#'   xlim = c(0, max_xy + 0.1),
#'   ylim = c(0, max_xy + 0.1),
#'   pch = 4
#' )
#' points(
#'   x = pr_IP_RG,
#'   y = as.vector(sim_in_fr[names(pr_IP_RG)]), col = "red"
#' )
#' abline(a = 0, b = 1)
#' legend("topleft",
#'   pch = c(1, 4), col = c(1:2),
#'   legend = c("well-specified", "miss-specified"),
#'   cex = 0.75, inset = 0.01
#' )
#'
#' @export
compute_pr_IP_RG <- function(IP, RG) {
  # Convert RG into a sparse matrix
  weight <- if (igraph::is.weighted(RG)) "weight" else NULL
  RG_matrix <- igraph::as_adjacency_matrix(RG, attr = weight)

  # Extract number of IBD clusters
  cluster_count <- length(IP)

  # Compute inter-cluster non-join prob
  if (cluster_count > 1) {
    # Create a matrix to store inter-cluster closest relatives
    inter_cluster_matrix <- array(0, dim = rep(cluster_count, 2))

    # Extract the relationships between IBD clusters
    for (i in 2:cluster_count) {
      for (j in 1:(i - 1)) {
        inter_cluster_matrix[i, j] <- extract_closest_relative(IP[[i]], IP[[j]], RG_matrix)
      }
    }

    # Compute inter-cluster non-join prob
    inter <- compute_inter_cluster_join(inter_cluster_matrix)
  } else {
    inter <- 1
  }

  # Compute intra-cluster join prob
  intra <- sapply(IP, compute_intra_cluster_join, RG_matrix)

  # Compute the product ofintra-cluster joins and inter-cluster non-joins
  prob <- prod(inter, intra)

  return(prob)
}


#' Extract the numerical representation of inter-cluster closest relationship
#'
#' Internal function on which compute_pr_IP_RG relies. Assumes clones, siblings and strangers
#' are represented by numerical values 1 > 0.5 > 0 respectively.
#'
#' @param IC1,IC2  Numeric vector of vertex ids for clusters one (IC1) and two (IC2).
#' @param RG_matrix Sparse matrix of relationships represented by numerical values
#' @noRd
extract_closest_relative <- function(IC1, IC2, RG_matrix) {
  relationships <- RG_matrix[IC1, IC2] # Inter-cluster relationships
  max(relationships) # Closest relationship
}

#' Compute inter-cluster probabilities
#'
#' Internal function on which compute_pr_IP_RG relies. Assumes siblings are IBD
#' with probability 0.5.
#'
#' @param inter_cluster_matrix Matrix of inter-cluster closest relatives
#' @noRd
compute_inter_cluster_join <- function(inter_cluster_matrix) {
  # Return zero probability if there are any inter-cluster clones
  if (any(inter_cluster_matrix == 1)) {
    return(0)
  } else {
    # Build a graph of inter-cluster sibling relationships
    inter_cluster_graph <- igraph::graph_from_adjacency_matrix(inter_cluster_matrix, mode = "lower", weighted = T)
    # Extract components of the inter-cluster sibling graph
    inter_cluster_components <- igraph::components(inter_cluster_graph)
    # Return zero probability if there are any components contain more than two disconnected sibling clusters
    if (any(inter_cluster_components$csize > 2)) {
      return(0)
    } else {
      # Compute probabilities within sibling components
      # (probabilities across sibling components are one in an outbred population)
      p <- (1 - 0.5)^sum(inter_cluster_components$csize == 2)
      return(p)
    }
  }
}

#' Compute intra-cluster probabilities
#'
#' Internal function on which compute_pr_IP_RG relies. Assumes the probability of IBD of given relationship is equal
#' to its numerical representation.
#'
#' @param IC numeric vector of intra-cluster vertex ids
#' @param inter_cluster_matrix Matrix of inter-cluster closest relatives
#'
#' @noRd
compute_intra_cluster_join <- function(IC, RG_matrix) {
  if (length(IC) > 1) { # If more than one genotype
    RC_matrix <- RG_matrix[IC, IC] # Within-cluster relationships
    RC_matrix[upper.tri(RC_matrix)] <- NA # Mask upper triangle for next step
    closest_relatives <- apply(RC_matrix, 1, max, na.rm = T)[-1] # Closest relative per added genotype relative to present
    p <- prod(closest_relatives) # Probability that genotypes were sequentially added
  } else {
    p <- 1
  }
  return(p)
}
