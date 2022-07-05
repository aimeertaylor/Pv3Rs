#' Sample a transitive relationship graph (RG)
#'
#' Samples one transitive graph of stranger and sibling relationships between
#' distinct parasite genotypes within infections and stranger, sibling and
#' clonal relationships between parasite genotypes across infection/s. Unlike
#' \code{\link{enumerate_RGs}}, \code{sample_RG} is not limited to six or fewer
#' genotypes among three or fewer infections, because it only generates one
#' graph and thus cannot overload memory. However, for more than six genotypes
#' among more than three infections, it could take a long time because
#' \code{sample_RG} is testing for transitivity while sampling among not
#' necessarily transitive graphs, of which there will be many; see the message
#' that \code{sample_RG} prints to the screen.
#'
#' @inheritParams enumerate_RGs
#'
#' @return Returns a transitive relationship graph (RG); see \code{\link{enumerate_RGs}} for details.
#'
#' @section Provenance: This function was adapted from
#'   \code{generate_random_Vivax_model} at
#'   \url{https://github.com/jwatowatson/RecurrentVivax/blob/master/Genetic_Model/iGraph_functions.R}.
#'
#'
#' @examples
#'
#' RG <- sample_RG(c(2,1,1))
#' print(RG)
#' plot_RG(RG)
#'
#' @export
sample_RG = function(MOIs){

  # Check MOIs are positive whole numbers
  if (!all(is.wholenumber(MOIs)) | any(MOIs < 1)) stop("MOIs should be positive integers")

  # Compute the number of not necessarily transitive graphs
  intra_edge_counts <- sapply(MOIs, choose, k = 2)
  inter_edge_count <- choose(sum(MOIs), 2) - sum(intra_edge_counts)
  RGs_to_eval_count <- prod(3^inter_edge_count, 2^intra_edge_counts)
  if(RGs_to_eval_count > 1000){
    writeLines(sprintf("\nThis might take a while: sampling among %s not necessarily transitive RGs",
                       RGs_to_eval_count))
  }

  # Hard code relationship types to satisfy the is.transitive function
  relationship_types <- c(stranger = 0, sibling = 0.5, clone = 1)
  intra_relationship_types <- relationship_types[setdiff(names(relationship_types), "clone")]
  infection_count <- length(MOIs) # Number of time points
  gs_count <- sum(MOIs) # Number of genotypes

  # Check if feasible to enumerate RGs
  if (gs_count <= 1) stop("Sorry, no RGs")

  # Enumerate indices for pairs of time points, etc.
  if (infection_count > 1) {
    infection_ijs <- gtools::combinations(n = infection_count,
                                           r = 2,
                                           v = 1:infection_count)
    infection_pair_count <- nrow(infection_ijs) # Number of pairs of time points
  } else {
    infection_pair_count <- 0 # Number of pairs of time points
  }

  gs <- paste0("g", 1:gs_count) # Genotype names
  ts_per_gs <- rep(1:infection_count, MOIs) # List time point indices per genotype

  # List genotypes per time point and time point pair
  gs_per_ts <- split(gs, ts_per_gs)
  if(infection_pair_count >= 1){
    gs_per_t_pairs <- lapply(1:infection_pair_count, function(t_pair) {
      list(gs_per_ts[[infection_ijs[t_pair, 1]]],
           gs_per_ts[[infection_ijs[t_pair, 2]]])
    })
  }


  # Allocate adjacency matrices
  # Preallocation adj_all as follows is preferable to
  # adj_all <- Matrix::bdiag(intra_adjs) because Matrix::bdiag does not inherit
  # dimnames and pre-populates inter-infection relationships with zero not NA
  # making it appear that relationships have been sampled when they haven't
  adj_all <- array(NA, dim = rep(gs_count, 2), dimnames = list(gs,gs))
  intra_adjs <- lapply(MOIs, function(MOI) array(NA, dim = rep(MOI, 2)))

  transitive <- FALSE
  while (!transitive) {

    # Populate adj_all with intra-relationships
    for(t in 1:infection_count) {
      intra_relationships <- sample(intra_relationship_types, choose(MOIs[t],2), replace = T)
      intra_adjs[[t]][lower.tri(intra_adjs[[t]])] <- intra_relationships
      adj_all[gs_per_ts[[t]], gs_per_ts[[t]]] <- intra_adjs[[t]]
    }

    # Populate lower tri of all_adj with between time-point relationships
    if (infection_pair_count >= 1) {
      for (t_pair in 1:infection_pair_count) {
        inter_relationship_count <- prod(MOIs[infection_ijs[t_pair,]])
        inter_relationships <- sample(relationship_types, inter_relationship_count, replace = T)
        adj_all[gs_per_t_pairs[[t_pair]][[2]],
                gs_per_t_pairs[[t_pair]][[1]]] <- inter_relationships
      }
    }

    # Convert adjacency matrix into an igraph item
    RG <- igraph::graph_from_adjacency_matrix(adjmatrix = adj_all,
                                              mode = "lower",
                                              diag = F,
                                              weighted = T)
    # Add a time-point attribute
    RG <- igraph::set_vertex_attr(RG, "group", value = ts_per_gs)

    # Test transitivity and store if transitive
    transitive <- is.transitive(RG)
  }
  return(RG)
}

