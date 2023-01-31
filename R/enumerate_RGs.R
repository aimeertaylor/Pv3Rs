#' Enumerate transitive relationship graphs (RGs)
#'
#' Enumerates all transitive graphs of stranger and sibling relationships
#' between distinct parasite genotypes within an infection and stranger, sibling
#' and clonal relationships between parasite genotypes across infections.
#' \code{enumerate_RGs} is limited to six or fewer genotypes among three or
#' fewer infections because the number of graphs grows exponentially with the
#' number of genotypes and the number of infections; see examples.
#'
#' @param MOIs A numeric vector specifying, for each infection, the number of
#'   distinct parasite genotypes, a.k.a the multiplicity of infection (MOI).
#'
#' @return Returns a list of transitive relationship graphs (RGs).
#'
#'   Each RG is an igraph graph, see \code{\link[igraph]{igraph}}. The seven
#'   character code after IGRAPH, is the first seven characters of the graph ID;
#'   see \code{\link[igraph]{graph_id}}.
#'
#'   RGs are undirected. They are weighted, except when they contain only stranger genotypes.
#'
#'   Vertices represent parasite genotypes. Vertices are named consecutively
#'   from one to the total genotype count, \code{sum(MOIs)}. Vertex names are
#'   used to index genotypes and inter-genotype relationships. Vertex IDs are
#'   not used to index genotypes and their relationships because vertex IDs are
#'   renumbered upon application of some igraph operations, e.g.
#'   \code{\link[igraph]{induced_subgraph}}. Genotypes are grouped into
#'   infections using a vertex attribute called group. The group attribute is
#'   used by \code{\link{plot_RG}} to group genotypes by infection.
#'
#'   Edges represent sibling or clonal relationships that are differentiated by
#'   weight: 0.5 encodes a sibling relationship, 1 encodes a clonal
#'   relationship. Strangers have zero weight and thus there are no edges
#'   between them.
#'
#' @section Provenance: This function was adapted from
#'   \code{generate_all_models_3Ts} and \code{generate_all_models_2Ts} at
#'   \url{https://github.com/jwatowatson/RecurrentVivax/blob/master/Genetic_Model/iGraph_functions.R}.
#'
#' @section To-do:
#' \enumerate{
#' \item Consider using RGs_to_eval_count as the too
#'   many RGs cut-off instead of gs_count > 6 | infection_count > 3
#' \item Read up about using save and load to preserve graph attributes; also
#' read_graph and write_graph
#' }
#'
#'
#' @examples
#' MOIs <- c(1, 2)
#' RGs <- enumerate_RGs(MOIs)
#' RGs[[1]]
#' igraph::vertex_attr(RGs[[1]])
#' par(mfrow = c(3, 3))
#' for (i in 1:length(RGs)) {
#'   plot_RG(RGs[[i]])
#'   box()
#' }
#'
#' # Compute by hand the number of not necessarily transitive graphs
#' intra_edge_counts <- sapply(MOIs, choose, k = 2)
#' inter_edge_count <- choose(sum(MOIs), 2) - sum(intra_edge_counts)
#' prod(3^inter_edge_count, 2^intra_edge_counts)
#'
#' @export
enumerate_RGs <- function(MOIs) {
  # Check MOIs are positive whole numbers
  if (!all(is.wholenumber(MOIs)) | any(MOIs < 1)) stop("MOIs should be positive integers")

  # Compute the number of not necessarily transitive graphs by hand
  intra_edge_counts <- sapply(MOIs, choose, k = 2)
  inter_edge_count <- choose(sum(MOIs), 2) - sum(intra_edge_counts)
  RGs_to_eval_count <- prod(3^inter_edge_count, 2^intra_edge_counts)

  # Hard code relationship types to satisfy the is.transitive function
  relationship_types <- c(stranger = 0, sibling = 0.5, clone = 1)
  intra_relationship_types <- relationship_types[setdiff(names(relationship_types), "clone")]
  infection_count <- length(MOIs) # Number of time points
  gs_count <- sum(MOIs) # Number of genotypes

  if (gs_count > 6 | infection_count > 3) stop("Sorry, too many RGs")
  if (gs_count <= 1) stop("Sorry, no RGs")

  # Enumerate indices for pairs of time points, etc.
  if (infection_count > 1) {
    infection_ijs <- gtools::combinations(
      n = infection_count,
      r = 2,
      v = 1:infection_count
    )
    infection_pair_count <- nrow(infection_ijs) # Number of pairs of time points
  } else {
    infection_pair_count <- 0 # Number of pairs of time points
  }

  gs <- paste0("g", 1:gs_count) # Genotype names
  ts_per_gs <- rep(1:infection_count, MOIs) # List time point indices per genotype

  # List genotypes per time point and time point pair
  gs_per_ts <- split(gs, ts_per_gs)
  if (infection_pair_count >= 1) {
    gs_per_t_pairs <- lapply(1:infection_pair_count, function(t_pair) {
      list(
        gs_per_ts[[infection_ijs[t_pair, 1]]],
        gs_per_ts[[infection_ijs[t_pair, 2]]]
      )
    })
  }

  # Enumerate all combinations of intra time-point relationships
  intra_relationships <- lapply(MOIs, function(MOI) {
    if (MOI > 1) {
      gtools::permutations(
        n = length(intra_relationship_types),
        r = choose(MOI, 2),
        v = intra_relationship_types,
        repeats.allowed = T
      )
    } else {
      matrix(NA, 1, 1)
    }
  })

  # Enumerate all combinations of inter time-point relationships
  if (infection_count > 1) {
    inter_relationships <- lapply(1:infection_pair_count, function(t_pair) {
      as.matrix(expand.grid(rep(list(relationship_types), prod(MOIs[infection_ijs[t_pair, ]]))))
    })
    inter_count <- sapply(inter_relationships, nrow) # Numbers of inter_relationships to permute
  } else {
    inter_count <- c()
  }

  # Compute the number of not necessarily transitive relationship graphs
  intra_count <- sapply(intra_relationships, nrow) # Numbers of intra_relationships to permute
  all_counts <- c(intra_count, inter_count)
  all_perms <- as.matrix(expand.grid(sapply(all_counts, function(x) 1:x))) # Matrix of permutations
  total_count <- prod(all_counts) # Total number of permutations
  if (nrow(all_perms) != total_count) {
    stop("Problem with the not necessarily transitive RG count")
  }
  writeLines(paste("\nnumber of not necessarily transitive graphs is", total_count))

  # Allocate adjacency matrices
  adj_all <- array(NA, dim = rep(gs_count, 2), dimnames = list(gs, gs))
  intra_adjs <- lapply(MOIs, function(MOI) array(NA, dim = rep(MOI, 2)))

  # Create progress bar and count of transitive RGs
  pbar <- txtProgressBar(min = 1, max = total_count)
  transitive_RGs <- list()
  transitive_i <- 1

  for (perm_i in 1:total_count) {
    setTxtProgressBar(pbar, perm_i)

    # Extract permutation indices
    intra_perm <- all_perms[perm_i, 1:infection_count]
    inter_perm <- all_perms[perm_i, -(1:infection_count)]

    # Populate adj_all with intra-relationships
    for (t in 1:infection_count) {
      intra_adjs[[t]][lower.tri(intra_adjs[[t]])] <- intra_relationships[[t]][intra_perm[t], ]
      adj_all[gs_per_ts[[t]], gs_per_ts[[t]]] <- intra_adjs[[t]]
    }

    # Populate lower tri of all_adj with between time-point relationships
    if (infection_pair_count >= 1) {
      for (t_pair in 1:infection_pair_count) {
        adj_all[
          gs_per_t_pairs[[t_pair]][[2]],
          gs_per_t_pairs[[t_pair]][[1]]
        ] <- inter_relationships[[t_pair]][inter_perm[t_pair], ]
      }
    }

    # Convert adjacency matrix into an igraph item
    RG <- igraph::graph_from_adjacency_matrix(
      adjmatrix = adj_all,
      mode = "lower",
      diag = F,
      weighted = T
    )
    # Add a time-point attribute
    RG <- igraph::set_vertex_attr(RG, "group", value = ts_per_gs)

    # Test transitivity and store if transitive
    if (is.transitive(RG)) {
      transitive_RGs[[transitive_i]] <- RG
      transitive_i <- transitive_i + 1
    }
  }

  writeLines(paste("\nnumber of transitive graphs is", length(transitive_RGs)))
  return(transitive_RGs)
}
