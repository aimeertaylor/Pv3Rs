# test_correct_graph() copied from the following link
# https://github.com/jwatowatson/RecurrentVivax/blob/master/Genetic_Model/iGraph_functions.R
# This function was written by James Watson and edited by Aimee Taylor

#' This implements the viable graph brute-force search algorithm described in the appendix of
#' Taylor & Watson et al. Nature communications 10.1 (2019).
#' As is, this function assumes there are only three relationships: strangers, siblings and clones
#' represented by values 0, 0.5 and 1, respectively. As such, all connected components of the relationship
#' graph must be cliques (in other words, the graph must be a cluster graph) and, among
#' cliques of size three, an edge weight sum of 2.5 is indicative of a intransitive clone + clone + sibling trio
test_transitive <- function(G){

  transitive <- TRUE # Presume correct until proven incorrect

  # First check that all connected components are cliques
  S <- igraph::components(graph = G)
  for (k in 1:S$no) {
    G_sub <- igraph::induced_subgraph(graph = G, vids = igraph::vertex_attr(G)$name[S$membership == k])
    MC <- igraph::max_cliques(graph = G_sub, min = igraph::vcount(G_sub))
    if (length(MC) == 0) transitive <- FALSE
  }

  # If all connected components are cliques, iterate through cliques of size 3
  if (transitive) {
    all_trios <- igraph::cliques(graph = G, min = 3, max = 3)
    K <- length(all_trios)
    if (K > 0) {
      i <- 1
      while (transitive & i <= K) { # stop searching when incorrect
        G_sub <- igraph::induced_subgraph(graph = G, vids = all_trios[[i]])
        score_sub <- sum(igraph::E(G_sub)$weight) # the sum of the edges tells us if it's correct
        if (score_sub == 2.5) transitive <- FALSE  # 2.5 is the only possible incorrect score so we reject
        i <- i + 1
      }
    }
  }

  return(transitive)
}
