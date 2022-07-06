#' Internal function to test if a relationship graph is transitive.
#'
#' \code{is.transitive} implements the brute-force algorithm described in the
#' appendix of Taylor & Watson et al. Nature communications (2019). It assumes
#' there are only three relationships: strangers, siblings and clones,
#' represented by weights 0, 0.5 and 1, respectively, such that strangers are
#' disconnected, whereas siblings and clones are connected. Consequently, all
#' components of a transitive relationship graph must be cliques. Meanwhile,
#' among cliques of size three, an edge weight sum of 2.5 is indicative of a
#' intransitive clone + clone + sibling trio.
#'
#' @param RG A relationship graph; see \code{\link{enumerate_RGs}}.
#'
#' @section Provenance: This function was adapted from
#' test_correct_graph at
#' \url{https://github.com/jwatowatson/RecurrentVivax/blob/master/Genetic_Model/iGraph_functions.R}
#'
#' @noRd
is.transitive <- function(RG){

  transitive <- TRUE # Presume correct until proven incorrect

  # First check that all connected components are cliques
  S <- igraph::components(graph = RG)
  for (k in 1:S$no) {
    gs <- igraph::vertex_attr(RG)$name[S$membership == k]
    RG_sub <- igraph::induced_subgraph(graph = RG, vids = gs)
    MC <- igraph::max_cliques(graph = RG_sub, min = igraph::vcount(RG_sub))
    if (length(MC) == 0) transitive <- FALSE
  }

  # If all connected components are cliques, iterate through cliques of size 3
  if (transitive) {
    all_trios <- igraph::cliques(graph = RG, min = 3, max = 3)
    K <- length(all_trios)
    if (K > 0) {
      i <- 1
      while (transitive & i <= K) { # stop searching when incorrect
        RG_sub <- igraph::induced_subgraph(graph = RG, vids = all_trios[[i]])
        score_sub <- sum(igraph::E(RG_sub)$weight) # the sum of the edges tells us if it's correct
        if (score_sub == 2.5) transitive <- FALSE  # 2.5 is the only possible incorrect score so we reject
        i <- i + 1
      }
    }
  }

  return(transitive)
}
