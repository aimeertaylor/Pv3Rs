#' Internal computation function
#'
#'
#' Notes copied from original:
#' Function to sum over vertex haploid genotype labels given edge kinship
#' labelled graph Returns: Pr(yn | Gnb) unnormalised = sum from a = 1 to A over
#' Pr(yn | Gnab) Unnormalised because  Pr(yn | Gnb) = sum from a = 1 to A over
#' Pr(yn | Gnab) * (1/A)
#'
#'
#' @section Provenance: This function was adapted from \code{Log_Pr_yn_Gnb_unnormalised} at
#'   \url{https://github.com/jwatowatson/RecurrentVivax/blob/master/Genetic_Model/iGraph_functions.R}.
#'

Log_Pr_yn_Gnb_unnormalised = function(G, Gabs, cn, Tn, log_Fs, MSs, alpha_terms){

  sum_cn = sum(cn)

  # Get adjacency matrix for Gnb
  if(is.null(igraph::edge_attr(G)$weight)){ # If all edges are strangers generate A afresh
    A = array(0, dim = c(sum_cn, sum_cn))
  } else { # Extract from A from graph
    A = as.matrix(igraph::get.adjacency(G, attr = 'w'))
  }
  A[upper.tri(A, diag = TRUE)] = NA # Work with lower triangular only

  # Sum over all Gna given adjacency matrix for Gnb
  log_Pr_yn_Gabs = sapply(Gabs, function(vertex_data_matrix){ # labelled_G_ind inc. compatible labelled graphs only
    Log_Pr_yn_Gab(G, log_Fs, MSs, Z = alpha_terms, A = A, vertex_data_matrix) # Pass Gn to Log_Pr_yn_Gab
  })

  log_pr_yn_Gb_unnormalised = logSumExp(log_Pr_yn_Gabs, na.rm = TRUE)
  return(log_pr_yn_Gb_unnormalised)
}


