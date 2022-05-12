#' Internal check function
#'
#' Check if Pr_yn_G > 0 due to a clonal edge between non-identical genotypes
#'
#' @section Provenance: This function was adapted from
#'   \code{test_cln_incompatible} at
#'   \url{https://github.com/jwatowatson/RecurrentVivax/blob/master/Genetic_Model/iGraph_functions.R}.
#'
test_cln_compatible <- function(A, vertex_data_matrix){

  cln_edges = which(A == 1, arr.ind = TRUE) # Which edges are clones
  num_cln_edges = nrow(cln_edges) # Number of edges that are clones

  if(num_cln_edges == 0){ # If there are no clones, return TRUE
    return(TRUE)
  } else { # For each clonal edge...
    for(i in 1:nrow(cln_edges)){
      alleles_1 <- vertex_data_matrix[cln_edges[i,'row'],]
      alleles_2 <- vertex_data_matrix[cln_edges[i,'col'],]
      score = all(alleles_1 == alleles_2, na.rm = T) # Doesn't count NAs as different
      if(!score){break()} # as soon as different, fail and break
    }
    return(score)
  }
}
