#' Internal computation function featuring in Log_Pr_yn_Gnb_unnormalised
#'
#' Notes copied from original:
#' MS where one or more NA are ignored when summing over m
#' G is assumed labelled
#' Z is a vector of alpha_terms (calculated above since not data dependent
#'
#' @section Provenance: This function was adapted from
#'   \code{Log_Pr_yn_Gab} at
#'   \url{https://github.com/jwatowatson/RecurrentVivax/blob/master/Genetic_Model/iGraph_functions.R}.
#'
#' @noRd
Log_Pr_yn_Gab <- function(G, log_Fs, MSs, Z, A, vertex_data_matrix){

  M <- length(MSs) # Number of microsatellites

  # Check if data have zero prob given G due to clones
  if (!is.clone.compatible(A, vertex_data_matrix)) { # If labelled G is incompatible

    log_Pr_yn_Gab = NA
    return(log_Pr_yn_Gab) # Function ends here

  } else {

    # Calculate non-zero Pr( yn | G )
    log_eij = array(NA, dim = dim(A))

    for(i in 2:nrow(A)){ # Over all edges
      for(j in 1:(i-1)){

        log_eijm = rep(NA, M) # Store
        hi = vertex_data_matrix[i,] # Extract h_im
        hj = vertex_data_matrix[j,] # Extract h_jm
        I_hj = hi == hj # Marker identity indicator

        # Stranger edge indicator
        if(A[i,j] == 0){
          for(m in 1:M){ # For each microsatellite m = 1:M
            log_fs = log_Fs[[ MSs[m] ]][c(hi[m], hj[m])] # Extract the log frequencies for h_im and h_jm
            LX = c(Z['log_alpha'] + log_fs[1], Z['log_1_alpha'] + sum(log_fs)) # log(x) vector for logSumExp
            if(Z['alpha'] == 0){ # Avoid logSumExp as log(alpha = 1) = -Inf will overwhelm
              log_eijm[m] = sum(log_fs)
            } else {
              log_eijm[m] = ifelse(I_hj[m], logSumExp(LX), sum(Z['log_1_alpha'], log_fs))
            }
          }
        }

        # Sibling edge indicator
        if(A[i,j] == .5){
          for(m in 1:M){ # For each microsatellite m = 1:M
            log_fs = log_Fs[[ MSs[m] ]][c(hi[m], hj[m])] # Extract the log frequencies for h_im and h_jm
            LX = c(Z['log_0.5alpha'] + log_fs[1], Z['log_0.5_alpha'] + sum(log_fs)) # log(x) vector for logSumExp
            log_eijm[m] = ifelse(I_hj[m], logSumExp(LX), sum(Z['log_0.5_alpha'],log_fs))
          }
        }

        # Clone edge indicator
        if(A[i,j] == 1){
          for(m in 1:M){ # For each microsatellite m = 1:M
            log_fs = log_Fs[[ MSs[m] ]][c(hi[m], hj[m])] # Extract the log frequencies for h_im and h_jm
            log_eijm[m] = ifelse(I_hj[m], log_fs[1], NA)
          }
        }
        log_eij[i,j] = sum(log_eijm, na.rm = TRUE) # Calculate Pr_hij_eij
      }
    }
    log_Pr_yn_Gab = sum(log_eij, na.rm = TRUE)
    return(log_Pr_yn_Gab)
  }
}


