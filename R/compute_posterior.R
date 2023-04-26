#' Compute the posterior distribution of recurrence states
#'
#' Entry point to Bayesian inference for \emph{P. vivax} recurrence states based
#' on genetic data. Specifically, this function finds the posterior
#' probabilities of relapse, reinfection and recrudescence from genetic data.
#'
#' We enumerate all possible relationship graphs between genotypes, where each
#' pair of genotypes may be clones, siblings, or strangers, each with a
#' different level of expected genetic relatedness. The likelihood of a
#' sequence of recurrence states can be determined from the likelihood of all
#' relationship graphs compatible with said sequence. More details on the
#' enumeration and likelihood calculation of relationship graphs
#' can be found in \code{\link{enumerate_RGs_alt}} and
#' \code{\link{RG_inference}} respectively.
#'
#' Model assumptions:
#' \itemize{
#'   \item{No genotyping errors or missing alleles}
#'   \item{Parasites are outbred}
#'   \item{Relationship graphs are equally likely given recurrence states}
#' }
#'
#' @param y Observed data in the form of a list of lists. The number of entries is the
#'   number of episodes. Each episode is in turn a list of observed alleles
#'   for each marker, or \code{NA} if not observed.
#' @param fs List of allele frequencies as vectors. Names of the list must
#'   match with the marker names in `y`.
#' @param prior Matrix of prior probabilities of the recurrence states for each
#'   recurrent episode. Each row corresponds to an episode. The column names
#'   must be C, L, and I for recrudescence, relapse and reinfection
#'   respectively. If prior is not provided, a uniform prior will be used.
#' @param return.RG Boolean for whether to return the relationship graphs,
#'   defaults to `FALSE`.
#'
#' @return List containing:
#'   \describe{
#'     \item{marg}{Matrix of marginal posterior probabilities of the possible
#'       recurrence states for each recurrent episode, one row per reinfection.}
#'     \item{joint}{Vector of joint posterior probabilities of each possible string
#'       of recurrence states.}
#'     \item{RGs}{List of relationship graphs with their log-likelihoods stored.
#'       Only returned if `return.RG` is `TRUE`. See
#'       \code{\link{enumerate_RGs_alt}}.}
#'   }
#'
#' @examples
#' # Data on an enrollment episode and a recurrence:
#' y <- list(list(m1 = c("A", "C"), m2 = c("G", "T")), list(m1 = "A", m2 = "G"))
#'
#' # Allele frequencies:
#' fs <- list(
#'   m1 = setNames(c(0.4, 0.6), c("A", "C")),
#'   m2 = setNames(c(0.2, 0.8), c("G", "T"))
#' )
#'
#' # Compute posterior probabilities using default uniform prior, note that
#' # since there is only one recurrence, the marginal probabilities are the same
#' # as the joint probabilities:
#' compute_posterior(y, fs)
#'
#'
#'
#' #===============================================================================
#' # compute_posterior() returns the prior when there are either no data whatsoever
#' # or no recurrent data
#' #===============================================================================
#'
#' # Allele frequencies:
#' fs <- list(m1 = setNames(c(0.25, 0.75), c("A", "Other")))
#'
#' # Data on an enrollment episode and two recurrences:
#' ys_no_data_whatsoever <- list(list(m1 = NA), list(m1 = NA), list(m1 = NA))
#' ys_no_recurrent_data <- list(list(m1 = "A"), list(m1 = NA), list(m1 = NA))
#'
#' # Prior on two recurrences:
#' prior <- array(c(0.4,0.4,0.2,0.6,0.4,0), dim = c(2,3),
#'                dimnames = list(NULL, c("C", "L", "I")))
#'
#' # Compute posterior probabilities:
#' post3Rs_no_data_whatsover <- compute_posterior(ys_no_data_whatsoever, fs, prior)
#' post3Rs_no_recurrent_data <- compute_posterior(ys_no_recurrent_data, fs, prior)
#'
#' # Compare prior with posterior marginal probabilities:
#' prior
#' post3Rs_no_data_whatsover$marg
#' post3Rs_no_recurrent_data$marg
#'
#' @export
compute_posterior <- function(y, fs, prior = NULL, return.RG = FALSE) {

  # Check y is a list of lists:
  if (class(y) != "list" | unique(unlist(lapply(y, class))) != "list") {
    stop("Data y must be a list of lists, even if only one marker is typed per infection")
  }

  # Check frequencies sum to one:
  if (!all(abs(1 - round(sapply(fs, sum))) < .Machine$double.eps^0.5)) {
    stop('For a given marker, allele frequencies must sum to one')
  }

  # Check all alleles have a named frequency
  ms <- unique(as.vector(sapply(y, names))) # Extract marker names
  # For each marker, list allele names:
  as <- lapply(ms, function(m) unique(as.vector(unlist(sapply(y, function(yt) yt[[m]])))))
  names(as) <- ms # Name list entries by marker
  # For each marker, check all alleles have a named frequency:
  all_got <- all(sapply(ms, function(m) all(as[[m]][!is.na(as[[m]])] %in% names(fs[[m]]))))
  if (!all_got) stop("Not all alleles have a named frequency")


  infection_count <- length(y)
  stopifnot(infection_count > 1)

  MOIs <- determine_MOIs(y)
  gs_count <- sum(MOIs)
  gs <- paste0("g", 1:gs_count)
  ts_per_gs <- rep(1:infection_count, MOIs)
  gs_per_ts <- split(gs, ts_per_gs)

  n_recur <- infection_count - 1
  causes <- c("C", "L", "I")

  # use uniform prior if prior not given
  if (is.null(prior)) {
    prior <- matrix(rep(1 / 3, 3 * (infection_count - 1)),
      ncol = 3
    )
    colnames(prior) <- causes
  } else {
    if (!identical(colnames(prior), causes)) stop('The prior should have colnames "C", "L" and "I"')
  }


  # For each infection we find the possible allele assignments
  # Because of permutation symmetry of the genotypes within an infection, we
  # fix the assignment of one marker's alleles, where the number of alleles for
  # that marker must be equal to the MOI
  # Result is list (over infections) of lists (over markers) of dataframes
  # Dataframe columns correspond to genotypes, rows correspond to assignments
  alleles_per_inf_per_m <- lapply(
    1:infection_count,
    function(t) { # for each infection
      enumerate_alleles(
        y[[t]],
        gs_per_ts[[t]]
      )
    }
  )

  # Take Cartesian product over infections to get allele assignments over all
  # genotypes, but still separated for each marker (list of dataframes)
  # No need to take the Cartesian product over markers due to outbred assumption
  alleles_per_m <- Reduce( # Reduce performs below repeatedly across infections
    # returns one list of merged dataframes given two lists of dataframes
    function(x, y) {
      mapply(merge, x, y,
        MoreArgs = list(by = NULL, all = T),
        SIMPLIFY = F
      )
    },
    alleles_per_inf_per_m
  )

  # enumerate all relationship graphs (RGs) and find their likelihood
  RGs <- RG_inference(MOIs, fs, alleles_per_m)

  # get all vectors of recurrence states
  n_rstrs <- 3^n_recur
  rstrs <- do.call(paste0, expand.grid(lapply(
    1:n_recur,
    function(x) causes
  )))

  # number of RGs consistent with each vector of recurrence states (R)
  n_rg_per_rstr <- setNames(rep(0, n_rstrs), rstrs)
  # for each R, stores log sum of p(y|RG) across RGs compatible with R
  logp_sum_per_rstr <- setNames(rep(-Inf, n_rstrs), rstrs)

  # for each RG, add p(y|RG) for recurrence states where RG is compatible
  RG_i <- 0
  n.RG <- length(RGs)
  pbar <- txtProgressBar(min = 0, max = n.RG)
  writeLines("Finding log-likelihood of each vector of recurrent states")
  for (RG in RGs) {
    prob_RG <- exp(RG$logp) # p(y|RG)
    for (rstr in compatible_rstrs(RG, gs_per_ts)) {
      n_rg_per_rstr[rstr] <- n_rg_per_rstr[rstr] + 1
      logp_sum_per_rstr[rstr] <- log(exp(logp_sum_per_rstr[rstr]) + prob_RG)
    }
    RG_i <- RG_i + 1
    setTxtProgressBar(pbar, RG_i)
  }
  writeLines("\n")

  # logp_per_rstr is now log of p(y|R) = sum of p(y|RG)*p(RG|R)
  # assume p(RG|R) is uniform
  logp_per_rstr <- logp_sum_per_rstr - log(n_rg_per_rstr)
  logp_per_rstr[is.nan(logp_per_rstr)] <- -Inf

  # add log prior so that logp_per_rstr is log of p(y,R) = p(y|R)p(R)
  logprior <- log(prior)
  for (rstr in rstrs) {
    logprior_rstr <- sum(sapply(
      1:n_recur,
      function(i) logprior[i, substr(rstr, i, i)]
    ))
    logp_per_rstr[rstr] <- logp_per_rstr[rstr] + logprior_rstr
  }
  # normalise p(y,R) to get posterior p(R|y)
  post_per_rstr <- exp(logp_per_rstr) / exp(matrixStats::logSumExp(logp_per_rstr))

  # find marginal probabilities
  marg <- matrix(NA, nrow = n_recur, ncol = 3)
  colnames(marg) <- causes
  for (i in 1:n_recur) { # for each recurrence
    # split returns a list of 3 vectors, each containing indices whose
    # corresponding R has i-th recurrence being C/L/I
    idx_list <- split(
      1:n_rstrs, # each int is an int representation of R
      # get i-th ternary digit and map to 1/2/3
      (1:n_rstrs - 1) %/% 3^(i - 1) %% 3 + 1
    )
    for (j in 1:3) {
      marg[i, j] <- sum(post_per_rstr[idx_list[[j]]])
    }
  }

  rownames(marg) <- names(y)[-1] # drop first infection (not reinfection)

  result <- list(marg = marg, joint = post_per_rstr)
  if (return.RG) result$RGs <- RGs
  result
}
