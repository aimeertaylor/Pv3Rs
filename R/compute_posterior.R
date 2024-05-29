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
#' enumeration and likelihood calculation of relationship graphs can be found in
#' \code{\link{enumerate_RGs}} and \code{\link{RG_inference}} respectively.
#'
#' Model assumptions:
#' \itemize{
#'   \item{No within-host mutations, genotyping errors, undetected alleles}
#'   \item{Parasites are outbred}
#'   \item{Relationship graphs are equally likely given recurrence states}
#' }
#'
#' \code{compute_posterior} supports missing data, which should be encoded as
#' \code{NA}s. However, to avoid estimates based entirely on the weakly
#' informative nature of multiple per-marker allele calls (see example below and
#' `vignette("missing_data", "Pv3Rs")` for more details), we recommend against
#' generating estimates for recurrences that have no paired data due to
#' missingness (see Microsatellite data example in `vignette("demo", "Pv3Rs")`).
#'
#' The data input expects each list of alleles (for an infection) to consist of
#' a set of distinct alleles. Providing the function with a multiset of alleles
#' may lead to undefined behaviour. Note that the function will automatically
#' produce all assignments with repeated alleles.
#'
#' @param y Observed data in the form of a list of lists. Alleles should be
#'   provided as a set of distinct alleles for each infection (per marker).
#'   The number of entries
#'   is the number of episodes in chronological order. Episode names can be
#'   specified, but they are not used. Each episode is in turn a list of
#'   observed alleles for each marker, which must be named, or \code{NA} if not
#'   observed. For a given marker, alleles are modeled as categorical random
#'   variables. As such, allele names are arbitrary, but must correspond with
#'   frequency names (see examples below). The same names can be used for
#'   alleles belonging to different markers. As such, frequencies must be
#'   specified per named allele per named marker.
#' @param fs List of allele frequencies as vectors. Names of the list must match
#'   with the marker names in `y`. Within lists (i.e., for each marker),
#'   frequencies must be specified per allele name.
#' @param prior Matrix of prior probabilities of the recurrence states for each
#'   recurrent episode. Each row corresponds to an episode in chronological
#'   order. The column names must be C, L, and I for recrudescence, relapse and
#'   reinfection respectively. Row names can be specified by they are not used.
#'   If prior is not provided, a uniform prior will be used.
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
#'       \code{\link{enumerate_RGs}}.}
#'   }
#'
#' @examples
#'
#' # ===========================================================================
#' # Example where alleles are named arbitrarily
#' # ===========================================================================
#' # Data on an enrolment episode and a recurrence:
#' y <- list(episode0 = list(marker1 = c("Tinky Winky", "Dipsy"),
#'                           marker2 = c("Laa-Laa", "Po")),
#'           episode1 = list(marker1 = "Tinky Winky",
#'                           marker2 = "Laa-Laa"))
#'
#' # Allele frequencies:
#' fs <- list(
#'   marker1 = setNames(c(0.4, 0.6), c("Tinky Winky", "Dipsy")),
#'   marker2 = setNames(c(0.2, 0.8), c("Laa-Laa", "Po"))
#' )
#'
#' # Compute posterior probabilities using default uniform prior, note that
#' # since there is only one recurrence, the marginal probabilities are the same
#' # as the joint probabilities:
#' compute_posterior(y, fs)
#'
#'
#'
#' # ===========================================================================
#' # Example where alleles are given numeric names: 1, 2, 3, 4, 5
#' # ===========================================================================
#' # Allele frequencies with numeric names
#' fs <- list(m1=c('1'=0.78, '2'=0.14, '3'=0.07, '4'=0.005, '5' = 0.005),
#'            m2=c('1'=0.27, '2'=0.35, '3'=0.38),
#'            m3=c('1'=0.55, '2'=0.45))
#'
#' # Data:
#' y <- list(enrol = list(m1=c('3','2'), m2=c('1','3'), m3=c('1','2')),
#'           recur1 = list(m1=c('1','4'), m2=c('1','2'), m3=c('1','2')),
#'           recur2 = list(m1=c('1','5'), m2=c('2','3'), m3=c('1')))
#'
#' compute_posterior(y, fs)
#'
#'
#'
#' # ===========================================================================
#' # Example demonstrating the cosmetic-only nature of episode names. Due to the
#' # cosmetic-only nature of episode names, input info (episode data and prior)
#' # and output results should be chronologically ordered and interpreted.
#' # ===========================================================================
#' y <- list(enrol = list(m1 = NA),
#'           recur2 = list(m1 = NA),
#'           recur1 = list(m1 = NA))
#' prior <- array(c(0.2,0.7,0.2,0.3,0.6,0), dim = c(2,3),
#'                dimnames = list(c("recur1", "recur2"), c("C", "L", "I")))
#'
#' # The first prior row named "recur1" is returned for the first recurrence
#' # despite it being named "recur2" and the the second prior row named "recur2"
#' # is returned for the second recurrence despite it being named "recur1":
#' compute_posterior(y, fs, prior)
#' prior
#'
#'
#'
#' #============================================================================
#' # compute_posterior() returns the prior when there are no data. However, the
#' # prior will be re-weighted if NAs encode MOIs that are incompatible with
#' # recrudescence. (Recrudescing parasites are clones of parasites in the
#' # preceding blood-stage infection. The Pv3R model assumes no within-host
#' # mutations, genotyping errors or undetected alleles. As such, recrudescence
#' # is incompatible with an MOI increase relative to the preceding infection.)
#' #============================================================================
#' # Allele frequencies:
#' fs <- list(m1 = setNames(c(0.25, 0.75), c("A", "Other")))
#'
#' # Data on an enrolment episode and two recurrences:
#' y_no_data_MOIs111 <- list(enrol = list(m1 = NA),
#'                           recur1 = list(m1 = NA),
#'                           recur2 = list(m1 = NA))
#'
#' y_no_data_MOIs121 <- list(enrol = list(m1 = NA),
#'                           recur1 = list(m1 = c(NA, NA)),
#'                           recur2 = list(m1 = NA))
#'
#' # Compute posterior probabilities:
#' post3Rs_no_data_MOIs111 <- compute_posterior(y_no_data_MOIs111, fs)
#' post3Rs_no_data_MOIs121 <- compute_posterior(y_no_data_MOIs121, fs)
#'
#' # Returns the default prior since MOIs 111 are compatible with all 3R
#' # sequences:
#' post3Rs_no_data_MOIs111
#'
#' # Returns the default prior re-weighted to the exclusion of 3R sequences with
#' # a recrudescence at the first recurrence:
#' post3Rs_no_data_MOIs121
#'
#'
#'
#' #============================================================================
#' # Example demonstrating the weakly informative nature of multiple per-marker
#' # allele calls; see vignette("missing_data", "Pv3Rs") for more details.
#' #============================================================================
#' fs = list(m1 = c('1' = 0.5, '2' = 0.5))
#' y <- list(enrol = list(m1 = c('1', '2')), recur = list(m1 = NA))
#'
#' # compute_posterior() does not return the prior despite there being no
#' # recurrent data:
#' compute_posterior(y, fs)
#'
#'
#'
#' #============================================================================
#' # Example of the small but undesirable effect on the posterior of prior on
#' # graphs: the marginal probability that the first recurrence is a
#' # recrudescence increases as the number of recurrences increases even though
#' # only the first recurrence has data (also see vignette [to-do - base on
#' # MyDevFiles/Graph_prior_bias_examples.R])
#' #============================================================================
#' # Allele frequencies:
#' fs <- list(m1 = setNames(c(0.25, 1-0.25), c("A", "Other")))
#'
#' # Data for different scenarios; scenarios where the number of recurrences
#' # increases but only the first recurrence has data
#' ys <- list(y1 = list(enrol = list(m1 = "A"), recur1 = list(m1 = "A")),
#'            y2 = list(enrol = list(m1 = "A"), recur1 = list(m1 = "A"), recur2 = list(m1 = NA)),
#'            y3 = list(enrol = list(m1 = "A"), recur1 = list(m1 = "A"), recur2 = list(m1 = NA), recur3 = list(m1 = NA)),
#'            y4 = list(enrol = list(m1 = "A"), recur1 = list(m1 = "A"), recur2 = list(m1 = NA), recur3 = list(m1 = NA), recur4 = list(m1 = NA)))
#'
#' # Compute posterior probabilities and extract marginal probabilities:
#' results <- lapply(ys, function(y) compute_posterior(y, fs)$marg)
#'
#' # Extract results for the first recurrence only:
#' first_recur <- sapply(results, function(result) result[1,])
#'
#' # Plot 2D simplex
#' n_recur <- max(sapply(ys, length)-1)
#' pardefault <- par()
#' par(mar = c(0,0,0,0))
#' plot_simplex(v_labels = rownames(first_recur))
#'
#' # Project probabilities onto 2D simplex coordinates
#' xy <- apply(first_recur, 2, project2D)
#'
#' # Plot divergence from one recurrence to four:
#' arrows(x0 = xy["x", 1], x1 = xy["x", n_recur],
#'        y0 = xy["y", 1], y1 = xy["y", n_recur],
#'        length = 0.05, col = "red")
#'
#' # Plot a point for each recurrence from one to four:
#' points(x = xy["x", ], y = xy["y", ], pch = ".")
#'
#' # Restore plotting margins
#' par(mar = pardefault$mar)
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

  # Check priors sum to one
  if (!is.null(prior)) {
    if (!all(abs(1 - rowSums(prior)) < .Machine$double.eps^0.5)) {
      stop('Prior probabilities for a given recurrence must sum to one')
    }
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
  pbar <- msg_progress_bar(n.RG)
  message("Finding log-likelihood of each vector of recurrent states")
  max.logp <- max(sapply(RGs, "[[", "logp"))
  for (RG in RGs) {
    # subtract maximum to avoid underflow
    prob_RG <- exp(RG$logp - max.logp) # p(y|RG)
    for (rstr in compatible_rstrs(RG, gs_per_ts)) {
      n_rg_per_rstr[rstr] <- n_rg_per_rstr[rstr] + 1
      logp_sum_per_rstr[rstr] <- log(exp(logp_sum_per_rstr[rstr]) + prob_RG)
    }
    RG_i <- RG_i + 1
    pbar$increment()
  }
  message("")

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
