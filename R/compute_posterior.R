#' Compute posterior probabilities of \emph{P. vivax} recurrent states
#'
#' @description
#' Compute posterior probabilities of \emph{P. vivax} recurrent states relapse,
#' reinfection and recrudescence using genetic data.
#'
#' Please note, the progress bar does not necessarily increment uniformly (see
#' details below); it may seem stuck when the code is still running.
#'
#' @details
#' \code{compute_posterior()} computes posterior probabilities proportional to the likelihood multiplied by the prior. The likelihood
#' sums over
#'
#' - ways to phase allelic data onto haploid genotypes
#' - graphs of relationships between haploid genotypes
#' - ways to partition alleles into clusters of identity-by-descent
#'
#' \code{compute_posterior()} expects each per-episode, per-marker allelic
#' vector to be a set of distinct alleles. Allele repeats at markers with
#' observed data, and \code{NA} repeats at markers with missing data, are
#' removed in a data pre-processing step. \code{NA}s in allelic vectors that
#' also contain non-\code{NA} values are removed in a data pre-processing step.
#'
#' We enumerate all possible relationship graphs between haploid genotypes,
#' where pairs of genotypes can either be clones, siblings, or strangers. The
#' likelihood of a sequence of recurrence states can be determined from the
#' likelihood of all relationship graphs compatible with said sequence. More
#' details on the enumeration and likelihood calculation of relationship graphs
#' can be found in \code{\link{enumerate_RGs}} and \code{\link{RG_inference}}
#' respectively. For each relationship graph, the model sums over all possible
#' identity-by-descent partitions. Because some relationship graphs are
#' compatible with more identity-by-descent partitions than others, the log
#' p(Y|RG) progress bar does not necessarily increment uniformly.
#'
#' Notable model assumptions and limitations:
#' \itemize{
#'   \item{Perfect detection of alleles (no genotyping error)}
#'   \item{No within-host \emph{de novo} mutations}
#'   \item{Parasites are outbred}
#'   \item{All siblings are regular siblings}
#'   \item{Relationship graphs compatible with a given sequence of recurrent
#'   states are equally likely \emph{a priori}}
#'   \item{We do not recommend running `compute_posterior() when the total
#'   genotype count (sum of per-episode multiplicities of infection) exceeds
#'   eight, because there are too many relationship graphs.}
#'   \item{Presently, \code{Pv3Rs} only supports prevalence data (categorical data that
#' signal the detection of alleles), not quantitative data (data that signal the proportional
#' abundance of the alleles detected).}
#' }
#'
#'
#' @param y Observed data in the form of a list of lists. The outer list is a
#'   list of episodes in increasing chronological order. The inner list is a list of named
#'   markers per episode. Episode names can be specified, but they are not used.
#'   Markers must be named. Each episode must list the same markers. If not all
#'   markers are typed per episode, data on untyped markers can be encoded as
#'   missing (see below). For each marker, one must specify an allelic vector: a
#'   set of distinct alleles detected at that marker. \code{NA}s encode missing
#'   per-marker data, i.e., when no alleles are observed for a given marker.
#'   \code{NA} entries in allelic vectors that contain both \code{NA} and
#'   non-\code{NA} entries are ignored. Allele names are arbitrary, but must
#'   correspond with frequency names (see examples below). The same names can be
#'   used for alleles belonging to different markers. As such, frequencies must
#'   be specified per named allele per named marker.
#' @param fs List of per-marker allele frequency vectors. Names of the list must match
#'   with the marker names in `y`. Within lists (i.e., for each marker),
#'   frequencies must be specified per allele name.
#' @param prior Matrix of prior probabilities of the recurrence states for each
#'   recurrent episode. Each row corresponds to an episode in increasing
#'   chronological order. The column names must be C, L, and I for
#'   recrudescence, relapse and reinfection respectively. Row names can be
#'   specified but they are not used. If `prior` is NULL (default), per-episode
#'   recurrent states are equally likely.
#' @param MOIs Multiplicity of infection for each episode. If MOIs are not
#'   provided, the most parsimonious MOIs compatible with the data will be used;
#'   see \code{\link{determine_MOIs}}.
#' @param return.RG Boolean for whether to return the relationship graphs,
#'   defaults to `FALSE`. If `return.logp` is set to `TRUE`, then `return.RG`
#'   is overridden to be `TRUE`, as log-probabilities are returned for each
#'   relationship graph.
#' @param return.logp Boolean for whether to return the log-likelihood for each
#'   relationship graph, defaults to `FALSE`. When setting `return.logp` to
#'   `TRUE`, `return.RG` should also be set to `TRUE`. Setting `return.logp` to
#'   `FALSE` allows for permutation symmetries to be exploited to save
#'   computational time, see \code{\link{enumerate_alleles}}. Setting this to
#'   `TRUE` will result in longer runtimes, especially when multiplicities of
#'   infection are large. Note that this argument does not affect the output of
#'   the posterior probabilities.
#'
#' @return List containing:
#'   \describe{
#'     \item{`marg`}{Matrix of marginal posterior probabilities of
#'       recurrent states for each recurrence, one row per recurrence with "C"
#'       for recrudescence, "L" for relapse, and "I" for reinfection. `marg` is
#'       a simple summary of `joint` (see next): each marginal probability of a
#'       recurrent state is a sum over a subset of joint probabilities of
#'       recurrent state sequences. For example, the marginal probability of "C"
#'       at the first of two recurrences is a sum over the joint probabilities
#'       of "CC",
#'       "CL", and "CI".}
#'     \item{`joint`}{Vector of joint posterior probabilities of each recurrent
#'       state sequence, where within a sequence "C" denotes recrudescence,
#'       "L" denotes relapse, and "I" denotes reinfection.}
#'     \item{`RGs`}{List of relationship graphs with their log-likelihoods stored.
#'       Only returned if `return.RG` is `TRUE`. See
#'       \code{\link{enumerate_RGs}}.}
#'   }
#'
#' @examples
#' # ===========================================================================
#' # Example where alleles are named numerically
#' # ===========================================================================
#' # Data
#' y <- list(enroll = list(m1 = c('3','2'), m2 = c('1','2')),
#'           recur1 = list(m1 = c('1','4'), m2 = c('1')),
#'           recur2 = list(m1 = c('1'), m2 = NA))
#'
#' # Allele frequencies
#' fs <- list(m1 = c('1' = 0.78, '2' = 0.14, '3' = 0.07, '4' = 0.01),
#'            m2 = c('1' = 0.27, '2' = 0.73))
#'
#' # Compute posterior probabilities using default prior
#' compute_posterior(y, fs)
#'
#'
#' # ===========================================================================
#' # Example where alleles are named arbitrarily and probabilities are plotted
#' # ===========================================================================
#' # Data
#' y <- list(episode0 = list(marker1 = c("Tinky Winky", "Dipsy"),
#'                           marker2 = c("Tinky Winky", "Laa-Laa", "Po")),
#'           episode1 = list(marker1 = "Tinky Winky",
#'                           marker2 = "Laa-Laa"))
#'
#' # Allele frequencies
#' fs <- list(marker1 = c("Tinky Winky" = 0.4, "Dipsy" = 0.6),
#'            marker2 = c("Tinky Winky" = 0.1, "Laa-Laa" = 0.1, "Po" = 0.8))
#'
#' # Compute posterior probabilities using default prior
#' posterior_probs <- compute_posterior(y, fs)
#'
#' # Plot posterior probabilities on the simplex
#' plot_simplex(c("Recrudescence", "Relapse", "Reinfection"), 0.5) # Simplex
#' xy <- project2D(posterior_probs$marg[1,]) # Project probabilities
#' points(x = xy["x"], y = xy["y"], pch = 20) # Plot projected probabilities
#'
#'
#' #============================================================================
#' # Demonstrating the return of the prior when all data are missing
#' #============================================================================
#' # Data
#' y_missing <- list(enroll = list(m1 = NA),
#'                   recur1 = list(m1 = NA),
#'                   recur2 = list(m1 = NA))
#'
#' # Return of the prior
#' suppressMessages(compute_posterior(y_missing, fs = list(m1 = c("A" = 1))))
#'
#' # Return of the prior re-weighted to the exclusion of recrudescence
#' suppressMessages(compute_posterior(y_missing, fs = list(m1 = c("A" = 1)),
#'                  MOIs = c(1,2,3)))
#'
#' # (Recrudescing parasites are clones of previous blood-stage parasites. The
#' # Pv3R model assumes no within-host de-novo mutations and perfect allele
#' # detection. As such, recrudescence is incompatible with an MOI increase on
#' # the preceding infection.)
#'
#'
#' # ===========================================================================
#' # Demonstrating the cosmetic-only nature of episode names
#' # ===========================================================================
#' # Data
#' y <- list(enroll = list(m1 = NA),
#'           recur2 = list(m1 = NA),
#'           recur1 = list(m1 = NA))
#'
#' # Use a non-uniform prior for the purpose of illustration
#' prior <- matrix(c(0.2,0.2,0.6,0.7,0.1,0.2), byrow = TRUE, nrow = 2,
#'                 dimnames = list(c("recur1_prior", "recur2_prior"),
#'                                 c("C", "L", "I")))
#'
#' # Print posterior and prior, noting that "recur1_prior" is returned for
#' # "recur2", and "recur2_prior" is returned for "recur1"
#' suppressMessages(compute_posterior(y, fs = list(m1 = c(a = 1)), prior))$marg
#' prior
#'
#'
#' #============================================================================
#' # Demonstrating the informative nature of non-recurrent data
#' #============================================================================
#' # Data and allele frequencies
#' y_het <- list(list(m1 = c('1', '2')), list(m1 = NA))
#' y_hom <- list(list(m1 = '1'), list(m1 = NA))
#' fs = list(m1 = c('1' = 0.5, '2' = 0.5))
#'
#' # The prior is not returned despite there being no recurrent data (see
#' # vignette XXX to understand why)
#' suppressMessages(compute_posterior(y = y_het, fs))$marg
#' suppressMessages(compute_posterior(y = y_hom, fs, MOIs = c(2,1)))$marg
#'
#'
#'
#' #============================================================================
#' # Demonstrating the effect of increasingly large relationship graphs: the
#' # marginal probability of the first recurrence changes slightly, albeit at a
#' # decreasing rate, as the number of additional recurrences (all without data)
#' # increases. The change is greatest when the observed allele is rare.
#' #============================================================================
#' # Data for different recurrence counts where only the 1st recurrence has data
#' ys <- list(scenario1 = list(enroll = list(m1 = "A"),
#'                             recur1 = list(m1 = "A")),
#'            scenario2 = list(enroll = list(m1 = "A"),
#'                             recur1 = list(m1 = "A"),
#'                             recur2 = list(m1 = NA)),
#'            scenario3 = list(enroll = list(m1 = "A"),
#'                             recur1 = list(m1 = "A"),
#'                             recur2 = list(m1 = NA),
#'                             recur3 = list(m1 = NA)),
#'            scenario4 = list(enroll = list(m1 = "A"),
#'                             recur1 = list(m1 = "A"),
#'                             recur2 = list(m1 = NA),
#'                             recur3 = list(m1 = NA),
#'                             recur4 = list(m1 = NA)))
#'
#' # Allele frequencies: smaller f_A leads to larger change
#' f_A <- 0.1; fs <- list(m1 = c("A" = f_A, "Other" = 1-f_A))
#'
#' # Compute posterior probabilities and extract marginal probabilities
#' results <- lapply(ys, function(y) compute_posterior(y, fs)$marg)
#'
#' # Extract results for the first recurrence
#' results_recur1 <- sapply(results, function(result) result[1,])
#' results_recur1 # Results are different for different scenarios
#'
#' # Visualise the change in the marginal probability of the first recurrence
#' plot_simplex(c("Recrudescence", "Relapse", "Reinfection")) # Plot simplex
#' xy <- apply(results_recur1, 2, project2D) # Project probabilities
#' points(x = xy["x", ], y = xy["y", ], pch = "-", col = 1:4) # Plot projections
#' legend("left", col = 1:4, pch = "-", pt.cex = 2, bty = "n", legend = 1:4,
#' title = "Recurrence \n count") # legend
#'
#' @export
compute_posterior <- function(y, fs, prior = NULL, MOIs = NULL,
                              return.RG = FALSE, return.logp = FALSE) {

  # Check y is a list of lists:
  if (!is.list(y) | !all(sapply(y, is.list))) {
    stop("Data y must be a list of lists,
         even if only one marker is typed per infection")
  }

  # Check there are at least 2 episodes
  infection_count <- length(y)
  if (!infection_count > 1) {
    stop("Need more than 1 episode")
  }

  # Check frequencies sum to one:
  if (!all(abs(1 - round(sapply(fs, sum))) < .Machine$double.eps^0.5)) {
    stop('For a given marker, allele frequencies must sum to one')
  }

  # Check prior
  causes <- c("C", "L", "I")
  if (!is.null(prior)) {
    if (!identical(colnames(prior), causes)) { # Column names
      stop('The prior should have colnames "C", "L" and "I"')
    }
    if (!all(abs(1 - rowSums(prior)) < .Machine$double.eps^0.5)) { # Sums to one
      stop('Prior probabilities for a given recurrence must sum to one')
    }
  } else {
    prior <- matrix(rep(1 / 3, 3 * (infection_count - 1)), ncol = 3)
    colnames(prior) <- causes
  }

  y <- prep_data(y) # Process data: remove repeats and NAs
  ms <- unique(do.call(c, lapply(y, names))) # Extract marker names
  names_y <- names(y) # NULL if no episode names
  names_prior <- rownames(prior) # Null if prior = NULL

  # Check all episodes have the same marker names
  if (!all(all(sapply(y, function(y_m) identical(sort(names(y_m)), sort(ms)))))){
    stop("Markers are inconsistent across episodes.
    NB: If not all markers are typed per episode, data on untyped markers can be encoded as missing using NAs.")
  }

  # Check all markers have allele frequencies
  if (!all(ms %in% names(fs))) {
    stop("Some markers are missing allele frequencies")
  }

  # Check all alleles have a named frequency
  # For each marker, list allele names
  as <- lapply(ms, function(m) unique(as.vector(unlist(sapply(y, function(yt) yt[[m]])))))
  names(as) <- ms # Name list entries by marker
  # For each marker, check all alleles have a named frequency:
  have_freq <- all(sapply(ms, function(m) all(as[[m]][!is.na(as[[m]])] %in% names(fs[[m]]))))
  if (!have_freq) {
    stop("Not all alleles have a named frequency")
  }

  # Check external MOIs are compatible with data if provided
  min_MOIs <- determine_MOIs(y)
  if (is.null(MOIs)) {
    MOIs <- min_MOIs
  } else {
    if(!all(MOIs >= min_MOIs)) {
      stop("MOIs provided do not support the observed allelic diversity")
    }
  }

  # PROCEED TO WARNINGS

  # Warn if episode names are liable to cause confusion
  if (!is.null(names_y) &
      !is.null(names_prior) &
      !all(names_y[-1] == names_prior)) {
    warning("Data and prior episode names disagree\n")
  }

  # Warn if there is only allele information for < 2 episodes (unpaired)
  for (m in ms) {
    # boolean vector for whether each episode has allele information for marker m
    has.allele <- sapply(y, function(y.epi) any(!is.na(y.epi[[m]])))
    if(sum(has.allele) < 2) {
      warning(paste("Marker", m, "has data on one episode only"))
    }
  }

  # Warn users if any episodes have no data
  na_episodes <- sapply(y, function(x) all(is.na(x)))
  if (sum(na_episodes) == 1) {
    warning(sprintf("Episode %s has no data\n", names(which(na_episodes))))
  } else if (sum(na_episodes) > 1) {
    warning(sprintf("Episodes %s have no data\n", paste(names(which(na_episodes)), collapse = " & ")))
  }

  gs_count <- sum(MOIs)
  gs <- paste0("g", 1:gs_count)
  ts_per_gs <- rep(1:infection_count, MOIs)
  gs_per_ts <- split(gs, ts_per_gs)
  n_recur <- infection_count - 1

  # Setting return.logp=T will enforce return.RG=T
  if (return.logp & !return.RG) {
    return.RG <- T
    message("return.RG is overridden to TRUE since return.logp=TRUE")
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
        y[[t]], gs_per_ts[[t]], !return.logp
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
  n_rg_per_rstr <- stats::setNames(rep(0, n_rstrs), rstrs)
  # for each R, stores log sum of p(y|RG) across RGs compatible with R
  logp_sum_per_rstr <- stats::setNames(rep(-Inf, n_rstrs), rstrs)

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

  rownames(marg) <- names_y[-1] # drop first infection (not reinfection)

  result <- list(marg = marg, joint = post_per_rstr)
  if (return.RG) {
    if (!return.logp) {
      for(i in 1:n.RG) {
        RGs[[i]]$logp <- NA
      }
    }
    result$RGs <- RGs
  }
  result
}
