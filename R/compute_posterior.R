#' Compute posterior probabilities of \emph{P. vivax} recurrence states
#'
#' @description
#' Computes posterior probabilities of \emph{P. vivax} recurrence states
#' relapse, reinfection and recrudescence using genetic data; for usage see
#' **Examples** below and
#' [vignette("demonstrate-usage")](../doc/demonstrate\-usage.html), for a more
#' complete understanding of the posterior output see
#' ["Understand posterior estimates"](https://aimeertaylor.github.io/Pv3Rs/articles/understand-posterior.html).
#'
#' Note: The progress bar increments non-uniformly (see
#' **Details**); it may appear stuck when computations are ongoing.
#'
#' @details
#' `compute_posterior()` computes posterior probabilities proportional to
#' the likelihood multiplied by the prior. The likelihood sums over:
#'
#' - ways to phase allelic data onto haploid genotypes
#' - graphs of relationships between haploid genotypes
#' - ways to partition alleles into clusters of identity-by-descent
#'
#' \code{compute_posterior()} expects each per-episode, per-marker allelic
#' vector to contain distinct alleles only. Allele repeats at markers with
#' observed data, and \code{NA} repeats at markers with missing data, are
#' removed in a data pre-processing step. \code{NA}s in allelic vectors that
#' also contain non-\code{NA} values are also removed in the data pre-processing
#' step.
#'
#' We enumerate all possible relationship graphs between haploid genotypes,
#' where pairs of genotypes can either be clones, siblings, or strangers. The
#' likelihood of a sequence of recurrence states can be determined from the
#' likelihood of all relationship graphs compatible with said sequence. More
#' details on the enumeration of relationship graphs can be found in
#' \code{\link{enumerate_RGs}}. For each relationship graph, the model sums over
#' all possible identity-by-descent partitions. Because some graphs
#' are compatible with more partitions than others, the log
#' p(Y|RG) progress bar may advance non-uniformly.
#'
#' We do not recommend running `compute_posterior() when the total genotype
#' count (sum of MOIs) exceeds eight because there are too many relationship
#' graphs.
#'
#' Notable model assumptions and limitations:
#' \itemize{
#'   \item{All siblings are regular siblings}
#'   \item{Recrudescent parasites derive only from the immediately preceding
#'   episode}
#'   \item{Recrudescence, relapse and reinfection are mutually exclusive}
#'   \item{Undetected alleles, genotyping errors, and \emph{de novo} mutations
#'   are not modelled}
#'   \item{Population structure and various other complexities that confound
#'   molecular correction are not modelled}
#' }
#'
#'
#' @param y List of lists encoding allelic data. The outer list contains
#'   episodes in increasing chronological order. The inner list contains named
#'   markers per episode. Marker names must be consistent across episodes. `NA`
#'   indicates missing marker data; otherwise, one must specify a per-marker
#'   vector of distinct alleles detected (presently, `compute_posterior` does
#'   not support data on the proportional abundance of detected alleles).
#'   \code{NA} entries within allelic vectors are ignored. Allele names are
#'   arbitrary, allowing for different data types, but must correspond with
#'   frequency names.
#' @param fs List of per-marker allele frequency vectors, with names matching
#'   marker names in `y`. Per-marker alleles frequencies mut contain one
#'   frequency per named allele, with names matching alleles in `y`. Per-marker
#'   frequencies must sum to one.
#' @param prior Matrix of prior probabilities of recurrence states per episode,
#'   with rows as episodes in chronological order, and columns named "C", "L",
#'   and "I" for recrudescence, relapse and reinfection, respectively. Row names
#'   are ignored. If `NULL` (default), per-episode recurrence states are assumed
#'   equally likely \emph{a priori}.
#' @param MOIs Vector of multiplicity of infection (MOI) per
#'   episode. If `NULL` (default), the most parsimonious MOIs compatible
#'   with the data are used; see \code{\link{determine_MOIs}}.
#' @param return.RG Logical; whether to return the relationship graphs,
#'   (default `FALSE`). Automatically set to `TRUE` if `return.logp = TRUE`.
#' @param return.logp Logical; whether to return the log-likelihood for each
#'   relationship graph (default `FALSE`). Setting `TRUE` disables permutation
#'   symmetry optimisation; increasing runtime, especially when MOIs are large.
#'   Does not affect the output of the posterior probabilities; see
#'   `vignette("demonstrate-usage")`, which also example of permutation
#'   symmetry.
#' @param progress.bar Logical; whether to print progress bars (default `TRUE`).
#'   Note that the progress bar may update non-uniformly and appear when the
#'   code is still running.
#'
#'
#' @return List containing:
#'   \describe{
#'     \item{`marg`}{Matrix of marginal posterior probabilities for each recurrence,
#'       with rows as recurrences and columns as "C" (recrudescence), "L"
#'       (relapse), and "I" (reinfection). Each marginal probability sums over a
#'       subset of joint probabilities. For example, the marginal probability of
#'       "C" at the first of two recurrences sums over the joint probabilities
#'       of "CC", "CL", and "CI".}
#'     \item{`joint`}{Vector of joint posterior probabilities for each recurrence
#'       state sequence; within a sequence "C", "L", and "I" are used as above}
#'     \item{`RGs`}{List of relationship graphs returned if
#'       `return.RG = TRUE`, with log-likelihoods if `return.logp = TRUE`;
#'       for more details on relationship graphs, see \code{\link{enumerate_RGs}}.}
#'   }
#'
#' @examples
#' # Numerically named alleles
#' y <- list(enrol = list(m1 = c('3','2'), m2 = c('1','2')),
#'           recur1 = list(m1 = c('1','4'), m2 = c('1')),
#'           recur2 = list(m1 = c('1'), m2 = NA))
#' fs <- list(m1 = c('1' = 0.78, '2' = 0.14, '3' = 0.07, '4' = 0.01),
#'            m2 = c('1' = 0.27, '2' = 0.73))
#' compute_posterior(y, fs, progress.bar = FALSE)
#'
#'
#' # Arbitrarily named alleles, plotting per-recurrence posteriors
#' y <- list(enrolment = list(marker1 = c("Tinky Winky", "Dipsy"),
#'                           marker2 = c("Tinky Winky", "Laa-Laa", "Po")),
#'           recurrence = list(marker1 = "Tinky Winky",
#'                           marker2 = "Laa-Laa"))
#' fs <- list(marker1 = c("Tinky Winky" = 0.4, "Dipsy" = 0.6),
#'            marker2 = c("Tinky Winky" = 0.1, "Laa-Laa" = 0.1, "Po" = 0.8))
#' plot_simplex(p.coords = compute_posterior(y, fs, progress.bar = FALSE)$marg)
#'
#'
#' # Episode names are cosmetic: "r1_prior" is returned for "r2"
#' y <- list(enrol = list(m1 = NA), r2 = list(m1 = NA), r1 = list(m1 = NA))
#' prior <- matrix(c(0.6,0.7,0.2,0.3,0.2,0), ncol = 3,
#'                 dimnames = list(c("r1_prior", "r2_prior"), c("C", "L", "I")))
#' suppressMessages(compute_posterior(y, fs = list(m1 = c(a = 1)), prior))$marg
#' prior
#'
#'
#' # Prior is returned when all data are missing
#' y_missing <- list(enrol = list(m1 = NA), recur = list(m1 = NA))
#' suppressMessages(compute_posterior(y_missing, fs = list(m1 = c("A" = 1))))
#'
#'
#' # Return of the prior re-weighted to the exclusion of recrudescence:
#' suppressMessages(compute_posterior(y_missing, fs = list(m1 = c("A" = 1)),
#'                  MOIs = c(1,2)))
#' # (Recrudescing parasites are clones of previous blood-stage parasites. The
#' # Pv3R model assumes no within-host de-novo mutations and perfect allele
#' # detection. As such, recrudescence is incompatible with an MOI increase on
#' # the preceding infection.)
#'
#'
#' # Beware provision of unpaired data: the prior is not necessarily returned;
#' # for more details, see link above to "Understand posterior estimates"
#' y_het <- list(list(m1 = c('1', '2')), list(m1 = NA))
#' y_hom <- list(list(m1 = '1'), list(m1 = NA))
#' fs = list(m1 = c('1' = 0.5, '2' = 0.5))
#' suppressMessages(compute_posterior(y = y_het, fs))$marg
#' suppressMessages(compute_posterior(y = y_hom, fs, MOIs = c(2,1)))$marg
#'
#' @export
compute_posterior <- function(y, fs, prior = NULL, MOIs = NULL,
                              return.RG = FALSE, return.logp = FALSE,
                              progress.bar = TRUE) {

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
  if (!all(abs(1 - sapply(fs, sum)) < .Machine$double.eps^0.5)) {
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
  na_markers <- sapply(ms, function(m)
    sum(sapply(y, function(y.epi) any(!is.na(y.epi[[m]])))) < 2
  )
  if (sum(na_markers) == 1) {
    warning(sprintf("Markers %s has data on one episode only\n",
                    names(which(na_markers))))
  } else if (sum(na_markers) > 1) {
    warning(sprintf("Markers %s have data on one episode only\n",
                    paste(names(which(na_markers)), collapse = " & ")))
  }

  # NB: The previous warning does necessarily trigger when the next warning does,
  #     because there can be a total of >2 episodes, with one episode without
  #     data, but every marker still has data on multiple episodes.

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
  RGs <- RG_inference(MOIs, fs, alleles_per_m, progress.bar)

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
  if (progress.bar) pbar <- msg_progress_bar(n.RG)
  message("Finding log-likelihood of each vector of recurrence states")
  max.logp <- max(sapply(RGs, "[[", "logp"))
  for (RG in RGs) {
    # subtract maximum to avoid underflow
    prob_RG <- exp(RG$logp - max.logp) # p(y|RG)
    for (rstr in compatible_rstrs(RG, gs_per_ts)) {
      n_rg_per_rstr[rstr] <- n_rg_per_rstr[rstr] + 1
      logp_sum_per_rstr[rstr] <- log(exp(logp_sum_per_rstr[rstr]) + prob_RG)
    }
    RG_i <- RG_i + 1
    if (progress.bar) pbar$increment()
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
