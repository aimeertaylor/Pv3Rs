#' Computes the posterior distribution of recurrence states
#'
#'
#' @param y Observed data in the form as a list. The number of entries is the
#'   number of episode. Each entry is in turn a list of observed alleles
#'   for each marker.
#' @param fs List of allele frequencies. Names of the list must match with the
#'   marker names in `y`.
#' @param prior Matrix of prior probabilities of the recurrence states for each
#'    recurrent episode. Each row corresponds to an episode, and the column names
#'    must be C, L, and I. If prior not provided, a uniform prior will be used.
#' @param return.RG Boolean for whether to return the relationship graphs, defaults
#'   to `FALSE`.
#'
#' @return List containing:
#'   \describe{
#'     \item{marg}{Matrix of marginal posterior probabilities of the possible
#'       recurrence states for each recurrent episode.}
#'     \item{joint}{Vector of joint posterior probabilities of each possible string
#'       of recurrence states.}
#'     \item{RGs}{List of relationship graphs with their log-likelihoods stored.
#'       Only returned if `return.RG` is `TRUE`.}
#'   }
#'
#' @details Assumptions:
#' \itemize{
#'   \item{No genotyping errors and missing data}
#'   \item{Parasites are outbred}
#'   \item{Relationship graphs are equally likely}
#'   \item{Uniform prior over recurrent states}
#' }
#'
#' @examples
#' y <- list(list(m1=c("A","C"), m2=c("G")), list(m1=c("A"), m2=c("G","T")))
#' fs <- list(m1=setNames(c(0.4, 0.6), c("A","C")),
#'            m2=setNames(c(0.2, 0.8), c("G","T")))
#' prior <- matrix(c(0.2, 0.3, 0.5), nrow=1)
#' post <- compute_posterior(y, fs, prior)
#' post$marg
#'
#' @export
compute_posterior <- function(y, fs, prior=NULL, return.RG=FALSE) {
  infection_count <- length(y)
  stopifnot(infection_count > 1)

  MOIs <- determine_MOIs(y)
  gs_count <- sum(MOIs)
  gs <- paste0("g", 1:gs_count)
  ts_per_gs <- rep(1:infection_count, MOIs)
  gs_per_ts <- split(gs, ts_per_gs)

  n_recur <- infection_count-1
  causes <- c("C","L","I")

  if(is.null(prior)) prior <- matrix(rep(1/3, 3*(infection_count-1)), ncol=3)
  if(is.null(colnames(prior))) colnames(prior) <- causes

  # for each infection we find the possible allele assignments for its genotypes
  # taking into account of permutation symmetry of the genotypes, we fix the
  # assignment of a marker's alleles, where the number of alleles for that marker
  # must equal to the MOI
  # each result is a list, with one dataframe for each marker
  # the rows of each dataframe correspond to a different assignment, whereas each
  # column corresponds to a genotype
  alleles_per_inf_per_m <- lapply(1:infection_count,
                                  function(t) enumerate_alleles(y[[t]],
                                                                gs_per_ts[[t]]))

  # take the cartesian product over infections to get allele assignments over all
  # genotypes, but still separated for each marker
  # we do not need to take the cartesian product over markers due to the outbred
  # assumption
  alleles_per_m <- Reduce(
    function(x, y) mapply(merge, x, y, MoreArgs=list(by=NULL, all=T), SIMPLIFY=F),
    alleles_per_inf_per_m
  )

  RGs <- RG_inference(MOIs, fs, alleles_per_m)

  # get all vectors of recurrence states
  n_rstrs <- 3^n_recur
  rstrs <- do.call(paste0, expand.grid(lapply(1:n_recur,
                                              function(x) causes)))
  n_rg_per_rstr <- setNames(rep(0, n_rstrs), rstrs)
  logp_sum_per_rstr <- setNames(rep(-Inf, n_rstrs), rstrs)

  # build up p(Y|R) = sum of p(y|RG) for all RG compatible with R
  RG_i <- 0
  n.RG <- length(RGs)
  pbar <- txtProgressBar(min=0, max=n.RG)
  writeLines("\nFinding log-likelihood of each vector of recurrent states")
  for(RG in RGs) {
    prob_RG <- exp(RG$logp) # p(y|RG)
    for(rstr in compatible_rstrs(RG, gs_per_ts)) {
      n_rg_per_rstr[rstr] <- n_rg_per_rstr[rstr] + 1
      logp_sum_per_rstr[rstr] <- log(exp(logp_sum_per_rstr[rstr])+prob_RG)
    }
    RG_i <- RG_i + 1
    setTxtProgressBar(pbar, RG_i)
  }
  writeLines("\n")

  # this step assumes that prior distribution of RGs is uniform, even if its
  # likelihood is zero
  logp_per_rstr <- logp_sum_per_rstr - log(n_rg_per_rstr)
  logp_per_rstr[is.nan(logp_per_rstr)] <- -Inf

  # add log prior to logp
  logprior <- log(prior)
  for(rstr in rstrs) {
    logprior_rstr <- sum(sapply(1:n_recur, function(i) logprior[i,substr(rstr,i,i)]))
    logp_per_rstr[rstr] <- logp_per_rstr[rstr] + logprior_rstr
  }
  post_per_rstr <- exp(logp_per_rstr) / exp(matrixStats::logSumExp(logp_per_rstr))

  marg <- matrix(NA, nrow=n_recur, ncol=3)
  colnames(marg) <- causes
  for(i in 1:n_recur) {
    idx_list <- split(1:n_rstrs, (1:n_rstrs-1) %/% 3^(i-1) %% 3 + 1)
    for(j in 1:3) {
      marg[i,j] <- sum(post_per_rstr[idx_list[[j]]])
    }
  }

  rownames(marg) <- names(y)[-1]

  result <- list(marg=marg, joint=post_per_rstr)
  if(return.RG) result$RGs <- RGs
  result
}

