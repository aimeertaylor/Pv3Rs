### Validate that NAs are handled correctly
#
# Example with MOIs = c(2, 2, 1), where data is missing for the first recurrence
# For each RG, check likelihood for data with NA = sum of likelihoods over all
# possible imputations for NA
#
# Takes about 2 minutes to run

library(Pv3Rs)
library(matrixStats)
library(MCMCpack)

alleles <- c("A", "B", "C")

get_logp <- function(post) {
  sapply(post$RGs, function(RG) RG$logp)
}

for(i in 1:10) {
  set.seed(i)
  fs <- list(m1=setNames(rdirichlet(1, rep(1, 3)), alleles)) # random frequencies

  # simulate data for initial episode
  a1 <- unique(sample(alleles, 2, replace=T, prob=fs$m1))
  # simulate data for second recurrence
  a3 <- sample(alleles, 1, replace=T, prob=fs$m1)
  # data with NA
  y <- list(initial=list(m1=a1), recur1=list(m1=NA), recur2=list(m1=a3))
  # likelihood calculation with NA
  post <- suppressMessages(compute_posterior(y, fs,
                                             return.RG=T, return.logp=T,
                                             MOIs=c(2,2,1)))
  logp_NA <- get_logp(post)

  # all possible imputations for first recurrence
  a2_list <- list("A", "B", "C", c("A","B"), c("B","C"), c("A","C"))
  # logp calculation for each imputation, each RG
  logp_mat <- sapply(a2_list, function(a2) {
    y_impute <- list(initial=list(m1=a1), recur1=list(m1=a2), recur2=list(m1=a3))
    post_impute <- suppressMessages(compute_posterior(y_impute, fs,
                                                      return.RG=T, return.logp=T,
                                                      MOIs=c(2,2,1)))
    get_logp(post_impute)
  })

  # check if logp with NA == logsumexp of logp with imputed values
  res <- all.equal(logp_NA, apply(logp_mat, 1, logSumExp))
  print(paste("Seed", i, ":", res))
}

