################################################################################
# Each test should cover a single unit of functionality; see
# (https://r-pkgs.org/testing-basics.html).
################################################################################
set.seed(1)
causes <- c("C", "L", "I")
alleles <- c("A", "C", "G", "T")
n.as <- length(alleles)
alphas <- rep(1, n.as) # Dirichlet param. vector
n.m <- 3
fmat <- gtools::rdirichlet(n.m, alphas)
colnames(fmat) <- alleles
fs <- setNames(lapply(1:n.m, function(x) fmat[x,]), paste0("m", 1:n.m))

testthat::test_that("No data returns the prior", {

  # Default prior, one recurrence
  y <- list(list(m1 = NA), list(m1 = NA))
  expect <- matrix(rep(1/3,3), ncol = 3, dimnames = list(NULL, causes))
  totest <- suppressMessages(compute_posterior(y, fs)$marg)
  testthat::expect_equal(totest, expect)

  # Default prior, multiple recurrence
  y <- list(list(m1 = NA), list(m1 = NA), list(m1 = NA))
  expect <- matrix(rep(1/3,6), nrow = 2, byrow = T, dimnames = list(NULL, causes))
  totest <- suppressMessages(compute_posterior(y, fs)$marg)
  testthat::expect_equal(totest, expect)

  # User-specified prior
  y <- list(list(m1 = NA), list(m1 = NA))
  prior <- gtools::rdirichlet(1, rep(1,3))
  colnames(prior) <- causes
  expect <- prior
  totest <- suppressMessages(compute_posterior(y, fs, prior)$marg)
  testthat::expect_equal(totest, expect)
})


# Verifies the examples presented in https://doi.org/10.1101/2022.11.23.22282669
testthat::test_that("Examples from the pre-print agree", {

  # Example 1
  y <- list(list(m1="A"), list(m1="T"))
  post <- suppressMessages(compute_posterior(y, fs))
  expect <- matrix(c(0, 1/3, 2/3), ncol=3, byrow=T, dimnames=list(NULL, causes))
  testthat::expect_equal(post$marg, expect)

  # Example 2
  y <- list(list(m1="A"), list(m1="A"))
  post <- suppressMessages(compute_posterior(y, fs))
  l.ntor <- fs$m1['A']+1
  i.ntor <- 2*fs$m1['A']
  c.ntor <- 2
  dtor <- sum(c(l.ntor, i.ntor, c.ntor))
  expect <- matrix(c(c.ntor, l.ntor, i.ntor)/dtor,
                   ncol=3, byrow=T, dimnames=list(NULL, causes))
  testthat::expect_equal(post$marg, expect)

  # Example 3
  y <- list(list(m1="A", m2="T", m3="G"),
            list(m1="T", m2="T", m3="G"))
  post <- suppressMessages(compute_posterior(y, fs))
  l.ntor <- (fs$m2['T']*fs$m3['G']+(fs$m2['T']+1)*(fs$m3['G']+1)/8)/3
  i.ntor <- fs$m2['T']*fs$m3['G']
  c.ntor <- 0
  dtor <- sum(c(l.ntor, i.ntor, c.ntor))
  expect <- matrix(c(c.ntor, l.ntor, i.ntor)/dtor,
                   ncol=3, byrow=T, dimnames=list(NULL, causes))
  testthat::expect_equal(post$marg, expect)


  # Example 4
  y <- list(list(m1="A", m2="T", m3="G"),
            list(m1="A", m2="T", m3="G"))
  post <- suppressMessages(compute_posterior(y, fs))
  l.ntor <- (fs$m1['A']+1)*(fs$m2['T']+1)*(fs$m3['G']+1)/8 + fs$m1['A']*fs$m2['T']*fs$m3['G'] + 1
  i.ntor <- 3*fs$m1['A']*fs$m2['T']*fs$m3['G']
  c.ntor <- 3
  dtor <- sum(c(l.ntor, i.ntor, c.ntor))
  expect <- matrix(c(c.ntor, l.ntor, i.ntor)/dtor,
                   ncol=3, byrow=T, dimnames=list(NULL, causes))
  testthat::expect_equal(post$marg, expect)


  # Example 5
  y <- list(list(m1="A"), list(m1="T"), list(m1="T"))
  post <- suppressMessages(compute_posterior(y, fs))
  l1.ntor <- 17/24*fs$m1['A']*fs$m1['T']^2 + 11/16*fs$m1['A']*fs$m1['T']
  i1.ntor <- 7/5*fs$m1['A']*fs$m1['T']^2 + 13/10*fs$m1['A']*fs$m1['T']
  c1.ntor <- 0
  dtor1 <- sum(c(l1.ntor, i1.ntor, c1.ntor))
  l2.ntor <- 73/120*fs$m1['A']*fs$m1['T']^2 + 39/80*fs$m1['A']*fs$m1['T']
  i2.ntor <- 3/2*fs$m1['A']*fs$m1['T']^2
  c2.ntor <- 3/2*fs$m1['A']*fs$m1['T']
  dtor2 <- sum(c(l2.ntor, i2.ntor, c2.ntor))
  expect <- matrix(c(c(c1.ntor, l1.ntor, i1.ntor)/dtor1,
                     c(c2.ntor, l2.ntor, i2.ntor)/dtor2),
                   ncol=3, byrow=T, dimnames=list(NULL, causes))
  testthat::expect_equal(post$marg, expect)


  # Example 6 (hand computed likelihood / Pv3R computed posterior)
  y <- list(list(m1="A"), list(m1="A"), list(m1="A"))
  post <- suppressMessages(compute_posterior(y, fs))
  ll_liklihood <- fs$m1['A']^3
  cc_liklihood <- fs$m1['A']
  testthat::expect_equal(ll_liklihood/post$joint["II"],
                         cc_liklihood/post$joint["CC"])


  # Example 7
  y <- list(list(m1=c("A","T")), list(m1="T"))
  post <- suppressMessages(compute_posterior(y, fs))
  l.ntor <- 5/18*fs$m1['A']*fs$m1['T']^2 + 1/4*fs$m1['A']*fs$m1['T']
  i.ntor <- 3/4*fs$m1['A']*fs$m1['T']^2
  c.ntor <- 3/8*fs$m1['A']*fs$m1['T']
  dtor <- sum(c(l.ntor, i.ntor, c.ntor))
  expect <- matrix(c(c.ntor, l.ntor, i.ntor)/dtor,
                   ncol=3, byrow=T, dimnames=list(NULL, causes))
  testthat::expect_equal(post$marg, expect)

  # Example 9
  y <- list(list(m1 = c("A","T"), m2 = "T", m3 = c("C", "G")),
            list(m1 = "T", m2 = "T", m3 = "C"))
  post <- suppressMessages(compute_posterior(y, fs))
  l_ntor <- (1/9)*  (2* fs$m1["T"]*fs$m2["T"]^2*fs$m3["C"] +
                       (3/8)* fs$m1["T"]*(fs$m2["T"]^2 + fs$m2["T"])*fs$m3["C"] +
                       (1/8)* fs$m1["T"]*(fs$m2["T"]^2 + fs$m2["T"])*(fs$m3["C"] + 1) +
                       (1/8)*(fs$m1["T"] + 1)*(fs$m2["T"]^2 + fs$m2["T"])*fs$m3["C"] +
                       (1/8)*(fs$m1["T"] + 1)*(fs$m2["T"]^2 + fs$m2["T"])*(fs$m3["C"] + 1) +
                       (1/32)*(3*fs$m2["T"] + 1) +
                       (1/8)*(9*fs$m2["T"] + 1))
  i_ntor <- (1/8)*fs$m1["T"]*fs$m2["T"]*fs$m3["C"]*(9*fs$m2["T"] + 1)
  c_ntor <- (1/32)*(9*fs$m2["T"] + 1)
  dtor <- sum(c(l_ntor, i_ntor, c_ntor))
  expect <- matrix(c(c_ntor, l_ntor, i_ntor)/dtor,
                   ncol=3, byrow=T, dimnames=list(NULL, causes))
  testthat::expect_equal(post$marg, expect)

})


# Check for under/overflow issues
testthat::test_that("Check that 1000 markers does not lead to under/overflow", {
  m <- 1000
  markers <- paste0("m", 1:m) # Marker names
  alleles <- letters # Alleles
  n_alleles <- length(alleles) # Number of alleles per marker

  # Sample allele frequencies
  fs <- sapply(markers, function(m) {
      fs_unnamed <- gtools::rdirichlet(1, alpha = rep(1, n_alleles))
      setNames(fs_unnamed, alleles)
    }, USE.NAMES = TRUE, simplify = FALSE
  )

  # Sample parental genotypes
  parent1 <- sapply(markers, function(t) {
    sample(alleles, size = 1, prob = 1-fs[[t]])}, simplify = F)
  parent2 <- sapply(markers, function(t) {
    sample(alleles, size = 1, prob = 1-fs[[t]])}, simplify = F)

  # Sample children genotypes (ensure all different)
  anyclones <- TRUE
  while (anyclones) {
    child1 <-  sapply(markers, function(t) sample(c(parent1[[t]], parent2[[t]]), 1), simplify = F)
    child2 <-  sapply(markers, function(t) sample(c(parent1[[t]], parent2[[t]]), 1), simplify = F)
    child3 <-  sapply(markers, function(t) sample(c(parent1[[t]], parent2[[t]]), 1), simplify = F)
    anyclones <- any(identical(child1, child2),
                     identical(child2, child3),
                     identical(child1, child3))
  }

  # Make parasite infection (incompatible with recrudescence)
  initial <- rbind(unlist(child1))
  relapse <- rbind(unlist(child2), unlist(child3))
  y <- list(apply(initial, 2, unique, simplify = F),
            apply(relapse, 2, unique, simplify = F))

  # Compute posterior
  post <- suppressMessages(compute_posterior(y, fs))

  expect <- matrix(c(0, 1, 0), ncol=3, byrow=T, dimnames=list(NULL, causes))
  testthat::expect_equal(post$marg, expect)
})


# Likelihood involving missing data should be equivalent to summing likelihoods
# over all possible imputations of missing data
testthat::test_that("Check NAs are handled correctly if no data for first recurrence", {
  alleles <- c("A", "B", "C")
  get_logp <- function(post) {
    sapply(post$RGs, function(RG) RG$logp)
  }
  fs <- list(m1=setNames(gtools::rdirichlet(1, rep(1, length(alleles))), alleles))

  # simulate data for initial episode, MOI = 2
  a1 <- unique(sample(alleles, 2, replace=T, prob=fs$m1))
  # data is NA for first recurrence, MOI = 2
  # simulate data for second recurrence, MOI = 1
  a3 <- sample(alleles, 1, replace=T, prob=fs$m1)

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
  testthat::expect_equal(logp_NA, apply(logp_mat, 1, matrixStats::logSumExp))
})
