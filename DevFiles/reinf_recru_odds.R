####################
# Investigate an example where the data the is compatible recrudescence, but the
# genetic evidence for recrudescence is weak. The evidence is weak in the sense
# that the recurrent allele has a high frequency.
#
# Reinfection gets a higher posterior probability, but this can be explained by
# the fact that the probabilistic distribution of alleles for the recurrent
# episode under recrudescence is uniform, instead of the `natural` allele
# frequencies, which is heavily biased towards the recurrent allele.
####################

library(Pv3Rs)

# returns a list of vectors of numbers, corresponding to relationship graphs
# that are compatible with recrudescence / relapse / refinection
list.compat <- function(MOI1, MOI2) {
  MOIs <- c(MOI1, MOI2)
  RGs <- suppressMessages(enumerate_RGs(MOIs=MOIs))
  gs_per_ts <- split(paste0("g", 1:sum(MOIs)), rep(1:length(MOIs), MOIs))
  compat <- lapply(RGs, Pv3Rs:::compatible_rstrs,  gs_per_ts)
  return(list(C=sapply(compat, function(rs) "C" %in% rs),
              L=sapply(compat, function(rs) "L" %in% rs),
              I=sapply(compat, function(rs) "I" %in% rs)))
}

# initial episode: common and rare alleles
# recurrence episode: common allele only
fs <- list(m1=c(A=0.99,B=0.01))
y <- list(enroll = list(m1=c("A","B")), recur = list(m1=c("A")))
MOIs <- c(2,1)

res <- suppressMessages(compute_posterior(y, fs, MOIs=MOIs, return.logp=T))

res$marg
res$marg[,"I"] / res$marg[,"C"] # posterior I:C odds close to 2

compat.lst <- list.compat(MOIs[1], MOIs[2])
likes <- sapply(res$RGs, function(RG) exp(RG$logp))
# ratio of highest graph likelihood compatible with I to highest graph likelihood compatible with C
max(likes[compat.lst$I]) / max(likes[compat.lst$C])


# repeat above for different MOIs
for (MOIs in list(c(2,1), c(3,1), c(3,2), c(4,1), c(4,2), c(4,3))) {
  cat("MOIs =", MOIs, "\n")
  cat("prediction based on MOIs", MOIs[1]/(MOIs[1]-MOIs[2]), "\n")

  res <- suppressMessages(compute_posterior(y, fs, MOIs=MOIs, return.logp=T))
  cat("posterior I:C odds", res$marg[,"I"] / res$marg[,"C"], "\n")

  compat.lst <- list.compat(MOIs[1], MOIs[2])
  likes <- sapply(res$RGs, function(RG) exp(RG$logp))
  cat("max p(y|g compat with I) / max p(y|g compat with C)",
      max(likes[compat.lst$I]) / max(likes[compat.lst$C]), "\n")
}

# the MOI-based prediction explains the ratio max p(y|g compat with I) / max p(y|g compat with C)
# the posterior odds is influenced also by relationship graphs with sibling edges
