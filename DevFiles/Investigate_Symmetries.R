################################################################################
# Investigate equivalence classes
################################################################################

library(Pv3Rs)

y <- list(enroll = list(m1=c("A","B")), recur = list(m1=c("A")))
MOIs <- c(3,1) # suppose user specifies greater MOI than observed diversity
alleles <- c("A","B","C")
gs_per_ts <- split(paste0("g", 1:sum(MOIs)), rep(1:length(MOIs), MOIs))
gs <- paste0("g", 1:sum(MOIs))
ts_per_gs <- rep(1:length(y), MOIs)

# Example 1
fs <- list(m1=c(A=0.6, B=0.1, C=0.3))
res <- suppressMessages(compute_posterior(y, fs, MOIs=MOIs, return.logp=T))
res$joint # C < L < I
RG_mode <- res$RGs[[which.max(sapply(res$RGs, function(RG) RG$logp))]]
compatible_rstrs(RG_mode, gs_per_ts) # not compatible with I
par(mar = rep(0.1,4))
plot_RG(RG_to_igraph(RG_mode, gs, ts_per_gs), edge.curved=0.25, vertex.size=20)
# this one is hard to interpret because posterior probabilities are close

# Example 2
fs <- list(m1=c(A=0.04, B=0.12, C=0.84))
res <- suppressMessages(compute_posterior(y, fs, MOIs=MOIs, return.logp=T))
res$joint # C > L > I
RG_mode <- res$RGs[[which.max(sapply(res$RGs, function(RG) RG$logp))]]
compatible_rstrs(RG_mode, gs_per_ts) # not compatible with C
par(mar = rep(0.1,4))
plot_RG(RG_to_igraph(RG_mode, gs, ts_per_gs), edge.curved=0.25, vertex.size=20)
# most likely RG has all siblings edges, only graph in its equivalence class
# compare this with the RG where g1 and g2 are siblings, g4 is a clone of g1
# this RG is 'intuitively' the most likely when all permutations are considered (there are 6)



# The phenomenon that the most likely RG is not necessarily compatible with the
# most likely recurrence sequence can also occur when MOIs match observed diversity
y <- list(enroll = list(m1=c("A","B")), recur1 = list(m1=c("A")), recur2 = list(m1=c("B")))
MOIs <- c(2,1,1)
alleles <- c("A","B","C")
gs_per_ts <- split(paste0("g", 1:sum(MOIs)), rep(1:length(MOIs), MOIs))
# Fix data, alter allele frequencies f
# Print seeds whose corresponding f results in the phenomenon above
# Have not analysed the reasons behind the incompatibility
for(s in 1:50) {
  set.seed(s)
  fs <- list(m1=setNames(gtools::rdirichlet(1, rep(0.5, length(alleles)))[1,], alleles))
  res <- suppressMessages(compute_posterior(y, fs, MOIs=MOIs, return.logp=T))
  post_r <- names(which.max(res$joint))
  RG_mode <- res$RGs[[which.max(sapply(res$RGs, function(RG) RG$logp))]]
  if(post_r %in% compatible_rstrs(RG_mode, gs_per_ts) == F) print(s)
}
