# Example of performing inference on recurrent states
# Inputs:
#   a) Observed data y in the form of a list, with one entry per infection.
#      Each infection is represent as a list, containing vectors of observed
#      alleles for every marker.
#   b) Allele frequencies in the form of a list, containing vectors of allele
#      frequencies for each marker. Allele names are found in the vector names.
# Output: Posterior probability of each vector of recurrent states.
#
# Assumptions:
#   - No genotyping errors and missing data
#   - Parasites are outbred
#   - Relationship graphs are equally likely
#   - Uniform prior over recurrent states

library(Pv3Rs)
library(partitions)
library(tictoc)

y <- readRDS("data/y_ICL_3212.rds") # observed data
fs <- readRDS("data/fs_ICL_3212.rds") # allele frequencies
ms <- rownames(fs) # marker names

infection_count <- length(y)
MOIs <- determine_MOIs(y)
gs_count <- sum(MOIs)
gs <- paste0("g", 1:gs_count)
ts_per_gs <- rep(1:infection_count, MOIs)
gs_per_ts <- split(gs, ts_per_gs)

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

# determine which allele assignments are viable given that two genotypes share
# the same allele
als_by_edge <- lapply(alleles_per_m, allele_filter)

tic("Compute log-likelihood of each RG")
RGs <- RG_inference(MOIs, fs, alleles_per_m, als_by_edge)
toc()

# get all vectors of recurrence states
n_recur <- infection_count-1
causes <- c("C","L","I")
n_rstrs <- 3^n_recur
rstrs <- do.call(paste0, expand.grid(lapply(1:n_recur,
                                            function(x) causes)))
n_rg_per_rstr <- setNames(rep(0, n_rstrs), rstrs)
logp_sum_per_rstr <- setNames(rep(-Inf, n_rstrs), rstrs)

# build up p(Y|R)
tic("Find log-likelihood of each vector of recurrent states")
for(RG in RGs) {
  prob_RG <- exp(RG$logp)
  for(rstr in compatible_rstrs(RG, gs_per_ts)) {
    n_rg_per_rstr[rstr] <- n_rg_per_rstr[rstr] + 1
    logp_sum_per_rstr[rstr] <- log(exp(logp_sum_per_rstr[rstr])+prob_RG)
  }
}
toc()

# this step assumes that prior distribution of RGs is uniform, even if its
# likelihood is zero
logp_per_rstr <- logp_sum_per_rstr - log(n_rg_per_rstr)
logp_per_rstr[is.nan(logp_per_rstr)] <- -Inf
# assume that prior of R is also uniform for now
post_per_rstr <- exp(logp_per_rstr) / exp(matrixStats::logSumExp(logp_per_rstr))
writeLines("\nPosterior probability of recurrent states")
print(post_per_rstr)

for(i in 1:n_recur) {
  writeLines(paste("Recurrence", i))
  idx_list <- split(1:n_rstrs, (1:n_rstrs-1) %/% 3^(i-1) %% 3 + 1)
  for(j in 1:3) {
    writeLines(paste0("Pr(R", i, "=", causes[j], "|Y)=",
                      sum(post_per_rstr[idx_list[[j]]])))
  }
}

# plot the RG that generated the observed data
par(mfrow=c(1,1))
RG0 <- readRDS("data/RG_ICL_3212.rds")
plot_RG(RG_to_igraph(RG0, gs, ts_per_gs), edge.curved=0.2)

# plot the top 9 RGs in terms of likelihood
par(mfrow=c(3,3), mar=c(1,1,1,1))
RG_order <- order(sapply(RGs, function(RG) RG$logp), decreasing=T)
for(RG_i in RG_order[1:9]) {
  plot_RG(RG_to_igraph(RGs[[RG_i]], gs, ts_per_gs), edge.curved=0.2)
}


