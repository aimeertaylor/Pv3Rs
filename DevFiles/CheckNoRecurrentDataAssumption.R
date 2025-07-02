################################################################################
# Script to empirically validate eq 4 of graph-prior.ltx (tex file)
################################################################################
devtools::load_all()

# Markers, Allele frequencies, Data
MOIs <- c(1,1,1)
ms <- paste0("m", 1:100)
fA <- 0.15
fs <- sapply(ms, function(x) c("A" = fA, "B" = 1-fA), simplify = FALSE)
As <- sapply(ms, function(x) "A", simplify = FALSE)
Bs <- sapply(ms, function(x) "B", simplify = FALSE)
NAs <- sapply(ms, function(x) NA, simplify = FALSE)
y <- list(enrol = As, recur1 = As, recur2 = NAs)

# Generate gs_per_ts for compatible_rstrs
gs <- paste0("g", 1:sum(MOIs)) # genotype names (graph vertices)
ts <- 1:length(MOIs) # episode indices
ts_per_gs <- rep(ts, MOIs) # episode index of each genotype
gs_per_ts <- split(gs, ts_per_gs) # genotypes grouped by episode

result <- compute_posterior(y, fs, MOIs = MOIs, return.RG = TRUE, return.logp = TRUE)
CIL_gvn_RGs <- sapply(result$RGs, Pv3Rs:::compatible_rstrs, gs_per_ts) # compatible states
CC_CI_CL_log <- t(sapply(CIL_gvn_RGs, function(x) c("CC", "CI", "CL") %in% x))
RGlogps <- sapply(result$RGs, function(x) x$logp)
unique(RGlogps)

# Likelihood sums over graph space are equal up to a normalising constant when
# all episodes are monoclonal
c(SumRGs_CC = sum(exp(RGlogps[CC_CI_CL_log[,1]])),
  SumRGs_CI = sum(exp(RGlogps[CC_CI_CL_log[,2]])),
  SumRGs_CL = sum(exp(RGlogps[CC_CI_CL_log[,3]])))

# This is because all log likelihoods consistent with the first recurrence are the same:
RGlogps[CC_CI_CL_log[,1]]
RGlogps[CC_CI_CL_log[,2]]
RGlogps[CC_CI_CL_log[,3]]

# Check we can discount summation over graph spaces inconsistent with first recrudescence
sum(exp(RGlogps[CC_CI_CL_log[,1]]))
sum(exp(RGlogps[sapply(CIL_gvn_RGs, function(x) ("LC" %in% x) & (!"CC" %in% x))]))
sum(exp(RGlogps[sapply(CIL_gvn_RGs, function(x) "IC" %in% x)]))

sum(exp(RGlogps[CC_CI_CL_log[,2]]))
sum(exp(RGlogps[sapply(CIL_gvn_RGs, function(x) ("LI" %in% x) & (!"CI" %in% x))]))
sum(exp(RGlogps[sapply(CIL_gvn_RGs, function(x) "II" %in% x)]))

sum(exp(RGlogps[CC_CI_CL_log[,3]]))
sum(exp(RGlogps[sapply(CIL_gvn_RGs, function(x) ("LL" %in% x) & (!"CL" %in% x))]))
sum(exp(RGlogps[sapply(CIL_gvn_RGs, function(x) "IL" %in% x)]))


# Equation 4:
result$joint["CC"]
1/(3 + 2/3 + 3/12)

result$joint # Posteriors on CC, CI, and CL are equal...

# As such the marginal probability is three times that of CC
result$marg["recur1","C"]
3/(3 + 2/3 + 3/12)
