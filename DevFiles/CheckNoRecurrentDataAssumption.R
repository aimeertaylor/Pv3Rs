
devtools::load_all()

# Markers, Allele frequencies, Data
MOIs <- c(2,1,1)
ms <- paste0("m", 1:100)
fs <- sapply(ms, function(x) c("A" = 0.01, "B" = 0.99), simplify = FALSE)
As <- sapply(ms, function(x) "A", simplify = FALSE)
Bs <- sapply(ms, function(x) "B", simplify = FALSE)
NAs <- sapply(ms, function(x) NA, simplify = FALSE)
y <- list(enrol = As, recur1 = As, recur2 = NAs)

# Generate gs_per_ts for compatible_rstrs
gs <- paste0("g", 1:sum(MOIs)) # genotype names (graph vertices)
ts <- 1:length(MOIs) # episode indices
ts_per_gs <- rep(ts, MOIs) # episode index of each genotype
gs_per_ts <- split(gs, ts_per_gs) # genotypes grouped by episode

result <- compute_posterior(y_hom, fs, MOIs = MOIs, return.RG = TRUE, return.logp = TRUE)
CIL_gvn_RGs <- sapply(result$RGs, compatible_rstrs, gs_per_ts) # compatible states
CC_CI_CL_log <- t(sapply(CIL_gvn_RGs, function(x) c("CC", "CI", "CL") %in% x))
RGlogps <- sapply(result$RGs, function(x) x$logp)

c(SumRGs_CC = sum(exp(RGlogps[CC_CI_CL_log[,1]])),
  SumRGs_CI = sum(exp(RGlogps[CC_CI_CL_log[,2]])),
  SumRGs_CL = sum(exp(RGlogps[CC_CI_CL_log[,3]])))


sum(exp(RGlogps[CC_CI_CL_log[,1]]))
sum(exp(RGlogps[sapply(CIL_gvn_RGs, function(x) ("LC" %in% x) & (!"CC" %in% x))]))
sum(exp(RGlogps[sapply(CIL_gvn_RGs, function(x) "IC" %in% x)]))

sum(exp(RGlogps[CC_CI_CL_log[,2]]))
sum(exp(RGlogps[sapply(CIL_gvn_RGs, function(x) ("LI" %in% x) & (!"CI" %in% x))]))
sum(exp(RGlogps[sapply(CIL_gvn_RGs, function(x) "II" %in% x)]))

sum(exp(RGlogps[CC_CI_CL_log[,3]]))
sum(exp(RGlogps[sapply(CIL_gvn_RGs, function(x) ("LL" %in% x) & (!"CL" %in% x))]))
sum(exp(RGlogps[sapply(CIL_gvn_RGs, function(x) "IL" %in% x)]))

result$joint
result$marg
