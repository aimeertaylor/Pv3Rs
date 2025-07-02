################################################################################
# Investigate posterior probability bounds for two recurrences
################################################################################

library(matrixStats) # for logSumExp
library(Pv3Rs)

marker_count <- 100 # Number of markers
ms <- paste0("m", 1:marker_count) # Marker names
all_As <- sapply(ms, function(t) "A", simplify = F) # As for all markers
all_AB <- sapply(ms, function(t) c("A","B"), simplify = F) # A+B for all markers
all_Bs <- sapply(ms, function(t) "B", simplify = F) # Bs for all markers
no_data <- sapply(ms, function(t) NA, simplify = F) # NAs for all markers
fs <- sapply(ms, function(m) c("A" = 0.01, "B" = 0.99), simplify = FALSE)
# NB: enrolment now has MOI=2, but still there i strong evidence of recrudescence
y <- list("enrol" = all_AB,
          "recur1" = all_As,
          "recur2" = no_data)

post <- compute_posterior(y, fs, return.RG = TRUE, return.logp = TRUE)
lliks <- sapply(post$RGs, function(RG) RG$logp)
gs <- paste0("g", 1:4)
MOIs <- determine_MOIs(y)
ts_per_gs <- rep(1:length(y), MOIs)
gs_per_ts <- split(paste0("g", 1:sum(MOIs)), rep(1:length(MOIs), MOIs))

rstrs_list <- lapply(post$RGs, Pv3Rs:::compatible_rstrs, gs_per_ts)
n_RGs <- length(post$RGs)


## Checking validity of (4)

rstrs <- names(post$joint)
idxs_list <- lapply(setNames(rstrs, rstrs) ,
                    function(rstr) which(sapply(1:n_RGs, function(i) rstr%in%rstrs_list[[i]])))
idxs_list[["CC"]] # indices of graphs compatible with CC

# (4) assumes sum of p(y|g) are not equal for g \in G_CC, g \in G_CL, g \in G_CI
# Sum of p(y|g) are not equal for g \in G_CC, g \in G_CL, g \in G_CI, there's about a 5-fold difference
logSumExp(sapply(idxs_list[["CC"]], function(i) post$RGs[[i]]$logp))
logSumExp(sapply(idxs_list[["CL"]], function(i) post$RGs[[i]]$logp))
logSumExp(sapply(idxs_list[["CI"]], function(i) post$RGs[[i]]$logp))

# Averages are closer, but still not exactly equal
logAvgExp <- function(x) logSumExp(x) - log(length(x))
logAvgExp(sapply(idxs_list[["CC"]], function(i) post$RGs[[i]]$logp))
logAvgExp(sapply(idxs_list[["CL"]], function(i) post$RGs[[i]]$logp))
logAvgExp(sapply(idxs_list[["CI"]], function(i) post$RGs[[i]]$logp))
# It's unclear what conditions would give close averages

# Check (4) itself
3/(
  3+
  length(idxs_list[["CC"]])/length(idxs_list[["LC"]])+
  length(idxs_list[["CL"]])/length(idxs_list[["LL"]])+
  length(idxs_list[["CI"]])/length(idxs_list[["LI"]])
)
post$marg["recur1","C"] # the result is close, but this is larger than (4) (no longer an upper bound)


## Check the 'converges to' result in Sec 1.1.2

(CC_bound <- 1/length(idxs_list[["CC"]]) / sum(
  1/length(idxs_list[["CC"]]),
  1/length(idxs_list[["LC"]]),
  1/length(idxs_list[["CL"]]),
  1/length(idxs_list[["LL"]])
))
(CL_bound <- 1/length(idxs_list[["CL"]]) / sum(
  1/length(idxs_list[["CL"]]),
  1/length(idxs_list[["LL"]])
))
(CI_bound <- 1/length(idxs_list[["CI"]]) / sum(
  1/length(idxs_list[["CI"]]),
  1/length(idxs_list[["LI"]]),
  1/length(idxs_list[["CL"]]),
  1/length(idxs_list[["LL"]])
))

# Upper bounds are satisfied, but the bounds are far from equality
post$joint[c("CC","CL","CI")]
