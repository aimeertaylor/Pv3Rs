# Run `Sample_example_ICL_3212.R` first to get data

library(Pv3Rs)
library(partitions)
library(tictoc)

y <- readRDS("data/y_ICL_3212.rds") # observed data
fs <- readRDS("data/fs_ICL_3212.rds") # allele frequencies
tic()
post <- compute_posterior(y, fs, return.RG=T) # defaults to uniform prior
toc()

writeLines("\nPosterior probability of recurrent states")
print(post$joint)

causes <- c("C","L","I")
for(i in 1:(length(y)-1)) {
  writeLines(paste("Recurrence", i))
  for(j in 1:3) {
    writeLines(paste0("Pr(R", i, "=", causes[j], "|Y)=",post$marg[i,j]))
  }
}

infection_count <- length(y)
MOIs <- determine_MOIs(y)
gs_count <- sum(MOIs)
gs <- paste0("g", 1:gs_count)
ts_per_gs <- rep(1:infection_count, MOIs)

# plot the RG that generated the observed data
par(mfrow=c(1,1))
RG0 <- readRDS("data/RG_ICL_3212.rds")
plot_RG(RG_to_igraph(RG0, gs, ts_per_gs), edge.curved=0.2)

# plot the top 9 RGs in terms of likelihood
par(mfrow=c(3,3), mar=c(1,1,1,1))
RG_order <- order(sapply(post$RGs, function(RG) RG$logp), decreasing=T)
for(RG_i in RG_order[1:9]) {
  plot_RG(RG_to_igraph(post$RGs[[RG_i]], gs, ts_per_gs), edge.curved=0.2)
}
