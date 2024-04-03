# ==============================================================================
# aside: underflow / overflow problems with large marker counts
# aside: how to plot a graph for a real sample where gs and ts_per_gs are unknown? - record
# ==============================================================================
library(MCMCpack) # For rdirichlet
set.seed(1)

for(n_markers in c(3, 5, 10, 50)) { # Number of markers
  markers <- paste0("m", 1:n_markers) # Marker names
  alleles <- letters # Alleles
  n_alleles <- length(alleles) # Number of alleles per marker

  # Sample allele frequencies
  fs <- sapply(markers, function(m) {
    setNames(rdirichlet(1, alpha = rep(1, n_alleles)), alleles)
  }, USE.NAMES = TRUE, simplify = FALSE)

  # Sample parental genotypes
  parent1 <- sapply(markers, function(t) sample(alleles, 1), simplify = F)
  parent2 <- sapply(markers, function(t) sample(alleles, 1), simplify = F)
  parent3 <- sapply(markers, function(t) sample(alleles, 1), simplify = F)

  # Sample filial genotypes
  child12 <-  sapply(markers, function(t) sample(c(parent1[[t]], parent2[[t]]), 1), simplify = F)
  child23 <-  sapply(markers, function(t) sample(c(parent2[[t]], parent3[[t]]), 1), simplify = F)
  child13 <-  sapply(markers, function(t) sample(c(parent1[[t]], parent3[[t]]), 1), simplify = F)

  # Format genotypes for plot_data
  gs <- list(parent1 = list(parent1),
             parent2 = list(parent2),
             parent3 = list(parent3),
             child12 = list(child12),
             child13 = list(child13),
             child23 = list(child23))

  # Plot genotypes
  # plot_data(gs, fs = fs, marker_alleles = lapply(fs, names), marker_annotate = F)

  # Make parasite infections and then list of per-marker list of observed alleles
  initial <- rbind(unlist(child12), unlist(child13))
  relapse <- rbind(unlist(child23))
  y1 <- list(initial = apply(initial, 2, unique, simplify = F),
             relapse = apply(relapse, 2, unique, simplify = F))

  initial <- rbind(unlist(child12), unlist(child13))
  relapse <- rbind(unlist(child13), unlist(child23))
  y2 <- list(initial = apply(initial, 2, unique, simplify = F),
             relapse = apply(relapse, 2, unique, simplify = F))

  # Plot data
  ys <- list(y1 = y1, y2 = y2)
  plot_data(ys, fs = fs, marker_alleles = lapply(fs, names), marker_annotate = F)

  # Compute posteriors
  post1 <- compute_posterior(y1, fs, return.RG = TRUE)
  post2 <- compute_posterior(y2, fs, return.RG = TRUE)

  # Marginal posterior probabilities
  print(post1$marg[,c("L", "I")])
  print(post2$marg[,c("L", "I")])

  # Plot most likely graph
  par(mar = c(0.5, 0.5, 0.5, 0.5), mfrow = c(2,1))

  RGlogp <- sapply(post2$RGs, function(RG) RG$logp) # log likelihoods of graphs
  RG <- post2$RGs[[which.max(RGlogp)]]
  gs <- paste0("g", 1:4)
  ts_per_gs <- c(1,1,2,2)
  plot_RG(RG_to_igraph(RG, gs, ts_per_gs), edge.curved=0.25, vertex.size=20)
  seqs_comp_MLE_RG <- compatible_rstrs(RG, split(gs, ts_per_gs))

  RGlogp <- sapply(post1$RGs, function(RG) RG$logp) # log likelihoods of graphs
  RG <- post1$RGs[[which.max(RGlogp)]]
  gs <- paste0("g", 1:3)
  ts_per_gs <- c(1,1,2)
  plot_RG(RG_to_igraph(RG, gs, ts_per_gs), edge.curved=0.25, vertex.size=20)
  seqs_comp_MLE_RG <- compatible_rstrs(RG, split(gs, ts_per_gs))
}
