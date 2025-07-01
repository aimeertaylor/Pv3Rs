# ==============================================================================
# Made up example, based on modification of Figure 1 of the 2022 medRxiv preprint
# ==============================================================================
rm(list = ls())
library(Pv3Rs)
library(MCMCpack) # For rdirichlet
library(tictoc) # For timing
recombine_parent_ids <- utils::getFromNamespace("recombine_parent_ids", "Pv3Rs")

# Magic numbers / quantities
set.seed(5) # For reproducibility
n_alleles <- 10 # Number of alleles per marker (marker cardinality)
n_markers <- 8 # Number of markers
n_strangers <- 3 # Number of stranger parasites
c_param <- 1 # Dirichlet concentration parameter

# Derived quantities
alleles <- letters[1:n_alleles]
markers <- paste0("m", 1:n_markers) # Marker names

# Sample allele frequencies
fs <- sapply(markers, function(m) {
  sort(setNames(MCMCpack::rdirichlet(1, rep(c_param, n_alleles)), alleles), decreasing = T)
  }, USE.NAMES = T, simplify = F)
rm(alleles)

# Sample strangers
strangers <- sapply(1:n_strangers, function(i) {
  sapply(markers, function(t) sample(names(fs[[t]]), size = 1, prob = fs[[t]]))
})

# Designate strangers
parents <- strangers[, 1:2]

# Map the markers to chromosomes. Assume equally sized chromosomes — reasonable
# providing we later assume an equal number of crossovers per chromosome
chrs_per_marker <- round(seq(0.51, 14.5, length.out = n_markers))

# Sample parental allocations dependently (generates meiotic siblings)
cs <- recombine_parent_ids(chrs_per_marker)

# Construct children from parental allocations
children <- sapply(1:n_markers, function(i) {
  sapply(1:ncol(cs), function(j) parents[i,cs[i,j]])
})
colnames(children) <- markers

# Make parasite infections
initial <- apply(children[3:4,], 2, unique, simplify = F)
relapse <- apply(rbind(children[1:2,], strangers[,3]), 2, unique, simplify = F)

# Make data
y <- list(initial, relapse)

# Plot the data
plot_data(ys = list(example = y), fs = fs, marker.annotate = F)

# Reduce the number of markers evaluated
for(n_markers_eval in 1:n_markers){
  y_eval <- sapply(y, function(x) x[1:n_markers_eval], USE.NAMES = T, simplify = F)

  # Compute results
  post <- compute_posterior(y_eval, fs, return.RG=TRUE, return.logp=TRUE)
  colnames(post$marg) <- c("Recrudesence", "Relapse", "Reinfection")

  # Extract results
  round(post$marg,4) # Marginal probabilities
  which.max(post$joint) # Most likely sequence

  # Plot result on the simplex:
  xy <- apply(post$marg, 1, project2D)
  plot_simplex(v.labels = c("", "", ""), 0.75, c("blue", "purple", "red")) # make it so can pass more information
  points(x = xy["x",], xy["y",], pch = 21, cex = 4, bg = c("purple", "blue", "red"))
  text(x = xy["x",], xy["y",], labels = rep(n_markers_eval, 3), col = "white")
}

# Most likely graph
RGlogp <- sapply(post$RGs, function(RG) RG$logp)
RG <- post$RGs[[which.max(RGlogp)]]
MOIs <- determine_MOIs(y)
gs <- paste0("g", 1:sum(MOIs))
ts_per_gs <- rep(1:length(MOIs), MOIs)
par(mar = c(0.5, 0.5, 0.5, 0.5))
plot_RG(RG_to_igraph(RG, gs, ts_per_gs),
        vertex_palette = "Greys",
        edge.curved = 0.08,
        vertex.size = 40,
        edge.width = 3,
        vertex.frame.width = 2,
        vertex.label.cex = 4,
        vertex.label.color = "black")
seqs_comp_MLE_RG <- compatible_rstrs(RG, split(gs, ts_per_gs))
