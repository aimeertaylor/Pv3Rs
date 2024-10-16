# ==============================================================================
# Made up example, based on modification of Figure 1 of the 2022 medRxiv preprint
# ==============================================================================
rm(list = ls())
# library(Pv3Rs)
# library(MCMCpack) # For rdirichlet
# library(tictoc) # For timing

# Magic numbers / quantities
set.seed(5) # For reproducibility
n_alleles <- 10 # Number of alleles per marker (marker cardinality)
n_markers <- 10 # Number of markers
n_strangers <- 4 # Number of stranger parasites
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
relapse <- apply(rbind(children[1,], strangers[,3]), 2, unique, simplify = F)
recrude <- as.list(strangers[,3])
reinfec <- as.list(strangers[,4])

# Make data
y <- list("Day 0 episode" = initial,
          "Day 40 recurrence" = relapse,
          "Day 70 recurrence" = recrude,
          "Day 200 recurrence" = reinfec)

# Make phased data
y_phased <- list(as.list(children[3,]),
                 as.list(children[4,]),
                 as.list(children[1,]),
                 as.list(strangers[,3]))

# Plot the data
plot_data(ys = list(example = y_phased), fs = fs, marker_annotate = F)
plot_data(ys = list(example = y), fs = fs, marker_annotate = F)

# Reduce the number of markers evaluated
for(n_markers_eval in 1:10){
  y_eval <- sapply(y, function(x) x[1:n_markers_eval], USE.NAMES = T, simplify = F)

  # Compute results
  post <- compute_posterior(y_eval, fs, return.RG=TRUE, return.logp=TRUE)
  colnames(post$marg) <- c("Recrudesence", "Relapse", "Reinfection")

  # Extract results
  round(post$marg,4) # Marginal probabilities
  which.max(post$joint) # Most likely sequence

  # Plot result on the simplex:
  xy <- apply(post$marg, 1, project2D)
  plot_simplex(v_labels = c("", "", ""), 0.75, c("blue", "purple", "red")) # make it so can pass more information
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
plot_RG(RG_to_igraph(RG, gs, ts_per_gs), edge.curved=0, vertex.size=20, vertex_palette = "Blues")
seqs_comp_MLE_RG <- compatible_rstrs(RG, split(gs, ts_per_gs))

y <- list(list(m1 = c("A", "B")), list(m1 = "A"))
fs <- list(m1 = c("A" = 0.9, "B" = 0.1))
ps <- compute_posterior(y, fs, return.logp = T)
ps$marg
plot_RG(RG_to_igraph(ps$RGs[[which.max(sapply(ps$RGs, function(x) x$logp))]],
             gs = paste0("g", 1:3), ts_per_gs = c(1,1,2)))
