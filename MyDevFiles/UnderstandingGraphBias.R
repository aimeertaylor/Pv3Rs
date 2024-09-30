################################################################################
# The output of compute posterior converges approximately to the prior-weighted
# relative fraction of relationship graphs compatible with the observed data.
#
# To do: add last two examples and understand why sometimes good, other times not
################################################################################
rm(list = ls())

# Simplex plot parameters
pardefault <- par()
vertex_names <- c(C = "Recrudescence", L = "Relapse", I = "Reinfection")
par(mar = c(0,0,0,0))

# Define function to get gs, RGs, states given RGs, and thus avoid repetition
get_graph_stuff <- function(x) {
  gs <- paste0("g", 1:sum(x)) # name genotypes
  ts <- 1:length(x) # episode indices
  ts_per_gs <- rep(ts, x) # episode index of each genotype
  gs_per_ts <- split(gs, ts_per_gs) # genotypes grouped by episode
  RGs <- enumerate_RGs(x) # generate all relationship graphs
  CIL_gvn_RGs <- sapply(RGs, compatible_rstrs, gs_per_ts) # list compatible states
  return(list(gs = gs, RGs = RGs, CIL_gvn_RGs = CIL_gvn_RGs))
}



#===============================================================================
# Homoallelic example without recurrent data using MOIs to change graphs.
#
# A rare homoallelic call has large effect (it limits the summation over
# relationship graphs to those with intra-episode sibling edges).
#
# Approx. based on graphs only is good because the likelihood surface is flat
# ===============================================================================
fs = list(m1 = c("A" = 0.001, "B" = 0.999)) # Allele frequencies
y <- list(enroll = list(m1 = c('A')), recur = list(m1 = NA)) # Data
MOIs <- list(c(2,1), c(3,1), c(2,2), c(3,2), c(3,3)) # MOIs

# ------------------------------------------------------------------------------
# Posterior based on the model
# ------------------------------------------------------------------------------
results <- sapply(MOIs, function(x) compute_posterior(y, fs, MOIs = x)$marg)
plot_simplex(c("Recrudescence", "Relapse", "Reinfection"), 0.5) # Plot simplex
xy <- apply(results, 1, project2D) # Project probabilities
points(x = xy["x", ], y = xy["y", ], pch = 20, cex = 1, col = 1:length(MOIs))
legend("left", col = 1:length(MOIs), pch = 20, title = "MOIs", inset = 0,
       legend = sapply(MOIs, paste, collapse = " & "), bty = "n")

# ------------------------------------------------------------------------------
# Posterior approximation based on graphs only
# ------------------------------------------------------------------------------
graph_frac_hom <- sapply(1:length(MOIs), FUN = function(i) {

  stuff <- get_graph_stuff(MOIs[[i]])
  RGs_C <- sapply(stuff$CIL_gvn_RGs, function(RG) "C" %in% RG)
  RGs_I <- sapply(stuff$CIL_gvn_RGs, function(RG) "I" %in% RG)

  # Identify RGs with intra-episode relatedness for the initial episode
  if(identical(MOIs[[i]], c(1,1))) {
    RGs_edges_first <- rep(TRUE, length(stuff$RGs)) # monoclonal case is trivial
  } else { # Make adjacency matrix for the first episode
    MOI1 <- MOIs[[i]][1]
    gs1 <- stuff$gs[1:MOI1]
    mat_first <- array(1, dim = c(MOI1, MOI1), dimnames = list(gs1, gs1))
    edges_first <- igraph::as_ids(igraph::E(igraph::graph_from_adjacency_matrix(mat_first, mode = "undirected", diag = F)))
    RGs_edges <- lapply(stuff$RGs, function(RG) igraph::as_ids(igraph::E(RG))) # Get edges for every RG
    RGs_edges_first <- sapply(RGs_edges, function(RG_E) all(edges_first %in% RG_E))
  }

  graph_frac_unnormalised <- c(frac_C = sum(RGs_edges_first*RGs_C)/sum(RGs_C),
                               frac_L = sum(RGs_edges_first)/length(stuff$RGs),
                               frac_I = sum(RGs_edges_first*RGs_I)/sum(RGs_I))
  graph_frac <- graph_frac_unnormalised/sum(graph_frac_unnormalised)
  return(graph_frac)
})
xy_graph_frac_hom <- apply(graph_frac_hom, 2, project2D)
points(x = xy_graph_frac_hom["x", ], y = xy_graph_frac_hom["y", ],
       col = 1:length(MOIs), cex = 2)
legend("bottom", pch = c(20, 1), pt.cex = c(1, 2), # Add legends
       legend = c("Posterior probability", "Relative weighted graphs"), bty = "n")

# ------------------------------------------------------------------------------
# Understanding approximation based on graphs only
# ------------------------------------------------------------------------------
MOIs <- c(3,3) # Illustrative example
gs <- paste0("g", 1:sum(MOIs))
ts <- 1:length(MOIs)
ts_per_gs <- rep(ts, MOIs)
gs_per_ts <- split(gs, ts_per_gs)
result33  <- compute_posterior(y, fs, MOIs = MOIs, return.logp = TRUE)
RGs <- sapply(result33$RGs, function(RG) RG_to_igraph(RG, gs, ts_per_gs))
logps <- sapply(result33$RGs, function(RG) RG$logp) # extract logp

# Make a vector of intra-infection edges by first creating a block diag. matrix
mat_first <- array(1, dim = c(3, 3), dimnames = list(paste0("g", 1:3), paste0("g", 1:3)))
edges_first <- igraph::as_ids(igraph::E(igraph::graph_from_adjacency_matrix(mat_first, mode = "undirected", diag = F)))
RGs_edges <- lapply(stuff$RGs, function(RG) igraph::as_ids(igraph::E(RG))) # Get edges for every RG
RGs_edges_first <- sapply(RGs_edges, function(RG_E) all(edges_first %in% RG_E))

# Inspect likelihoods:
if (max(logps[!(RGs_edges_first)]) < min(logps[RGs_edges_first])) {
  writeLines("Expected-most-likely graphs are most likely ")
  for(i in which(RGs_edges_first)) plot_RG(RGs[[i]])
  high_logps <- table(logps[RGs_edges_first])
  writeLines(sprintf("%s likelihood values among %s graphs", length(high_logps), sum(high_logps)))
  print(high_logps)
}











#===============================================================================
# Heteroallelic example without recurrent data using MOIs to change graphs.
#
# A Heteroallelic has a minor effect that doesn't change with allele frequency
# (it limits the summation over IBD partitions, but not the summation over
# relationship graphs)
#
# Graph expectation poor
#===============================================================================
fs = list(m1 = c("A" = 0.001, "B" = 0.999)) # Allele frequencies
y <- list(enroll = list(m1 = c('A','B')), recur = list(m1 = NA)) # Data
MOIs <- list(c(3,1), c(2,2), c(3,2)) # MOIs

# ------------------------------------------------------------------------------
# Posterior based on model
# ------------------------------------------------------------------------------
results <- sapply(MOIs, function(x) compute_posterior(y, fs, MOIs = x)$marg)
plot_simplex(c("Recrudescence", "Relapse", "Reinfection"), 0.5) # Plot simplex
xy <- apply(results, 2, project2D) # Project probabilities
points(x = xy["x", ], y = xy["y", ], pch = 20, cex = 1, col = 1:length(MOIs))
legend("left", col = 1:length(MOIs), pch = 20, title = "MOIs", inset = 0,
       legend = sapply(MOIs, paste, collapse = " & "), bty = "n")

# ------------------------------------------------------------------------------
# Posterior approximation based on graphs only
# ------------------------------------------------------------------------------
graph_frac_het <- sapply(1:length(MOIs), FUN = function(i) {
  stuff <- get_graph_stuff(MOIs[[i]])
  RGs_C <- sapply(stuff$CIL_gvn_RGs, function(RG) "C" %in% RG)
  RGs_I <- sapply(stuff$CIL_gvn_RGs, function(RG) "I" %in% RG)
  graph_frac_unnormalised <- c(frac_C = sum(RGs_C)/sum(RGs_C),
                               frac_L = 1,
                               frac_I = sum(RGs_I)/sum(RGs_I))
  graph_frac <- graph_frac_unnormalised/sum(graph_frac_unnormalised)
  return(graph_frac)
})
xy_graph_frac_het <- apply(graph_frac_het, 2, project2D)
points(x = xy_graph_frac_het["x", ], y = xy_graph_frac_het["y", ],
       col = (1:length(MOIs))[-1], cex = 2)
legend("bottom", pch = c(20, 1), pt.cex = c(1, 2), # Add legends
       legend = c("Posterior probability", "Relative weighted graphs"), bty = "n")





#===============================================================================
# Homoallelic example with recurrent data using recurrences to change graphs.
#
# The marginal probability that the first recurrence is a recrudescence
# increases at a decreasing rate as the graph grows without data. Decreasing the
# frequency of the observed allele decreases the probability of reinfection.
#
# Graph expectation poor
#===============================================================================
# Data for different recurrence counts where only the 1st recurrence has data
ys_match <- list(one_recurrence = list(enroll = list(m1 = "A"),
                                       recur_match = list(m1 = "A")),
                 two_recurrences = list(enroll = list(m1 = "A"),
                                        recur_match = list(m1 = "A"),
                                        recur2 = list(m1 = NA)),
                 three_recurrences = list(enroll = list(m1 = "A"),
                                          recur_match = list(m1 = "A"),
                                          recur2 = list(m1 = NA),
                                          recur3 = list(m1 = NA)),
                 four_recurrences = list(enroll = list(m1 = "A"),
                                         recur_match = list(m1 = "A"),
                                         recur2 = list(m1 = NA),
                                         recur3 = list(m1 = NA),
                                         recur4 = list(m1 = NA)))
fs <- list(m1 = c("A" = 0.001, "B" = 0.999)) # Allele frequencies
results <- lapply(ys_match, function(y) suppressMessages(compute_posterior(y, fs)$marg))
results_first <- sapply(results, function(result) result[1,]) # Extract 1st recurrence results
plot_simplex(vertex_names, 0.5) # Plot 2D simplex
xy <- apply(results_first, 2, project2D)
points(x = xy["x", ], y = xy["y", ], pch = "-", col = 1:length(ys_match), cex = 1)

# ------------------------------------------------------------------------------
# Add graph expectation for a rare repeat allele
# ------------------------------------------------------------------------------
MOIs <- lapply(ys_match, determine_MOIs) # First, get MOIs
graph_frac_first <- sapply(1:length(MOIs), FUN = function(i) {
  stuff <- get_graph_stuff(MOIs[[i]])
  # RGs compatible with recrudescence at the first recurrence
  RGs_C_first <- sapply(stuff$CIL_gvn_RGs, function(x) "C" %in% do.call(rbind, strsplit(x, split = ""))[,1])
  graph_frac_unnormalised <- c(frac_C = 1/sum(RGs_C_first),
                               frac_L = 1/length(stuff$RGs),
                               frac_I = 0)
  graph_frac <- graph_frac_unnormalised/sum(graph_frac_unnormalised)
  return(graph_frac)
})
xy <- apply(graph_frac_first, 2, project2D) # project onto the simplex
points(x = xy["x", ], y = xy["y", ], col = 1:length(ys_match), cex = 1) # plot
legend("right", col = 1:length(ys_match), pch = "-", pt.cex = 2, bty = "n",
       legend = 1:length(ys_match), title = "Recurrence \n count")












#===============================================================================
# Reinfection versus relapse using MOIs to change graphs.
#
# When data are informative, model returns the prior on graphs re-weighted
# ===============================================================================
set.seed(1)
marker_count <- 100; allele_count <- 2 # Data informativeness
ms <- paste0("m", 1:marker_count)
alleles <- LETTERS[1:allele_count]
fs <- sapply(ms, function(m) {
  setNames(MCMCpack::rdirichlet(1, rep(1, allele_count)), alleles)
}, USE.NAMES = TRUE, simplify = FALSE)
stranger1 <- as.list(sapply(ms, function(t) sample(alleles, 1, prob = fs[[t]])))
stranger2 <- as.list(sapply(ms, function(t) sample(alleles, 1, prob = fs[[t]])))
y <- list(enrol = stranger1, recur = stranger1)
MOIs <- list(c(1,1), c(2,1), c(3,1), c(2,2), c(3,2), c(3,3), c(4,2))

# ------------------------------------------------------------------------------
# Posterior based on full model
# ------------------------------------------------------------------------------
results <- sapply(MOIs, function(x) compute_posterior(y, fs, MOIs = x)$marg)
plot_simplex(vertex_names, 0.5); xy <- apply(results, 2, project2D)
points(x = xy["x", ], y = xy["y", ], col = 1:length(MOIs), pch = 20)
legend("left", bty = "n", col = 1:length(MOIs), pch = 20, title = "MOIs",
       legend = sapply(MOIs, paste, collapse = " & "))

# ------------------------------------------------------------------------------
# Posterior approximation based on graphs only: prior on graphs re-weighted
# ------------------------------------------------------------------------------
graph_frac_IL <- sapply(1:length(MOIs), FUN = function(i) {
  stuff <- get_graph_stuff(MOIs[[i]])
  RGs_I <- sapply(stuff$CIL_gvn_RGs, function(RG) "I" %in% RG)
  x <- c(frac_C = 0, frac_L = 1/length(stuff$RGs), frac_I = 1/sum(RGs_I))
  z <- x/sum(x)
})
xy_graph_frac <- apply(graph_frac_IL, 2, project2D)
points(x = xy_graph_frac["x", ], y = xy_graph_frac["y", ],
       col = 1:length(MOIs), pch = 1, cex = 2)
legend("bottom", bty = "n", pch = c(20, 1), pt.cex = c(1,2),
       legend = c("Posterior probability", "Relative graph proportion"), )



#===============================================================================
# Recrudescence vs relapse using MOIs to change graphs
#
# When data are informative, model returns the prior on graphs re-weighted
# ===============================================================================
set.seed(1)
marker_count <- 100; allele_count <- 10 # Data informativeness
ms <- paste0("m", 1:marker_count)
alleles <- LETTERS[1:allele_count]
fs <- sapply(ms, function(m) {
  setNames(MCMCpack::rdirichlet(1, rep(1, allele_count)), alleles)
}, USE.NAMES = TRUE, simplify = FALSE)
clone <- as.list(sapply(ms, function(t) sample(alleles, 1, prob = fs[[t]])))
y <- list(enrol = clone, recur = clone)
MOIs <- list(c(1,1), c(2,1), c(3,1), c(2,2), c(3,2), c(3,3))

# ------------------------------------------------------------------------------
# Posterior based on full model
# ------------------------------------------------------------------------------
results <- sapply(MOIs, function(x) compute_posterior(y, fs, MOIs = x)$marg)
xy <- apply(results, 2, project2D)
plot_simplex(vertex_names, 0.5)
points(x = xy["x", ], y = xy["y", ], pch = 20, col = 1:length(MOIs))
legend("left", bty = "n", col = 1:length(MOIs), pch = 20, title = "MOIs",
       legend = sapply(MOIs, paste, collapse = " & "))

# -------------------------------------------------------------------------------
# Posterior approximation based on graphs only: prior on graphs re-weighted
# -------------------------------------------------------------------------------
graph_frac_CL <- sapply(1:length(MOIs), FUN = function(i) {
  stuff <- get_graph_stuff(MOIs[[i]])
  RGs_C <- sapply(stuff$CIL_gvn_RGs, function(RG) "C" %in% RG)
  x <- c(frac_C = 1/sum(RGs_C), frac_L = 1/length(stuff$RGs), frac_I = 0)
  z <- x/sum(x)
})
xy_graph_frac <- apply(graph_frac_CL, 2, project2D)
points(x = xy_graph_frac["x", ], y = xy_graph_frac["y", ],
       col = 1:length(MOIs), pch = 1, cex = 2)
legend("bottom", bty = "n", pch = c(20, 1), pt.cex = c(1,2),
       legend = c("Posterior probability", "Relative graph proportion"), )


# ------------------------------------------------------------------------------
# Graphs with notable likelihood in the MOI c(3,3) case are limited to those
# with intra-infection relatedness compatible with recrudescence or relapse but
# not reinfection; there are multiple likelihood values among them.
# ------------------------------------------------------------------------------
MOIs <- c(3,3) # Illustrative example
gs <- paste0("g", 1:sum(MOIs))
ts <- 1:length(MOIs)
ts_per_gs <- rep(ts, MOIs)
gs_per_ts <- split(gs, ts_per_gs)
result33  <- compute_posterior(y, fs, MOIs = MOIs, return.logp = TRUE)
RGs <- sapply(result33$RGs, function(RG) RG_to_igraph(RG, gs, ts_per_gs))
CIL_gvn_RGs <- sapply(RGs, compatible_rstrs, gs_per_ts = gs_per_ts) # states
RGs_C_or_L <- sapply(CIL_gvn_RGs, function(CIL_gvn_RG) !("I" %in% CIL_gvn_RG))

# Make a vector of intra-infection edges by first creating a block diag. matrix
mat_within <- Matrix::bdiag(lapply(MOIs,function(x) matrix(1, ncol=x, nrow=x)))
colnames(mat_within) <- gs; rownames(mat_within) <- gs
edges_within <- igraph::as_ids(igraph::E(igraph::graph_from_adjacency_matrix(mat_within, mode = "undirected", diag = F)))
RGs_edges <- lapply(RGs, function(RG) igraph::as_ids(igraph::E(RG)))
RGs_edges_within <- sapply(RGs_edges, function(RG_E) all(edges_within %in% RG_E))
logps <- sapply(result33$RGs, function(RG) RG$logp) # extract logp

# Inspect likelihoods:
if (max(logps[!(RGs_C_or_L & RGs_edges_within)]) < min(logps[RGs_C_or_L & RGs_edges_within])) {
  writeLines("Expected-most-likely graphs are most likely ")
  for(i in which(RGs_C_or_L & RGs_edges_within)) plot_RG(RGs[[i]])
  high_logps <- table(logps[RGs_C_or_L & RGs_edges_within])
  writeLines(sprintf("%s likelihood values among %s graphs", length(high_logps), sum(high_logps)))
  print(high_logps)
}











