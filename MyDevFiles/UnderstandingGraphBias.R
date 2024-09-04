################################################################################
# General to-dos:
# Does it make sense to have both return.RG and return.logp?
# Add a warning for unpaired data
#
# Examples of graph prior effects on posterior
# Make this into a vignette: need to compute RG after ruling out a state
#
# When genetic data can only rule out a state, the output of compute_posterior
# is approximately equal to the relative prior proportion of compatible graphs
# given MOIs.
#
# explain multiple recurrences with relative prior proportion
################################################################################
pardefault <- par()

# Simplex plot parameters
par(mfrow = c(1,2), mar = c(0,0,0,0))
vertex_names <- c(C = "Recrudescence", L = "Relapse", I = "Reinfection")

#============================================================================
# Example for data on a single episode using MOIs to change graph structure
#
# Pv3Rs is designed to analyse data on two or more episodes.
#
# The output of compute posterior is approximately equal to the relative prior
# proportion of relationship graphs compatible with the observed data. Is this
# ever not the case?
#
# A homoallelic call has large effect that is amplified when the allele of the
# homoallelic call is rare.
#
# A heteroallelic call has a minor effect that doesn't change with the allele
# frequency of the rare allele.
#
# A single heteroallelic observation changes the
# summation over IBD partitions but not the summation over relationship graphs
# and this has a smaller effect.
#============================================================================
f_rare <- 0.001
fs = list(m1 = c('1' = f_rare, '2' = 1-f_rare))
MOIs <- list(c(1,1), c(2,1), c(3,1), c(2,2), c(3,2), c(3,3))

# Data: a homoallelic or heteroallelic call for a single unpaired episode
y_hom <- list(enroll = list(m1 = c('1')), recur = list(m1 = NA))
y_het <- list(enroll = list(m1 = c('1','2')), recur = list(m1 = NA))

# Compute posterior probabilities and extract marginal probabilities:
results_hom <- do.call(rbind, lapply(MOIs, function(x) suppressMessages(compute_posterior(y_hom, fs, MOIs = x)$marg)))
results_het <- do.call(rbind, lapply(MOIs[-1], function(x) suppressMessages(compute_posterior(y_het, fs, MOIs = x)$marg)))

# Project probabilities onto 2D simplex coordinates
xy_hom <- apply(results_hom, 1, project2D)
xy_het <- apply(results_het, 1, project2D)

# Plot 2D simplex:het
plot_simplex(v_labels = vertex_names[colnames(results_het)], classifcation_threshold = 0.5)
points(x = xy_het["x", ], y = xy_het["y", ], pch = 20, cex = 0.25)
title(main = "Heteroallelic", line = -5)

# Plot 2D simplex:hom
plot_simplex(v_labels = vertex_names[colnames(results_hom)], classifcation_threshold = 0.5)
points(x = xy_hom["x", ], y = xy_hom["y", ], pch = 20, cex = 0.25)
title(main = "Homoallelic", line = -5)

# ------------------------------------------------------------------------------
# Understanding the notable impact of the homoallelic call
# ------------------------------------------------------------------------------
graph_frac_hom <- sapply(1:length(MOIs), FUN = function(i) {

  gs <- paste0("g", 1:sum(MOIs[[i]])) # name genotypes
  ts <- 1:length(MOIs[[i]]) # episode indices
  ts_per_gs <- rep(ts, MOIs[[i]]) # episode index of each genotype
  gs_per_ts <- split(gs, ts_per_gs) # genotypes grouped by episode
  RGs <- enumerate_RGs(MOIs[[i]]) # generate all relationship graphs
  CIL_gvn_RGs <- sapply(RGs, compatible_rstrs, gs_per_ts) # list compatible states
  RGs_C <- sapply(CIL_gvn_RGs, function(RG) "C" %in% RG) # RGs compatible with recrudescence
  RGs_I <- sapply(CIL_gvn_RGs, function(RG) "I" %in% RG) # RGs compatible with reinfection

  # Since a rare homoallelic call is best explained by intra-episode relatedness
  # Find RGs with intra-episode relatedness
  if(identical(MOIs[[i]], c(1,1))) {
    RGs_edges_within <- rep(TRUE, length(RGs)) # monoclonal case is trivial
  } else {
    # Make adjecency matrix with intra-episode edges for the first episode
    mat_within <- matrix(0, ncol = sum(MOIs[[i]]), nrow = sum(MOIs[[i]]), dimnames = list(gs, gs)) # initiate matrix
    mat_within[1:MOIs[[i]][1], 1:MOIs[[i]][1]] <- 1 # populate matrix
    edges_within <- igraph::as_ids(igraph::E(igraph::graph_from_adjacency_matrix(mat_within, mode = "undirected", diag = F)))
    RGs_edges <- lapply(RGs, function(RG) igraph::as_ids(igraph::E(RG))) # Get edges for every RG
    RGs_edges_within <- sapply(RGs_edges, function(RG_E) all(edges_within %in% RG_E))
  }

  graph_frac_unnormalised <- c(frac_C = sum(RGs_edges_within*RGs_C)/sum(RGs_C),
                                  frac_L = sum(RGs_edges_within)/length(RGs),
                                  frac_I = sum(RGs_edges_within*RGs_I)/sum(RGs_I))

  graph_frac <- graph_frac_unnormalised/sum(graph_frac_unnormalised)

  return(graph_frac)
})

# Add to relative prior proportion to simplex
xy_graph_frac_hom <- apply(graph_frac_hom, 2, project2D)
points(x = xy_graph_frac_hom["x", ], y = xy_graph_frac_hom["y", ], pch = 1)

# ------------------------------------------------------------------------------
# Understanding the negligible impact of the heteroallelic call
# ------------------------------------------------------------------------------
graph_frac_het <- sapply(1:length(MOIs), FUN = function(i) {

  gs <- paste0("g", 1:sum(MOIs[[i]]))
  ts <- 1:length(MOIs[[i]])
  ts_per_gs <- rep(ts, MOIs[[i]])
  gs_per_ts <- split(gs, ts_per_gs)
  RGs <- enumerate_RGs(MOIs[[i]])
  CIL_gvn_RGs <- sapply(RGs, compatible_rstrs, gs_per_ts = gs_per_ts) # states
  RGs_C <- sapply(CIL_gvn_RGs, function(RG) "C" %in% RG)
  RGs_I <- sapply(CIL_gvn_RGs, function(RG) "I" %in% RG)

  graph_frac_unnormalised <- c(frac_C = sum(RGs_C)/sum(RGs_C),
                               frac_L = 1,
                               frac_I = sum(RGs_I)/sum(RGs_I))

  graph_frac <- graph_frac_unnormalised/sum(graph_frac_unnormalised)
  return(graph_frac)
})

# Almost identical
results_het
t(graph_frac_het)

#===============================================================================
# Example for data on a single episode using recurrences to change graph structure
#
# Changing the graph structure by adding recurrences has very little impact when
# data are limited to a single episode.
#===============================================================================
# Allele frequencies:
f_rare <- 0.25
fs <- list(m1 = setNames(c(f_rare, 1-f_rare), c("A", "Other")))

# The number of recurrences increases but only the first recurrence has data
ys_homo <- list(list(list(m1 = "A"), list(m1 = NA)), # 1 recurrence
                list(list(m1 = "A"), list(m1 = NA), list(m1 = NA)), # 2 recurrences
                list(list(m1 = "A"), list(m1 = NA), list(m1 = NA), list(m1 = NA)), # etc.
                list(list(m1 = "A"), list(m1 = NA), list(m1 = NA), list(m1 = NA), list(m1 = NA)))

# The number of recurrences increases and all have the same data
ys_het <- list(list(list(m1 = c("A", "Other")), list(m1 = NA)),  # 1 recurrence
               list(list(m1 = c("A", "Other")), list(m1 = NA), list(m1 = NA)), # 2 recurrences
               list(list(m1 = c("A", "Other")), list(m1 = NA), list(m1 = NA), list(m1 = NA)), # etc.
               list(list(m1 = c("A", "Other")), list(m1 = NA), list(m1 = NA), list(m1 = NA), list(m1 = NA)))

# Compute posterior probabilities and extract marginal probabilities:
results_hom <- lapply(ys_homo, function(y) compute_posterior(y, fs)$marg)
results_het <- lapply(ys_het, function(y) compute_posterior(y, fs)$marg)

# Extract results for the first recurrence only:
first_recur_hom <- sapply(results_hom, function(result) result[1,])
first_recur_het <- sapply(results_het, function(result) result[1,])

# Project probabilities onto 2D simplex coordinates
xy_hom <- apply(first_recur_homo, 2, project2D)
xy_het <- apply(first_recur_het, 2, project2D)

# Plot 2D simplex:het
plot_simplex(v_labels = vertex_names[colnames(results_het[[1]])], classifcation_threshold = 0.5)
points(x = xy_hom["x", ], y = xy_hom["y", ], pch = 20, cex = 0.25, col = "red")
points(x = xy_het["x", ], y = xy_het["y", ], pch = 20, cex = 0.25, col = "blue")
title(main = "Heteroallelic", line = -5)
legend("topright", inset = 0.1, col = c("red", "blue"), pch = 20, bty = "n",
       legend = c("Homoallelic", "Heteroallelic"))

#===============================================================================
# Example with recurrent data using multiple recurrences to change graph structure
#
# The marginal probability that the first recurrence is a recrudescence
# increases as the graph grows without data (red arrow on simplex plot).
#
# The marginal probability that the first recurrence is a reinfection increases
# as the graph grows with repeat data (blue arrow on simplex plot).
#
# In both cases, the  marginal probability that the first recurrence is a
# relapse decreases (both arrows on the simplex plot).
#
# These general observations are not very sensitive to the observed allele
# frequency (try comparing Obs_allele equal to 0.05, 0.5 and 0.95).
#
# Also of note, the marginal probabilities of recurrences with no data diverge
# from the prior at a decreasing rate (compare rows 2 to 4 of results0[[4]]).
#
# Also of note, the marginal probabilities of different recurrences differ even
# when all recurrences have the same data (compare rows 2 to 4 of
# results1[[4]]).
#===============================================================================

# Allele frequencies:
f_rare <- 0.25
fs <- list(m1 = setNames(c(f_rare, 1-f_rare), c("A", "Other")))

# The number of recurrences increases but only the first recurrence has data
ys_one <- list(list(list(m1 = "A"), list(m1 = "A")), # 1 recurrence
            list(list(m1 = "A"), list(m1 = "A"), list(m1 = NA)), # 2 recurrences
            list(list(m1 = "A"), list(m1 = "A"), list(m1 = NA), list(m1 = NA)), # etc.
            list(list(m1 = "A"), list(m1 = "A"), list(m1 = NA), list(m1 = NA), list(m1 = NA)))

# The number of recurrences increases and all have the same data
ys_all <- list(list(list(m1 = "A"), list(m1 = "A")),  # 1 recurrence
            list(list(m1 = "A"), list(m1 = "A"), list(m1 = "A")), # 2 recurrences
            list(list(m1 = "A"), list(m1 = "A"), list(m1 = "A"), list(m1 = "A")), # etc.
            list(list(m1 = "A"), list(m1 = "A"), list(m1 = "A"), list(m1 = "A"), list(m1 = "A")))

# Compute posterior probabilities and extract marginal probabilities:
results_one <- lapply(ys_one, function(y) compute_posterior(y, fs)$marg)
results_all <- lapply(ys_all, function(y) compute_posterior(y, fs)$marg)

# Extract results for the first recurrence only:
first_recur_one <- sapply(results_one, function(result) result[1,])
first_recur_all <- sapply(results_all, function(result) result[1,])

# ------------------------------------------------------------------------------
# Plot divergence on 2D simplex
# ------------------------------------------------------------------------------
n_recur <- ncol(first_recur0)
plot_simplex(v_labels = vertex_names, classifcation_threshold = 0.5)
legend("topright", inset = 0.1, col = c("red", "blue"), pch = 20, bty = "n",
       legend = c("Only the first recurrence has data",
                  "All recurrences have the same data"))

# Project and plot first_recur0
xy <- apply(first_recur_one, 2, project2D)
arrows(x0 = xy["x", 1], x1 = xy["x", n_recur],
       y0 = xy["y", 1], y1 = xy["y", n_recur],
       length = 0.05, col = "red")
points(x = xy["x", ], y = xy["y", ], pch = ".")

# Project and plot first_recur1
xy <- apply(first_recur_all, 2, project2D)
arrows(x0 = xy["x", 1], x1 = xy["x", n_recur],
       y0 = xy["y", 1], y1 = xy["y", n_recur],
       length = 0.05, col = "blue")
points(x = xy["x", ], y = xy["y", ], pch = ".")


# The marginal probabilities of recurrences with no data diverge from the prior
# at a decreasing rate, at least for relapse and reinfection:
results_one

# The marginal probabilities of different recurrences differ when they all
# recurrences have the same data:
results_all

# ------------------------------------------------------------------------------
# Understanding the notable impact of the homoallelic call
# ------------------------------------------------------------------------------
graph_frac_hom <- sapply(1:length(MOIs), FUN = function(i) {

  gs <- paste0("g", 1:sum(MOIs[[i]])) # name genotypes
  ts <- 1:length(MOIs[[i]]) # episode indices
  ts_per_gs <- rep(ts, MOIs[[i]]) # episode index of each genotype
  gs_per_ts <- split(gs, ts_per_gs) # genotypes grouped by episode
  RGs <- enumerate_RGs(MOIs[[i]]) # generate all relationship graphs
  CIL_gvn_RGs <- sapply(RGs, compatible_rstrs, gs_per_ts) # list compatible states
  RGs_C <- sapply(CIL_gvn_RGs, function(RG) "C" %in% RG) # RGs compatible with recrudescence
  RGs_I <- sapply(CIL_gvn_RGs, function(RG) "I" %in% RG) # RGs compatible with reinfection

  # Since a rare homoallelic call is best explained by intra-episode relatedness
  # Find RGs with intra-episode relatedness
  if(identical(MOIs[[i]], c(1,1))) {
    RGs_edges_within <- rep(TRUE, length(RGs)) # monoclonal case is trivial
  } else {
    # Make adjecency matrix with intra-episode edges for the first episode
    mat_within <- matrix(0, ncol = sum(MOIs[[i]]), nrow = sum(MOIs[[i]]), dimnames = list(gs, gs)) # initiate matrix
    mat_within[1:MOIs[[i]][1], 1:MOIs[[i]][1]] <- 1 # populate matrix
    edges_within <- igraph::as_ids(igraph::E(igraph::graph_from_adjacency_matrix(mat_within, mode = "undirected", diag = F)))
    RGs_edges <- lapply(RGs, function(RG) igraph::as_ids(igraph::E(RG))) # Get edges for every RG
    RGs_edges_within <- sapply(RGs_edges, function(RG_E) all(edges_within %in% RG_E))
  }

  graph_frac_unnormalised <- c(frac_C = sum(RGs_edges_within*RGs_C)/sum(RGs_C),
                               frac_L = sum(RGs_edges_within)/length(RGs),
                               frac_I = sum(RGs_edges_within*RGs_I)/sum(RGs_I))

  graph_frac <- graph_frac_unnormalised/sum(graph_frac_unnormalised)

  return(graph_frac)
})

#============================================================================
# Example with recurrent data using MOIs to change graph structure
#
# Having ruled out reinfection, posterior odds of relapse versus recrudescence
# depends on graphs. Having ruled out recrudescence, posterior odds of relapse
# verus reinfection depends on graphs. This is easiest to see for relapse versus
# recrudescence using one very rare allele (as informative as many alleles and
# quicker); see example below.
#============================================================================
# Frequencies, data and MOIs
f_rare <- 0.001 # Make reinfection unlikely, recrudescence likely
fs = list(m1 = c('1' = f_rare, '2' = 1-f_rare))
y <- list(enrol = list(m1 = "1"), recur = list(m1 = c("1")))
MOIs <- list(c(1,1), c(2,1), c(3,1), c(2,2), c(3,2), c(3,3), c(4,2))

# Compute posterior probabilities, extract marginal probabilities and project
prior <- matrix(rep(1/3,3), nrow = 1, dimnames = list(NULL,c("C", "L", "I")))
results <- do.call(rbind, lapply(MOIs, function(x) suppressMessages(compute_posterior(y, fs, MOIs = x, prior = prior)$marg)))
xy <- apply(results, 1, project2D) # Project probabilities onto 2D simplex coordinates

# Plot 2D simplex
pardefault <- par()
par(mar = c(0,0,0,0))
vertex_names <- c(C = "Recrudescence", L = "Relapse", I = "Reinfection")
plot_simplex(v_labels = vertex_names[colnames(results)], classifcation_threshold = 0.5)
points(x = xy["x", ], y = xy["y", ], pch = as.character(1:length(MOIs)), cex = 0.25)
par(mar = pardefault$mar) # Restore plotting margins

# ------------------------------------------------------------------------------
# Understanding the relapse vs recrudescence results...
#
# The relative proportion of intra-infection related graphs compatible with
# recrudescence versus relapse but not reinfection is a function of the MOI
# ------------------------------------------------------------------------------
graph_frac_CL <- sapply(1:length(MOIs), FUN = function(i) {

  gs <- paste0("g", 1:sum(MOIs[[i]]))
  ts <- 1:length(MOIs[[i]])
  ts_per_gs <- rep(ts, MOIs[[i]])
  gs_per_ts <- split(gs, ts_per_gs)
  RGs <- enumerate_RGs(MOIs[[i]])
  CIL_gvn_RGs <- sapply(RGs, compatible_rstrs, gs_per_ts = gs_per_ts) # states
  RGs_C_or_L <- sapply(CIL_gvn_RGs, function(RG) !("I" %in% RG))
  RGs_C <- sapply(CIL_gvn_RGs, function(RG) "C" %in% RG)

  if(identical(MOIs[[i]], c(1,1))) {
    RGs_edges_within <- rep(TRUE, length(RGs))
  } else {
    # Make a vector of intra-infection edges by first creating a block diag. matrix
    mat_within <- Matrix::bdiag(lapply(MOIs[[i]],function(x) matrix(1, ncol=x, nrow=x)))
    colnames(mat_within) <- gs; rownames(mat_within) <- gs
    edges_within <- igraph::as_ids(igraph::E(igraph::graph_from_adjacency_matrix(mat_within, mode = "undirected", diag = F)))
    RGs_edges <- lapply(RGs, function(RG) igraph::as_ids(igraph::E(RG)))
    RGs_edges_within <- sapply(RGs_edges, function(RG_E) all(edges_within %in% RG_E))
  }

  graph_frac_CL_unnormalised <- c(frac_c = sum(RGs_edges_within*RGs_C)/sum(RGs_C),
                                  frac_l = sum(RGs_edges_within*RGs_C_or_L)/length(RGs))

  graph_frac_CL <- graph_frac_CL_unnormalised/sum(graph_frac_CL_unnormalised)

  return(graph_frac_CL)
})

# Compare tabulated results
t(graph_frac_CL)
results

# Compare plotted results
plot(x = graph_frac_CL["frac_c",], y = results[,c("C")],
     pch = as.character(1:length(MOIs)),
     xlim = c(0,1), ylim = c(0,1), bty = "n",
     xlab = c("Relative prior proportion"),
     ylab = c("Posterior probability"))
abline(a = 0, b = 1)

# ------------------------------------------------------------------------------
# Understanding the relapse vs recrudescence results...
#
# The likelihood of the graphs in the MOI c(3,3) case is indeed limited to the
# those with intra-infection relatedness compatible with recrudescence or relapse
# but not reinfection.
# ------------------------------------------------------------------------------
# For the MOIs c(3,3) case
gs <- paste0("g", 1:sum(c(3,3)))
ts <- 1:length(c(3,3))
ts_per_gs <- rep(ts, c(3,3))
gs_per_ts <- split(gs, ts_per_gs)
result33  <- compute_posterior(y, fs, MOIs = c(3,3), return.RG = TRUE, return.logp = T)

RGs <- sapply(result33$RGs, function(RG) RG_to_igraph(RG, gs, ts_per_gs))
CIL_gvn_RGs <- sapply(RGs, compatible_rstrs, gs_per_ts = gs_per_ts) # states
RGs_C_or_L <- sapply(CIL_gvn_RGs, function(CIL_gvn_RG) !("I" %in% CIL_gvn_RG))

# Make a vector of intra-infection edges by first creating a block diag. matrix
mat_within <- Matrix::bdiag(lapply(c(3,3),function(x) matrix(1, ncol=x, nrow=x)))
colnames(mat_within) <- gs; rownames(mat_within) <- gs
edges_within <- igraph::as_ids(igraph::E(igraph::graph_from_adjacency_matrix(mat_within, mode = "undirected", diag = F)))
RGs_edges <- lapply(RGs, function(RG) igraph::as_ids(igraph::E(RG)))
RGs_edges_within <- sapply(RGs_edges, function(RG_E) all(edges_within %in% RG_E))

# extract logp
x <- sapply(result33$RGs, function(RG) RG$logp)
x[x == -Inf] <- NA # Mask -Inf
x <- x - min(x, na.rm = T) # Re-scale before exponentiating (otherwise all 0)
x <- exp(x)/sum(exp(x), na.rm = T) # Exponentiate and normalise

# Inspect probabilities
unique(x[RGs_C_or_L & RGs_edges_within]) # some more likely than others
max(x[!(RGs_C_or_L & RGs_edges_within)]) # all others are unlikely

# Plot the most likely
for(i in which(RGs_C_or_L & RGs_edges_within)){
  plot_RG(RGs[[i]])
}
