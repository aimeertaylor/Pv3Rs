################################################################################
# The output of compute posterior converges (exactly?) to the prior-weighted
# relative fraction of relationship graphs compatible with the observed data.
#
# This is one when the data suggest there are inter-episode siblings XXX
################################################################################
rm(list = ls())

# Simplex plot parameters
pardefault <- par()
vertex_names <- c(C = "Recrudescence", L = "Relapse", I = "Reinfection")

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


# ==============================================================================
# Marginal probabilities of different recurrences with the same data differ
# slightly
# ==============================================================================
fs <- list(m1 = setNames(c(0.25, 1-0.25), c("A", "Other")))
y <- list(list(m1 = "A"), list(m1 = "A"), list(m1 = "A"), list(m1 = "A"))
compute_posterior(y, fs)$marg


#===============================================================================
# Example for data on a single episode (no recurrent data) using MOIs to change
# graph structure.
#
# A rare homoallelic call has large effect (it limits the summation over
# relationship graphs to those with intra-episode sibling edges). A
# heteroallelic has a minor effect that doesn't change with allele frequency (it
# limits the summation over IBD partitions, but not the summation over
# relationship graphs)
#===============================================================================
# Allele frequencies and MOIs
f_rare <- 0.001
fs = list(m1 = c('A' = f_rare, 'Other' = 1-f_rare))
MOIs <- list(c(1,1), c(2,1), c(3,1), c(2,2), c(3,2), c(3,3))

# Data: a homoallelic or heteroallelic call for a single unpaired episode
y_hom <- list(enroll = list(m1 = c('A')), recur = list(m1 = NA))
y_het <- list(enroll = list(m1 = c('A','Other')), recur = list(m1 = NA))

# Expectation based on graphs for the het case: no limit on graphs
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

# Expectation based on graphs for the hom case: limit on graphs
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

# Compute posterior probabilities and extract marginal probabilities
results_hom <- do.call(rbind, lapply(MOIs, function(x) suppressMessages(compute_posterior(y_hom, fs, MOIs = x)$marg)))
results_het <- do.call(rbind, lapply(MOIs[-1], function(x) suppressMessages(compute_posterior(y_het, fs, MOIs = x)$marg)))

# Project probabilities onto 2D simplex coordinates
xy_hom <- apply(results_hom, 1, project2D)
xy_het <- apply(results_het, 1, project2D)
xy_graph_frac_het <- apply(graph_frac_het, 2, project2D)
xy_graph_frac_hom <- apply(graph_frac_hom, 2, project2D)

# Plotting
par(mfrow = c(1,2))

# Plot 2D simplex:hom
plot_simplex(vertex_names[colnames(results_hom)], classifcation_threshold = 0.5)
points(x = xy_graph_frac_hom["x", ], y = xy_graph_frac_hom["y", ], col = 1:length(MOIs), cex = 2)
points(x = xy_hom["x", ], y = xy_hom["y", ], pch = 20, cex = 1, col = 1:length(MOIs))

# Add legends
legend("bottom", pch = c(20, 1), pt.cex = c(1, 2),
       legend = c("Posterior probability", "Relative weighted graphs"), bty = "n")
legend("left", col = 1:length(MOIs), pch = 20, title = "MOIs", inset = 0,
       legend = sapply(MOIs, paste, collapse = " & "), bty = "n")

# Plot 2D simplex:het
plot_simplex(vertex_names[colnames(results_het)], classifcation_threshold = 0.5)
points(x = xy_graph_frac_het["x", ], y = xy_graph_frac_het["y", ], col = (1:length(MOIs))[-1], cex = 2)
points(x = xy_het["x", ], y = xy_het["y", ], col = (1:length(MOIs))[-1], pch = 20)



#===============================================================================
# Example for data on a single episode (no recurrent data) using recurrences to
# change relationship graph structure.
#
# Adding recurrences has very little impact when data are limited to a single
# episode: no limits are imposed on relationship graphs.
#===============================================================================
# Allele frequencies:
f_rare <- 0.001
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
xy_hom <- apply(first_recur_hom, 2, project2D)
xy_het <- apply(first_recur_het, 2, project2D)

# Plot 2D simplex:het
par(mfrow = c(1,1))
plot_simplex(vertex_names[colnames(results_het[[1]])],  0.5)
points(x = xy_hom["x", ], y = xy_hom["y", ], pch = 20, col = "hotpink")
points(x = xy_het["x", ], y = xy_het["y", ], pch = 20, col = "blue")
legend("topright", inset = 0.1, col = c("hotpink", "blue"), pch = 20, bty = "n",
       legend = c("Homoallelic", "Heteroallelic"))


#===============================================================================
# Example with recurrent data using multiple recurrences to change graph
# structure.
#
# The marginal probability that the first recurrence is a recrudescence
# increases at a decreasing rate as the graph grows without data. Decreasing the
# frequency of the observed allele decreases the probability of reinfection.
# (try changing f_rare).
#===============================================================================
# Allele frequencies:
f_rare <- 0.01
fs <- list(m1 = setNames(c(f_rare, 1-f_rare), c("A", "Other")))

# The number of recurrences increases but only the first recurrence has data
ys <- list(list(list(m1 = "A"), list(m1 = "A")), # 1 recurrence
           list(list(m1 = "A"), list(m1 = "A"), list(m1 = NA)), # 2 recurrences
           list(list(m1 = "A"), list(m1 = "A"), list(m1 = NA), list(m1 = NA)), # etc.
           list(list(m1 = "A"), list(m1 = "A"), list(m1 = NA), list(m1 = NA), list(m1 = NA)))

# Expectation for a rare repeat allele
MOIs <- lapply(ys, determine_MOIs) # First, get MOIs
graph_frac_first <- sapply(1:length(MOIs), FUN = function(i) {
  stuff <- get_graph_stuff(MOIs[[i]])
  # RGs compatible with recrudescence at the first recurrence
  RGs_C_first <- sapply(stuff$CIL_gvn_RGs, function(x) "C" %in% do.call(rbind, strsplit(x, split = ""))[,1])
  # RGs compatible with recrudescence or relapse at the first recurrence
  RGs_C_or_L_first <- sapply(stuff$CIL_gvn_RGs, function(x) !("I" %in% do.call(rbind, strsplit(x, split = ""))[,1]))
  RGs_edges <- lapply(stuff$RGs, function(RG) igraph::as_ids(igraph::E(RG)))
  RGs_g1g2_edge <- sapply(RGs_edges, function(RG_E) all("g1|g2" %in% RG_E)) # A clone between the first and second infection
  graph_frac_unnormalised <- c(frac_C = sum(RGs_g1g2_edge*RGs_C_first)/sum(RGs_C_first),
                               frac_L = sum(RGs_g1g2_edge*RGs_C_or_L_first)/length(stuff$RGs),
                               frac_I = 0)
  graph_frac <- graph_frac_unnormalised/sum(graph_frac_unnormalised)
  return(graph_frac)
})

# Compute posterior probabilities, extract marginal probabilities:
results <- lapply(ys, function(y) compute_posterior(y, fs)$marg)
results_first <- sapply(results, function(result) result[1,]) # Extract 1st recurrence results

# Project onto simplex
xy <- apply(results_first, 2, project2D)
xy <- apply(graph_frac_first, 2, project2D)

# Plot divergence on 2D simplex
plot_simplex(vertex_names, 0.5)
points(x = xy["x", ], y = xy["y", ], pch = "-", col = 1:length(ys), cex = 1)
points(x = xy["x", ], y = xy["y", ], col = 1:length(ys), cex = 1)
legend("right", col = 1:length(ys), pch = "-", pt.cex = 2, bty = "n",
       legend = 1:length(ys), title = "Recurrence \n count")


#===============================================================================
# Example with recurrent data using MOIs to change graph structure.
#===============================================================================
#-------------------------------------------------------------------------------
# Understanding reinfection versus relapse: Can rule out recrudescence using no-matching
# alleles (because errors are not accounted for by the current version of the
# model); however many markers are required to converge to graph expectation.
#-------------------------------------------------------------------------------
# Expectation based on summation over graphs without inter-episode edges
graph_frac_IL <- sapply(1:length(MOIs), FUN = function(i) {
  stuff <- get_graph_stuff(MOIs[[i]])
  RGs_I_or_L <- sapply(stuff$CIL_gvn_RGs, function(RG) !("C" %in% RG))
  RGs_I <- sapply(stuff$CIL_gvn_RGs, function(RG) "I" %in% RG)
  mat_within <- Matrix::bdiag(lapply(MOIs[[i]], function(x) matrix(1, ncol=x, nrow=x)))
  colnames(mat_within) <- stuff$gs; rownames(mat_within) <- stuff$gs
  mat_across <- as.matrix(abs(mat_within-1)) # Convert to across matrix
  edges_across <- igraph::as_ids(igraph::E(igraph::graph_from_adjacency_matrix(mat_across, mode = "undirected", diag = F)))
  RGs_edges <- lapply(stuff$RGs, function(RG) igraph::as_ids(igraph::E(RG)))
  RGs_strangers_across <- sapply(RGs_edges, function(RG_E) !any(edges_across %in% RG_E))
  graph_frac_unnormalised <- c(frac_C = 0,
                               frac_L = sum(RGs_strangers_across*RGs_I_or_L)/length(stuff$RGs),
                               frac_I = sum(RGs_strangers_across*RGs_I)/sum(RGs_I))
  graph_frac <- graph_frac_unnormalised/sum(graph_frac_unnormalised)
  return(graph_frac)
})

# Frequencies, data and MOIs
set.seed(1)
marker_count <- 100
allele_count <- 2
alleles <- LETTERS[1:allele_count]
ms <- paste0("m", 1:marker_count)
fs <- sapply(ms, function(m) {
  setNames(MCMCpack::rdirichlet(1, rep(1, allele_count)), alleles)
}, USE.NAMES = TRUE, simplify = FALSE)
y <- list(enrol = as.list(sapply(ms, function(t) sample(alleles, size = 1, prob = fs[[t]]))),
          recur = as.list(sapply(ms, function(t) sample(alleles, size = 1, prob = fs[[t]]))))
MOIs <- list(c(1,1), c(2,1), c(3,1), c(2,2), c(3,2), c(3,3), c(4,2))
results <- do.call(rbind, lapply(MOIs, function(x) suppressMessages(compute_posterior(y, fs, MOIs = x)$marg)))

# Project posterior probabilities onto 2D simplex coordinates
xy <- apply(results, 1, project2D)
xy_graph_frac <- apply(graph_frac_IL, 2, project2D)

# Plot simplex
par(mfrow = c(1,1))
plot_simplex(vertex_names[colnames(results)], 0.5)
points(x = xy["x", ], y = xy["y", ], col = 1:length(MOIs), pch = 20)
points(x = xy_graph_frac["x", ], y = xy_graph_frac["y", ], col = 1:length(MOIs), pch = 1, cex = 2)

legend("bottom", bty = "n", pch = c(20, 1), pt.cex = c(1,2),
       legend = c("Posterior probability", "Relative graph proportion"), )
legend("left", bty = "n", col = 1:length(MOIs), pch = 20, title = "MOIs",
       legend = sapply(MOIs, paste, collapse = " & "))

#-------------------------------------------------------------------------------
# Recrudescence vs relapse: can rule out reinfection with either many alleles or
# one very rare allele.
#-------------------------------------------------------------------------------
# Frequencies, data and MOIs
f_rare <- 0.001 # Make reinfection unlikely, recrudescence likely
fs = list(m1 = c('1' = f_rare, '2' = 1-f_rare))
y <- list(enrol = list(m1 = "1"), recur = list(m1 = c("1")))
MOIs <- list(c(1,1), c(2,1), c(3,1), c(2,2), c(3,2), c(3,3))

# Expectation based on weighted graphs for rare allele
graph_frac_CL <- sapply(1:length(MOIs), FUN = function(i) {

  stuff <- get_graph_stuff(MOIs[[i]])
  RGs_C_or_L <- sapply(stuff$CIL_gvn_RGs, function(RG) !("I" %in% RG))
  RGs_C <- sapply(stuff$CIL_gvn_RGs, function(RG) "C" %in% RG)

  if(identical(MOIs[[i]], c(1,1))) {
    RGs_edges_within <- rep(TRUE, length(stuff$RGs))
  } else {
    # Make a vector of intra-infection edges by first creating a block diag. matrix
    mat_within <- Matrix::bdiag(lapply(MOIs[[i]],function(x) matrix(1, ncol=x, nrow=x)))
    colnames(mat_within) <- stuff$gs; rownames(mat_within) <- stuff$gs
    edges_within <- igraph::as_ids(igraph::E(igraph::graph_from_adjacency_matrix(mat_within, mode = "undirected", diag = F)))
    RGs_edges <- lapply(stuff$RGs, function(RG) igraph::as_ids(igraph::E(RG)))
    RGs_edges_within <- sapply(RGs_edges, function(RG_E) all(edges_within %in% RG_E))
  }

  graph_frac_unnormalised <- c(frac_C = sum(RGs_edges_within*RGs_C)/sum(RGs_C),
                               frac_L = sum(RGs_edges_within*RGs_C_or_L)/length(stuff$RGs),
                               frac_I = 0)
  graph_frac <- graph_frac_unnormalised/sum(graph_frac_unnormalised)
  return(graph_frac)
})

# Compute posterior probabilities, extract marginal probabilities and project
results <- do.call(rbind, lapply(MOIs, function(x) suppressMessages(compute_posterior(y, fs, MOIs = x)$marg)))

# Project probabilities onto 2D simplex coordinates
xy <- apply(results, 1, project2D)
xy_graph_frac <- apply(graph_frac_CL, 2, project2D)

# Plot 2D simplex
par(mfrow = c(1,1))
plot_simplex(vertex_names[colnames(results)], 0.5)
points(x = xy["x", ], y = xy["y", ], pch = 20, col = 1:length(MOIs))
points(x = xy_graph_frac["x", ], y = xy_graph_frac["y", ], col = 1:length(MOIs), pch = 1, cex = 2)

legend("bottom", bty = "n", pch = c(20, 1), pt.cex = c(1,2),
       legend = c("Posterior probability", "Relative graph proportion"), )
legend("left", bty = "n", col = 1:length(MOIs), pch = 20, title = "MOIs",
       legend = sapply(MOIs, paste, collapse = " & "))

# ------------------------------------------------------------------------------
# Understanding the relapse vs recrudescence results extra....
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







