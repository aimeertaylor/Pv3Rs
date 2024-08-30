################################################################################
# Does it make sense to have both return.RG and return.logp?.
# Examples of graph prior effects on posterior
# Make this into a vignette: need to compute RG after ruling out a state
# If the data can only rule out a state, then the result reflects the prior on graphs given MOIs.
# Add an warning for recurrences with no data can have non-zero marginal probabilities
# What is going on with the homoallelic exampled with no recurrent data?
################################################################################
#============================================================================
# Examples without recurrent data
#
# The prior is only returned when the initial infection has a homoallelic call
# and an MOI of one and.
#
# A homoallelic call has a systematic effect that is amplified when the allele
# of the homoallelic call is rare.
#
# A heteroallelic call has a minor effect that doesn't change with the allele
# frequency of the rare allele.
#============================================================================
f_rare <- 0.001
fs = list(m1 = c('1' = f_rare, '2' = 1-f_rare))
MOIs <- list(c(1,1), c(2,1), c(1,3), c(3,1), c(2,2), c(3,2), c(3,3))

# Data
y_hom <- list(enroll = list(m1 = c('1')), recur = list(m1 = NA))
y_het <- list(enroll = list(m1 = c('1','2')), recur = list(m1 = NA))

# Compute posterior probabilities and extract marginal probabilities:
results_hom <- do.call(rbind, lapply(MOIs, function(x) suppressMessages(compute_posterior(y_hom, fs, MOIs = x)$marg)))
results_het <- do.call(rbind, lapply(MOIs[-1], function(x) suppressMessages(compute_posterior(y_het, fs, MOIs = x)$marg)))

# Print results
results_hom # Prior is only returned for MOI = c(1,1)
results_het

# Project probabilities onto 2D simplex coordinates
xy_hom <- apply(results_hom, 1, project2D)
xy_het <- apply(results_het, 1, project2D)

pardefault <- par()
par(mfrow = c(1,2))

# Plot 2D simplex:hom
par(mar = c(0,0,0,0))
vertex_names <- c(C = "Recrudescence", L = "Relapse", I = "Reinfection")
plot_simplex(v_labels = vertex_names[colnames(results)], classifcation_threshold = 0.5)
points(x = xy_hom["x", ], y = xy_hom["y", ], pch = as.character(1:length(MOIs)), cex = 0.25)
title(main = "Homologous initial", line = -2)

# Plot 2D simplex:het
par(mar = c(0,0,0,0))
vertex_names <- c(C = "Recrudescence", L = "Relapse", I = "Reinfection")
plot_simplex(v_labels = vertex_names[colnames(results)], classifcation_threshold = 0.5)
points(x = xy_het["x", ], y = xy_het["y", ], pch = as.character(1:length(MOIs)), cex = 0.25)
title(main = "Heterologous initial", line = -2)

# Restore plotting margins
par(mar = pardefault$mar)





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
ys0 <- list(list(list(m1 = "A"), list(m1 = "A")), # 1 recurrence
            list(list(m1 = "A"), list(m1 = "A"), list(m1 = NA)), # 2 recurrences
            list(list(m1 = "A"), list(m1 = "A"), list(m1 = NA), list(m1 = NA)), # etc.
            list(list(m1 = "A"), list(m1 = "A"), list(m1 = NA), list(m1 = NA), list(m1 = NA)))

# The number of recurrences increases and all have the same data
ys1 <- list(list(list(m1 = "A"), list(m1 = "A")),  # 1 recurrence
            list(list(m1 = "A"), list(m1 = "A"), list(m1 = "A")), # 2 recurrences
            list(list(m1 = "A"), list(m1 = "A"), list(m1 = "A"), list(m1 = "A")), # etc.
            list(list(m1 = "A"), list(m1 = "A"), list(m1 = "A"), list(m1 = "A"), list(m1 = "A")))

# Compute posterior probabilities and extract marginal probabilities:
results0 <- lapply(ys0, function(y) compute_posterior(y, fs)$marg)
results1 <- lapply(ys1, function(y) compute_posterior(y, fs)$marg)

# Extract results for the first recurrence only:
first_recur0 <- sapply(results0, function(result) result[1,])
first_recur1 <- sapply(results1, function(result) result[1,])

# As expected, posterior probabilities agree when there is only one recurrence
# (first columns match), otherwise they diverge
first_recur0
first_recur1


# ------------------------------------------------------------------------------
# Plot divergence on 2D simplex
# ------------------------------------------------------------------------------
n_recur <- ncol(first_recur0)
pardefault <- par()
par(mar = c(0,0,0,0))
plot_simplex(v_labels = rownames(first_recur0), classifcation_threshold = 0.5)
legend("topright", inset = 0.1, col = c("red", "blue"), pch = 20, bty = "n",
       legend = c("Only the first recurrence has data",
                  "Recurrences have repeat data"))

# Project and plot first_recur0
xy <- apply(first_recur0, 2, project2D)
arrows(x0 = xy["x", 1], x1 = xy["x", n_recur],
       y0 = xy["y", 1], y1 = xy["y", n_recur],
       length = 0.05, col = "red")
points(x = xy["x", ], y = xy["y", ], pch = ".")

# Project and plot first_recur1
xy <- apply(first_recur1, 2, project2D)
arrows(x0 = xy["x", 1], x1 = xy["x", n_recur],
       y0 = xy["y", 1], y1 = xy["y", n_recur],
       length = 0.05, col = "blue")
points(x = xy["x", ], y = xy["y", ], pch = ".")


# The marginal probabilities of recurrences with no data diverge from the prior
# at a decreasing rate, at least for relapse and reinfection:
results0[[4]] - array(1/3, dim = dim(results0[[4]]))

# The marginal probabilities of different recurrences differ when they all
# recurrences have the same data:
results1[[4]] - matrix(results1[[4]][1,], byrow = T,
                       nrow = nrow(results0[[4]]),
                       ncol = ncol(results0[[4]]))

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
f_rare <- 0.0001 # Make reinfection unlikely
fs = list(m1 = c('1' = f_rare, '2' = 1-f_rare))
y <- list(enrol = list(m1 = "1"), recur = list(m1 = c("1")))
MOIs <- list(c(1,1), c(2,1), c(3,1), c(2,2), c(3,2), c(3,3))

# Compute posterior probabilities, extract marginal probabilities and project
results <- do.call(rbind, lapply(MOIs, function(x) suppressMessages(compute_posterior(y, fs, MOIs = x)$marg)))
xy <- apply(results, 1, project2D) # Project probabilities onto 2D simplex coordinates

# Plot 2D simplex
pardefault <- par()
par(mar = c(0,0,0,0))
vertex_names <- c(C = "Recrudescence", L = "Relapse", I = "Reinfection")
plot_simplex(v_labels = vertex_names[colnames(results)], classifcation_threshold = 0.5)
points(x = xy["x", ], y = xy["y", ], pch = as.character(1:length(MOIs)), cex = 0.25)
par(mar = pardefault$mar) # Restore plotting margins

# ------------------------------------------------------------------------------
# Understanding the relapse vs recrudescence results
#
# The fraction of intra-infection related graphs compatible with
# recrudescence versus relapse is a function of the MOI
# ------------------------------------------------------------------------------
graph_frac_CL <- sapply(2:length(MOIs), FUN = function(i) {

  gs <- paste0("g", 1:sum(MOIs[[i]]))
  ts <- 1:length(MOIs[[i]])
  ts_per_gs <- rep(ts, MOIs[[i]])
  gs_per_ts <- split(gs, ts_per_gs)
  RGs <- enumerate_RGs(MOIs[[i]])
  CIL_gvn_RGs <- sapply(RGs, compatible_rstrs, gs_per_ts = gs_per_ts) # states
  RGs_C_or_L <- sapply(CIL_gvn_RGs, function(RG) !("I" %in% RG))
  RGs_C <- sapply(CIL_gvn_RGs, function(RG) "C" %in% RG)

  # Make a vector of intra-infection edges by first creating a block diag. matrix
  mat_within <- Matrix::bdiag(lapply(MOIs[[i]],function(x) matrix(1, ncol=x, nrow=x)))
  colnames(mat_within) <- gs; rownames(mat_within) <- gs
  edges_within <- igraph::as_ids(igraph::E(igraph::graph_from_adjacency_matrix(mat_within, mode = "undirected", diag = F)))
  RGs_edges_within <- sapply(RGs, function(RG) all(edges_within %in% igraph::as_ids(igraph::E(RG))))

  graph_frac_CL_unnormalised <- c(frac_c = sum(RGs_edges_within*RGs_C)/sum(RGs_C),
                                  frac_l = sum(RGs_edges_within*RGs_C_or_L)/length(RGs))

  graph_frac_CL <- graph_frac_CL_unnormalised/sum(graph_frac_CL_unnormalised)

  return(graph_frac_CL)
})

# Compare tabulated results
t(graph_frac_CL)
results

# Compare plotted results
plot(x = graph_frac_CL["frac_c",], y = results[-1,c("C")],
     xlim = c(0,1), ylim = c(0,1), pch = 20, bty = "n",
     xlab = c("Prior prevalence"),
     ylab = c("Posterior probability"))
abline(a = 0, b = 1)

# TO SORT FROM HERE:
# For the MOIs c(3,3) case
gs <- paste0("g", 1:6)
ts <- 1:length(c(3,3))
ts_per_gs <- rep(ts, c(3,3))
X <- compute_posterior(y, fs, MOIs = c(3,3), return.RG = TRUE, return.logp = T)

# extract logp
x <- sapply(X$RGs, function(RG) RG$logp)
x[x == -Inf] <- NA # Mask -Inf
x <- x - min(x, na.rm = T) # Re-scale before exponentiating (otherwise all 0)
x <- exp(x)/sum(exp(x), na.rm = T) # Exponentiate and normalise

unique(x[as.logical(RGs_C_or_L * RGs_edges_within)]) # some more likely than others
max(x[!(RGs_C_or_L * RGs_edges_within)]) # all others are unlikely
