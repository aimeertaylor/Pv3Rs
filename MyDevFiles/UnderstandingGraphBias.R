################################################################################
# To-do list:
# Examples of unwanted effects on posterior of prior on graphs
# Make this into a vignette
# Add an warning for recurrences with no data can have non-zero marginal probabilities
# Add examples from Thur 22nd August
################################################################################

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
Obs_allele <- 0.25
fs <- list(m1 = setNames(c(Obs_allele, 1-Obs_allele), c("A", "Other")))

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
plot_simplex(v_labels = rownames(first_recur0))
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


# ------------------------------------------------------------------------------
# An alternative visualisation of the same phenomena:
# ------------------------------------------------------------------------------
par(mar = pardefault$mar)
plot(NULL, ylim = c(0,1), xlim = c(1,4), xaxt = "n", las = 1, bty = "n",
     xlab = "Number of recurrences", ylab = "Posterior probabilities for the first recurrence")
axis(side = 1, at  = 1:4)
abline(h = first_recur0[,1], lty = "dotted")
text(x = rep(1,3), y = first_recur0[,1], labels = rownames(first_recur0), pos = rep(3,1))
legend("topright", inset = 0.1, col = c("red", "blue"), pch = 20, bty = "n",
       legend = c("Only the first recurrence has data",
                  "Recurrences have repeat data"))

for(i in 1:3) lines(x = 1:4, y = first_recur0[i,], type = "b", col = "red", pch = 19)
for(i in 1:3) lines(x = 1:4, y = first_recur1[i,], type = "b", col = "blue", pch = 17)


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
# One very rare allele is as informative as many alleles (both force the
# posterior away from reinfection), however many alleles slower (example
# removed). Having ruled out reinfection, graph structure influences posterior
# odds of relapse versus recrudescence.
#============================================================================
# Markers, alleles and allele frequencies:
f_rare <- 0.0001
fs = list(m1 = c('1' = f_rare, '2' = 1-f_rare))
y <- list(enrol = list(m1 = "1"), recur = list(m1 = "1"))

MOIs <- list(c(1,1), c(2,1), c(3,1), c(2,2), c(3,2), c(3,3))

# Compute posterior probabilities and extract marginal probabilities:
results <- do.call(rbind, lapply(MOIs, function(x) suppressMessages(compute_posterior(y, fs, MOIs = x)$marg)))
xy <- apply(results, 1, project2D) # Project probabilities onto 2D simplex coordinates

# Plot 2D simplex
pardefault <- par()
par(mar = c(0,0,0,0))
vertex_names <- c(C = "Recrudescence", L = "Relapse", I = "Reinfection")
plot_simplex(v_labels = vertex_names[colnames(results)], classifcation_threshold = 0.5)
points(x = xy["x", ], y = xy["y", ], pch = as.character(1:length(MOIs)), cex = 0.25)

# Restore plotting margins
par(mar = pardefault$mar)

#============================================================================
# Examples without recurrent data
#
# The prior is only returned when the initial infection has a homoallelic call
# and an MOI of one and.
#
# A homoallelic call has a systematic effect that is amplified when the allele
# of the homoallelic call is rare — what is going on?
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
