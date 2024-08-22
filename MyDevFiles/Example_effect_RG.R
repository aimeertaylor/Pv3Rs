#============================================================================
# Example of the effect on the posterior of the prior on relationship graphs
# graphs using MOIs
#============================================================================
# Allele frequencies:
fs <- list(m1 = setNames(c(0.25, 1-0.25), c("A", "Other")))

# Data (toggle between recur being "A" and NA)
y <- list(enroll = list(m1 = "A"), recur = list(m1 = "A"))

MOIs <- list(c(1,1), c(2,1), c(3,1), c(2,2), c(3,2), c(3,3))

# Compute posterior probabilities and extract marginal probabilities:
results <- do.call(rbind, lapply(MOIs, function(x) compute_posterior(y, fs, MOIs = x)$marg))

# Plot 2D simplex
n_recur <- max(sapply(ys, length)-1)
pardefault <- par()
par(mar = c(0,0,0,0))
vertex_names <- c(C = "Recrudescence", L = "Relapse", I = "Reinfection")
plot_simplex(v_labels = vertex_names[colnames(results)], classifcation_threshold = 0.5)

# Project probabilities onto 2D simplex coordinates
xy <- apply(results, 1, project2D)

# # Plot divergence from one recurrence to four:
# arrows(x0 = xy["x", 1], x1 = xy["x", length(MOIs)],
#        y0 = xy["y", 1], y1 = xy["y", length(MOIs)],
#        length = 0.05, col = "red")

# Plot a point for each recurrence from one to four:
points(x = xy["x", ], y = xy["y", ], pch = as.character(1:length(MOIs)), cex = 0.2)

# Restore plotting margins
par(mar = pardefault$mar)
