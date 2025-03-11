################################################################################
# Script to understand the effect on bounds of the summation over cliques of
# three or more intra-episode siblings
################################################################################
all_MOIs <- sapply(colnames(maxima), function(x) as.numeric(strsplit(x, split = "")[[1]]))
maxMOIs <- sapply(all_MOIs, function(x) max(x))

# ==============================================================================
# First, sanity check for cases that should not change
# ==============================================================================
# Compare theory with and without for cases that should not change (sanity check)
plot(NULL, xlim = c(0.5, 1), ylim = c(0.5,1),
     main = "Exact bounds for graphs with MOI 2 max",
     ylab = "without summation", xlab = "with summation"); abline(a = 0, b = 1)
points(x = maxima["theory_C_with", maxMOIs < 3],
       y = maxima["theory_C_wout", maxMOIs < 3], pch = 21, bg = "yellow")
points(x = maxima["theory_I_with", maxMOIs < 3],
       y = maxima["theory_I_wout", maxMOIs < 3], pch = 21, bg = "red")
legend("topleft", pt.bg = c("yellow", "red"), pch = 21, bty = "n",
       legend = c("Recrudescence", "Reinfection"))

# Stochastic variation for those that should not change (sanity check)
plot(NULL, xlim = c(0.5, 1), ylim = c(0.5,1),
     main = "Simulated bounds for graphs with MOI 2 max",
     ylab = "without summation", xlab = "with summation"); abline(a = 0, b = 1)
points(x = maxima["sim_C_with", maxMOIs < 3],
       y = maxima["sim_C_wout", maxMOIs < 3], pch = 21, bg = "yellow")
points(x = maxima["sim_I_with", maxMOIs < 3],
       y = maxima["sim_I_wout", maxMOIs < 3], pch = 21, bg = "red")

# Compute maximum difference due to stochastic simulation
max_stoch_diff_sanity <- max(abs(c(maxima["theory_C_with", maxMOIs < 3] - maxima["sim_C_with", maxMOIs < 3],
                                   maxima["theory_I_with", maxMOIs < 3] - maxima["sim_I_with", maxMOIs < 3],
                                   maxima["theory_C_wout", maxMOIs < 3] - maxima["sim_C_wout", maxMOIs < 3],
                                   maxima["theory_I_wout", maxMOIs < 3] - maxima["sim_I_wout", maxMOIs < 3])),
                             na.rm = TRUE) # Delete this line when NAs removed from maxima

# ==============================================================================
# Second, explore results for cases that should change
# ==============================================================================
# First, compare exact with approximate: looks reasonable
plot(NULL, xlim = c(0.5, 1), ylim = c(0.5,1),
     ylab = "simulation", xlab = "theory"); abline(a = 0, b = 1)
points(x = maxima["theory_C_with", maxMOIs > 2], y = maxima["sim_C_with", maxMOIs > 2], pch = 21, bg = "yellow")
points(x = maxima["theory_C_wout", maxMOIs > 2], y = maxima["sim_C_wout", maxMOIs > 2], pch = 22, bg = "yellow")
points(x = maxima["theory_I_with", maxMOIs > 2], y = maxima["sim_I_with", maxMOIs > 2], pch = 21, bg = "red")
points(x = maxima["theory_I_wout", maxMOIs > 2], y = maxima["sim_I_wout", maxMOIs > 2], pch = 22, bg = "red")
legend("topleft", pt.bg = rep(c("yellow", "red"), each = 2), pch = rep(21:22, 2), bty = "n",
       legend = c("Recrudescence with summation", "Recrudescence without summation",
                  "Reinfection with summation", "Reinfection without summation"))

# Second, compare with and without summation
plot(NULL, xlim = c(0.5, 1), ylim = c(0.5,1),
     ylab = "without summation", xlab = "with summation"); abline(a = 0, b = 1)
# points(x = maxima["sim_C_with", maxMOIs > 2], y = maxima["sim_C_wout", maxMOIs > 2], pch = 22, bg = "yellow")
# points(x = maxima["sim_I_with", maxMOIs > 2], y = maxima["sim_I_wout", maxMOIs > 2], pch = 22, bg = "red")
points(x = maxima["theory_C_with", maxMOIs > 2], y = maxima["theory_C_wout", maxMOIs > 2], pch = 21, bg = "yellow")
points(x = maxima["theory_I_with", maxMOIs > 2], y = maxima["theory_I_wout", maxMOIs > 2], pch = 21, bg = "red")
legend("topleft", pt.bg = rep(c("yellow", "red"), 2), pch = rep(21:22, each = 2), bty = "n",
       legend = c("Recrudescence in theory", "Reinfection in theory",
                  "Recrudescence by simulation", "Reinfection by simulation"))


# Is the theoretical difference with and without summation on the order of the
# difference between theory and simulation?
diffs_withwout <- c(maxima["theory_C_with", maxMOIs > 2] - maxima["theory_C_wout", maxMOIs > 2],
                    maxima["theory_I_with", maxMOIs > 2] - maxima["theory_I_wout", maxMOIs > 2])
mean(diffs_withwout < 0, na.rm = T)
max_theory_diff_withwout <- max(abs(diffs_withwout), na.rm = T)
max_stoch_diff_interest <- max(abs(c(maxima["theory_C_with", maxMOIs > 2] - maxima["sim_C_with", maxMOIs > 2],
                                     maxima["theory_I_with", maxMOIs > 2] - maxima["sim_I_with", maxMOIs > 2],
                                     maxima["theory_C_wout", maxMOIs > 2] - maxima["sim_C_wout", maxMOIs > 2],
                                     maxima["theory_I_wout", maxMOIs > 2] - maxima["sim_I_wout", maxMOIs > 2])),
                               na.rm = T)

# Order of difference is comparable to stochastic variation
max_stoch_diff_sanity
max_stoch_diff_interest
max_theory_diff_withwout
