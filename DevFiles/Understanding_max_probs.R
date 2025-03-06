################################################################################
# Script to explore: 1) how maximum probabilities vary with MOI vectors; and 2)
# effect of summation over cliques of 3 or more intra-episode siblings To-do:
# finish 2) and Take-away 2): including / excluding sibling cliques makes little
# difference
# To-do: elaborate on all ones graph
# Change: point characters on other graphs
################################################################################
rm(list = ls())
load("../data/maxima.rda")
all_MOIs <- sapply(colnames(max_probs), function(x) as.numeric(strsplit(x, split = "")[[1]]))
maxMOIs <- sapply(all_MOIs, function(x) max(x))

# ==============================================================================
# General understanding of maximum probabilities
# ==============================================================================
nMOI <- sapply(all_MOIs, function(x) length(x)) # Number of recurrences
MOI1 <- sapply(all_MOIs, function(x) x[1]) # MOI of episode one
MOI2 <- sapply(all_MOIs, function(x) x[2]) # MOI of episode two
MOI12diff <- sapply(all_MOIs, function(x) x[1] - x[2]) # Difference
MOI2 <- sapply(all_MOIs, function(x) x[2]) # MOI of episode two
Episode_counts <- sort(unique(nMOI)) #
MOI1s_log <- sapply(all_MOIs, function(x) all(x == 1)) # MOI of episode two

# Imbalance between episodes one and two
plot(x = MOI12diff, y = max_probs["theory_C_with", ], ylim = c(0.5, 1),
     col = MOI1, pch = as.character(nMOI), main = "Recrudescence",
     xlab = "MOI episode 1 - MOI episode 2", ylab = "Theoretical maximum probability")
text(y = 0.5, x = 0, adj = 1,
     labels = paste0("Number of episodes: ", paste(Episode_counts, collapse = " ")))
legend("left", title = "MOI of episode 1", bty = "n",
       legend = unique(MOI1), fill = unique(MOI1))
plot(x = MOI12diff, y = max_probs["theory_I_with", ], ylim = c(0.5, 1),
     col = MOI1, pch = as.character(nMOI), main = "Reinfection",
     xlab = "MOI episode 1 - MOI episode 2", ylab = "Theoretical maximum probability")
text(y = 0.5, x = 0, adj = 1,
     labels = paste0("Number of episodes: ", paste(Episode_counts, collapse = " ")))
legend("right", title = "MOI of episode 1", bty = "n",
       legend = unique(MOI1), fill = unique(MOI1))

# Focusing in on cases where the MOI is always one
plot(x = nMOI[MOI1s_log], y = max_probs["theory_C_with", MOI1s_log],
     type = "b", ylim = c(0.5, 1), main = "", bg = "yellow", pch = 21,
     xlab = "Number of episodes", ylab = "Theoretical maximum probability")
lines(x = nMOI[MOI1s_log], y = max_probs["theory_I_with", MOI1s_log],
      type = "b", bg = "red", pch = 21)


# Focusing in on cases where MOI episode 1 - MOI episode a2 = 0
# Trend is the same for both recrudescence and reinfection
# (number of potentially independent inter-episode edges for episodes 1 & 2)
plot(x = MOI2[MOI12diff == 0], bty = "n",
     y = max_probs["theory_C_with", ][MOI12diff == 0],
     ylim = c(0.5, 1), main = "Recrudescence",
     col = MOI1[MOI12diff == 0], pch = nMOI[MOI12diff == 0],
     xlab = "MOI of episodes 1 = 2", ylab = "Theoretical maximum probability")
legend("bottomright", title = "Number of episodes", bty = "n",
       legend = unique(nMOI[MOI12diff == 0]),
       pch = unique(nMOI[MOI12diff == 0]))
legend("bottom", title = "MOI of episode 1", bty = "n",
       legend = unique(nMOI[MOI12diff == 0]),
       fill = unique(nMOI[MOI12diff == 0]))
plot(x = MOI2[MOI12diff == 0], bty = "n",
     y = max_probs["theory_I_with", ][MOI12diff == 0],
     ylim = c(0.5, 1), main = "Reinfection",
     col = MOI1[MOI12diff == 0], pch = nMOI[MOI12diff == 0],
     xlab = "MOI of episodes 1 = 2", ylab = "Theoretical maximum probability")
legend("bottomright", title = "Number of episodes", bty = "n",
       legend = unique(nMOI[MOI12diff == 0]),
       pch = unique(nMOI[MOI12diff == 0]))
legend("bottom", title = "MOI of episode 1", bty = "n",
       legend = unique(nMOI[MOI12diff == 0]),
       fill = unique(nMOI[MOI12diff == 0]))

# What is the smallest maximum probability? 0.57
ind0 <- max_probs[grep("theory", rownames(max_probs)), ] == 0
min_max_prob <- min(max_probs[grep("theory", rownames(max_probs)), ][!ind0])

# Which MOI vector has the smallest maximum probability?
min_max_prob_ind <- which(max_probs == min_max_prob, arr.ind = T)
rownames(max_probs)[min_max_prob_ind[,"row"]]
colnames(max_probs)[min_max_prob_ind[,"col"]] # 1,1 case with 3 episodes


# ==============================================================================
# Summation over cliques of 3 or more intra-episode siblings
# ==============================================================================

# ------------------------------------------------------------------------------
# First, sanity check for cases that should not change
# ------------------------------------------------------------------------------
# Compare theory with and without for cases that should not change (sanity check)
plot(NULL, xlim = c(0.5, 1), ylim = c(0.5,1), main = "Theory for graphs with MOI 2 max",
     ylab = "without summation", xlab = "with summation"); abline(a = 0, b = 1)
points(x = max_probs["theory_C_with", maxMOIs < 3], y = max_probs["theory_C_wout", maxMOIs < 3], pch = 21, bg = "yellow")
points(x = max_probs["theory_I_with", maxMOIs < 3], y = max_probs["theory_I_wout", maxMOIs < 3], pch = 21, bg = "red")

# Stochastic variation for those that should not change (sanity check)
plot(NULL, xlim = c(0.5, 1), ylim = c(0.5,1), main = "Simulation for graphs with MOI 2 max",
     ylab = "without summation", xlab = "with summation"); abline(a = 0, b = 1)
points(x = max_probs["sim_C_with", maxMOIs < 3], y = max_probs["sim_C_wout", maxMOIs < 3], pch = 21, bg = "yellow")
points(x = max_probs["sim_I_with", maxMOIs < 3], y = max_probs["sim_I_wout", maxMOIs < 3], pch = 21, bg = "red")

# Compute maximum difference due to stochastic simulation
max_stoch_diff_sanity <- max(abs(c(max_probs["theory_C_with", maxMOIs < 3] - max_probs["sim_C_with", maxMOIs < 3],
                                   max_probs["theory_I_with", maxMOIs < 3] - max_probs["sim_I_with", maxMOIs < 3],
                                   max_probs["theory_C_wout", maxMOIs < 3] - max_probs["sim_C_wout", maxMOIs < 3],
                                   max_probs["theory_I_wout", maxMOIs < 3] - max_probs["sim_I_wout", maxMOIs < 3])))

#-------------------------------------------------------------------------------
# Second, explore results for cases that should change
#-------------------------------------------------------------------------------
# First, compare theory with simulation: looks reasonable
plot(NULL, xlim = c(0.5, 1), ylim = c(0.5,1),
     ylab = "simulation", xlab = "theory"); abline(a = 0, b = 1)
points(x = max_probs["theory_C_with", maxMOIs > 2], y = max_probs["sim_C_with", maxMOIs > 2], pch = 21, bg = "yellow")
points(x = max_probs["theory_C_wout", maxMOIs > 2], y = max_probs["sim_C_wout", maxMOIs > 2], pch = 22, bg = "yellow")
points(x = max_probs["theory_I_with", maxMOIs > 2], y = max_probs["sim_I_with", maxMOIs > 2], pch = 21, bg = "red")
points(x = max_probs["theory_I_wout", maxMOIs > 2], y = max_probs["sim_I_wout", maxMOIs > 2], pch = 22, bg = "red")
legend("topleft", pt.bg = rep(c("yellow", "red"), each = 2), pch = rep(21:22, 2), bty = "n",
       legend = c("Recrudescence with summation",
                  "Recrudescence without summation",
                  "Reinfection with summation",
                  "Reinfection without summation"))

# Second, compare with and without summation
plot(NULL, xlim = c(0.5, 1), ylim = c(0.5,1),
     ylab = "without summation", xlab = "with summation"); abline(a = 0, b = 1)
points(x = max_probs["sim_C_with", maxMOIs > 2], y = max_probs["sim_C_wout", maxMOIs > 2], pch = 22, bg = "yellow")
points(x = max_probs["sim_I_with", maxMOIs > 2], y = max_probs["sim_I_wout", maxMOIs > 2], pch = 22, bg = "red")
points(x = max_probs["theory_C_with", maxMOIs > 2], y = max_probs["theory_C_wout", maxMOIs > 2], pch = 21, bg = "yellow")
points(x = max_probs["theory_I_with", maxMOIs > 2], y = max_probs["theory_I_wout", maxMOIs > 2], pch = 21, bg = "red")
legend("topleft", pt.bg = rep(c("yellow", "red"), 2), pch = rep(21:22, each = 2), bty = "n",
       legend = c("Recrudescence in theory",
                  "Reinfection in theory",
                  "Recrudescence by simulation",
                  "Reinfection by simulation"))


# Is the theoretical difference with and without summation on the order of the
# difference between theory and simulation?
diffs_withwout <- c(max_probs["theory_C_with", maxMOIs > 2] - max_probs["theory_C_wout", maxMOIs > 2],
                    max_probs["theory_I_with", maxMOIs > 2] - max_probs["theory_I_wout", maxMOIs > 2])
mean(diffs_withwout < 0)
max_theory_diff_withwout <- max(abs(diffs_withwout))

max_stoch_diff_interest <- max(abs(c(max_probs["theory_C_with", maxMOIs > 2] - max_probs["sim_C_with", maxMOIs > 2],
                                     max_probs["theory_I_with", maxMOIs > 2] - max_probs["sim_I_with", maxMOIs > 2],
                                     max_probs["theory_C_wout", maxMOIs > 2] - max_probs["sim_C_wout", maxMOIs > 2],
                                     max_probs["theory_I_wout", maxMOIs > 2] - max_probs["sim_I_wout", maxMOIs > 2])))

# Order of difference is comparable to stochastic variation
max_stoch_diff_sanity
max_stoch_diff_interest
max_theory_diff_withwout




