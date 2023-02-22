###############################################################################
# Examples of unwanted effect on posterior of prior on graphs
#
# Conclusion: the marginal probability that the first recurrence is a
# reinfection increases as the graph grows only if there are data (example 4);
# otherwise, it decreases (examples 1,2,3). Also, marginal probabilities of the
# different recurrences differ even when they all have the same data (example 4).
#
# Problem: NA seems to negate recrudescence (see first and second recurrence of
# example 2 and 3, respectively). Is there are clause somewhere that compares
# an allele character to NA and negates recrudescence because not equal
#
# To-do list:
# Add an warning for recurrences with no data having non-zero marginal probabilities
# Integrate examples into compute_posterior documentation
################################################################################

# When not in development mode
# devtools::install_github("aimeertaylor/Pv3Rs", build_vignettes = TRUE)
library(Pv3Rs)
fs <- list(m1 = setNames(0.25, "A")) # List of allele frequencies


#===============================================================================
# Example 1: the size of the graph changes and no data whatsoever
#===============================================================================
ys <- list(list(list(m1 = NA), list(m1 = NA)), # 1 recurrence
           list(list(m1 = NA), list(m1 = NA), list(m1 = NA)), # 2 recurrences
           list(list(m1 = NA), list(m1 = NA), list(m1 = NA), list(m1 = NA))) # 3 recurrences
lapply(ys, function(y) compute_posterior(y, fs, return.RG = T)$marg)


#===============================================================================
# Example 2: the size of the graph changes and no recurrent data
#===============================================================================
ys <- list(list(list(m1 = "A"), list(m1 = NA)), # 1 recurrence
           list(list(m1 = "A"), list(m1 = NA), list(m1 = NA)), # 2 recurrences
           list(list(m1 = "A"), list(m1 = NA), list(m1 = NA), list(m1 = NA))) # 3 recurrences
lapply(ys, function(y) compute_posterior(y, fs, return.RG = T)$marg)


#===============================================================================
# Example 3: the size of the graph changes but only one recurrence with data
#===============================================================================
ys <- list(list(list(m1 = "A"), list(m1 = "A")), # 1 recurrence
           list(list(m1 = "A"), list(m1 = "A"), list(m1 = NA)), # 2 recurrences
           list(list(m1 = "A"), list(m1 = "A"), list(m1 = NA), list(m1 = NA))) # 3 recurrences
lapply(ys, function(y) compute_posterior(y, fs, return.RG = T)$marg)


#===============================================================================
# Example 4: the size of the graph changes with repeat data
#===============================================================================
ys <- list(list(list(m1 = "A"), list(m1 = "A")),  # 1 recurrence
           list(list(m1 = "A"), list(m1 = "A"), list(m1 = "A")), # 2 recurrences
           list(list(m1 = "A"), list(m1 = "A"), list(m1 = "A"), list(m1 = "A"))) # 3 recurrences
lapply(ys, function(y) compute_posterior(y, fs)$marg)


