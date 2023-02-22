# Marginals are not the same give same data
# When not in development mode
# devtools::install_github("aimeertaylor/Pv3Rs", build_vignettes = TRUE)
library(Pv3Rs)
library(gtools) # For rdirichlet
set.seed(1) # Given rdirichlet


# Set up
causes <- c("C", "L", "I") # Possible causes of recurrence
alleles <- c("A", "T", "C", "G") # Possible alleles
return_freqs <- function(num_markers, alleles) {
  n.as <- length(alleles) # Number of possible alleles
  alphas <- rep(1, n.as) # Dirichlet param. vector
  setNames(lapply(1:num_markers, function(x) setNames(rdirichlet(1, alphas)[1,], alleles)),
           paste0("m", 1:num_markers)) # List of allele frequencies
}

#===============================================================================
# Example: the size of the graph changes the marginal posterior of the first recurrence
# Format problem in first instance
# in second
#===============================================================================
# List of allele frequencies
fs <- return_freqs(num_markers = 1, alleles)

# Data for multiple scenarios: different numbers of recurrences but data only on the first
ys <- list(list(list(m1 = "T"), list(m1 = "T")),  # Single recurrence
           list(list(m1 = "T"), list(m1 = "T"), list(m1 = NA)), # Two recurrences
           list(list(m1 = "T"), list(m1 = "T"), list(m1 = NA), list(m1 = NA))) # Three recurrences

# Computation
lapply(ys, function(y) {post <- compute_posterior(y, fs); post$marg})

# Data for multiple scenarios: different numbers of recurrences but data only on the first
ys <- list(list(m1 = "T", m1 = "T"),  # Single recurrence
           list(m1 = "T", m1 = "T", m1 = NA), # Two recurrences
           list(m1 = "T", m1 = "T", m1 = NA, m1 = NA))  # Three recurrences

# Computation
lapply(ys, function(y) {post <- compute_posterior(y, fs); post$marg})


#===============================================================================
# Example showing how the marginals are not the same given the same data
#===============================================================================
# Data for a single scenario: four recurrences, each with the same data
y <- list(list(m1 = "T", m2 = "T", m3 = "C"), #0 initial infection
          list(m1 = "T", m2 = "T", m3 = "C"), #1 first recurrence
          list(m1 = "T", m2 = "T", m3 = "C"), #2 second recurrence
          list(m1 = "T", m2 = "T", m3 = "C"), #3 third recurrence
          list(m1 = "T", m2 = "T", m3 = "C")) #4 fourth recurrence

# List of allele frequencies
fs <- return_freqs(num_markers = length(y[[1]]), alleles)

# Computation
post <- compute_posterior(y, fs)
post$marg

