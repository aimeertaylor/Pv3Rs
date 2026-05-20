################################################################################
# Generate results for a initial infection of two siblings and
# a recurrent sibling, where siblings are half siblings.
#
# For consistency with theoretical results, when siblings are half siblings, the
# first marker is forced to have three different alleles.
#
# Explore two scenarios. One where parents draw from the same
# allele distribution. Another with admixture where the first parent who parents
# intra-episode parasites draws alleles disproportionally:
#
# Half siblings:
# initial <- rbind(child12, child13) # migrant is 1
# relapse <- rbind(child23)
#
# Generate results for data on all marker counts when alleles are
# equifrequent, and for a subset of marker counts when alleles are not.
################################################################################
rm(list = ls())
library(Pv3Rs)

#===============================================================================
# Question one: what is the impact of the prior?
#===============================================================================
marker_count <- 4 # Number of markers
ms <- paste0("m", 1:marker_count) # Marker names
all_As <- sapply(ms, function(t) "A", simplify = F) # As for all markers
all_Bs <- sapply(ms, function(t) "B", simplify = F) # Bs for all markers
BAAAAA <- c(m1 = "B", sapply(ms[-1], function(t) "A", simplify = F)) # Bs for all markers

fA <- 0.1 # Frequency of rare allele
fB <- 1 - fA # Frequency of common allele
fs <- sapply(ms, function(m) c("A"=fA, "B"=fB), simplify = FALSE)

hom <- list(enrol = all_As, recur = BAAAAA)
het <- list(enrol = all_As, recur = all_Bs)
prior_low_I <- array(c(0, 0.99, 0.01), dim = c(1,3), dimnames = list(NULL, c("C", "L", "I")))
prior_low_L <- array(c(0, 0.01, 0.99), dim = c(1,3), dimnames = list(NULL, c("C", "L", "I")))

# Prior has zero impact when data are strongly suggestive of relapse
suppressMessages(compute_posterior(hom, fs))$marg
suppressMessages(compute_posterior(hom, fs, prior_low_L))$marg

# Prior does impact when data are consistent with reinfection
suppressMessages(compute_posterior(het, fs))$marg
suppressMessages(compute_posterior(het, fs, prior_low_I))$marg

#===============================================================================
# Question two: can the prior offset misspecification?
#===============================================================================

#===============================================================================
# Magic numbers / quantities
#===============================================================================
n_repeats <- 10 # Number of simulations per parameter combination
max_n_markers <- 150 # 150 # Number of markers for which RG likelihood returned
c_cutoff <- 99 # Switch from Dirichlet r.v. to 1/n_alleles above c_cutoff
seed <- 10 # For reproducibility
prior <- array(c(0.05, 0.9, 0.05), dim = c(1,3), dimnames = list(NULL, c("C", "L", "I")))

#===============================================================================
# Stores for results
#===============================================================================
ps_store <- list() # p for posterior
ps_store_all_ms_uniform <- list() # all_ms for all marker counts
ps_store_all_ms_admix_rare <- list() # for admixed case with some rare freqs

#===============================================================================
# Load frequencies and data
#===============================================================================
load("output_HalfSib.PCLikeSib.rda")
fs_store <- output_HalfSib.PCLikeSib$Half$fs_store
ys_store <- output_HalfSib.PCLikeSib$Half$ys_store
m_rorder <- output_HalfSib.PCLikeSib$Half$m_rorder

#=============================================================================
# Generate results for markers 1:max_n_markers
#=============================================================================
c <- 100 # For uniform allele frequencies and all markers
admixture <- FALSE # For parents from the same population
fs <- fs_store[[as.character(c)]] # Extract frequencies
for(i in 1:n_repeats){
  y_all_markers <- ys_store[[as.character(c)]][[sprintf("admixture_%s", admixture)]][[i]]
  for(j in 1:max_n_markers){
    marker_subset <- m_rorder[1:j]
    y <- sapply(y_all_markers, function(x) x[marker_subset], simplify = FALSE)
    ps <- suppressMessages(compute_posterior(y, fs, prior))
    ps_store_all_ms_uniform[[as.character(i)]][[j]] <- ps$marg
  }
}

c <- 0.5 # For admixed case with imbalanced allele frequencies
admixture <- TRUE # For parents from different populations
fs <- fs_store[[as.character(c)]] # Extract frequencies
for(i in 1:n_repeats){
  y_all_markers <- ys_store[[as.character(c)]][[sprintf("admixture_%s", admixture)]][[i]]
  for(j in 1:max_n_markers){
    marker_subset <- m_rorder[1:j]
    y <- sapply(y_all_markers, function(x) x[marker_subset], simplify = FALSE)
    ps <- suppressMessages(compute_posterior(y, fs, prior))
    ps_store_all_ms_admix_rare[[as.character(i)]][[j]] <- ps$marg
  }
}

#=============================================================================
# Bundle magic numbers, data and results
#=============================================================================
output <- list(
  n_repeats = n_repeats,
  max_n_markers = max_n_markers,
  m_rorder = m_rorder,
  fs_store = fs_store,
  c_cutoff = c_cutoff,
  seed = seed,
  prior = prior,
  ps_store = ps_store,
  ps_store_all_ms_uniform = ps_store_all_ms_uniform,
  ps_store_all_ms_admix_rare = ps_store_all_ms_admix_rare)

save(output, file = "output_HalfSib_relapse_prior.rda")


#=============================================================================
# Results
#=============================================================================
cols <- RColorBrewer::brewer.pal(n = output$n_repeats, "Paired") # For repeats
par(mfrow = c(2,1))

# Compute effective cardinality cumulatively in the order of markers genotyped
cum_card_eff <- sapply(output$fs_store, function(fs) {
  cumsum(sapply(fs[output$m_rorder], function(x) 1/sum(x^2)))})

# Determine which fs in cum_card_eff are equifrequent
equifs <- which(as.numeric(colnames(cum_card_eff)) > output$c_cutoff)

# ------------------------------------------------------------------------------
# Plot the posterior relapse probability trajectories
# ------------------------------------------------------------------------------
plot(NULL, bty = "n", las = 1, xaxt = "n", xlim = c(1,max_n_markers),
     xlab = "Marker count (effective cardinality)",
     ylim = c(0,1), ylab = "Posterior relapse probability")
legend("right", lwd = 2, col = cols, legend = 1:output$n_repeats, bty = "n",
       cex = 0.5, title = "Repeat")

# Add horizontal axis
axis_at <- c(1, seq(0, max_n_markers, 50)[-1])
axis(side = 1, at = axis_at, cex.axis = 0.7) # Marker count
axis(side = 1, line = 1, at = axis_at, cex.axis = 0.6, # Effective cardinality
     tick = F, labels = sprintf("(%s)", cum_card_eff[,equifs][axis_at]))

# Add trajectories
for(i in 1:output$n_repeats){
  lines(x = 1:max_n_markers,
        y = sapply(output$ps_store_all_ms_uniform[[as.character(i)]],
                   function(x) x[,"L"]),
        col = cols[i], lwd = 2)
}

# ------------------------------------------------------------------------------
# Plot the posterior relapse probability trajectories
# ------------------------------------------------------------------------------
plot(NULL, xlim = c(1,max_n_markers), ylim = c(0,1), bty = "n", las = 1,
     xaxt = "n", xlab = "Marker count (effective cardinality)",
     ylab = "Posterior relapse probability")

legend("right", col = cols, lwd = 2, inset = 0.1, legend = (1:output$n_repeats)-1,
       bty = "n", cex = 0.5, title = "Repeats")

# Add horizontal axis
axis(side = 1, at = axis_at, cex.axis = 0.7) # Marker count
axis(side = 1, line = 1, at = axis_at, cex.axis = 0.6, # Effective cardinality
     tick = F, labels = sprintf("(%s)", round(cum_card_eff[,"0.5"][axis_at])))

# Add trajectories
for(i in 1:output$n_repeats){
  y <- sapply(output$ps_store_all_ms_admix_rare[[as.character(i)]],
              function(x) x[,"L"])
  lines(x = 1:max_n_markers, y = y, col = cols[i], lwd = 2)
}
