################################################################################
# Simulate data and generate results for a initial infection of two siblings and
# a recurrent sibling, where siblings are either parent child-like siblings or
# half siblings.
#
# In both cases, explore two scenarios. One where parents draw from the same
# allele distribution. Another with admixture where the first parent who parents
# intra-episode parasites draws alleles disproportionally:
#
# Half siblings:
# initial <- rbind(child12, child13) # migrant is 1
# relapse <- rbind(child23)
#
# Parent child-like siblings:
# initial <- rbind(child1, child12) # migrant is 1
# relapse <- rbind(child2)
#
# The alternative scenario parent-child like scenario with two selfed parasites
# in the initial infection is covered by the meiotic sibling case with MOIs 3 &
# 1 and no extenal MOI specification (see vignette on understanding posterior
# estimates).
#
# For each case, generate results for data on all marker counts when alleles are
# equifrequent, and for a subset of marker counts when alleles are not.
################################################################################
rm(list = ls())
library(Pv3Rs)
library(MCMCpack) # For rdirichlet

#===============================================================================
# Magic numbers / quantities
#===============================================================================
cases <- c("PCLike", "Half")
n_alleles <- 5 # Number of alleles per marker (marker cardinality)
n_repeats <- 10 # Number of simulations per parameter combination
n_markers <- c(10, 50, 100, 150) # Number of markers for which RG likelihood returned
c_params <- c(0.5, 1, 10, 100) # Dirichlet concentration parameter
c_cutoff <- 99 # Switch from Dirichlet r.v. to 1/n_alleles above c_cutoff
seed <- 1 # For reproducibility

#===============================================================================
# Stores for data, frequencies & results
#===============================================================================
output_HalfSib.PClikeSib <- list()
ys_store <- list() # y for data
fs_store <- list() # f for frequency
ps_store <- list() # p for posterior
ps_store_all_ms <- list() # all_ms for all marker counts
ps_store_all_ms_admix_rare <- list() # for admixed case with some rare freqs

#===============================================================================
# Set the seed, name markers and get alleles, maker subsets etc.
#===============================================================================
set.seed(seed) # Set the seed
min_n_markers <- min(n_markers)
max_n_markers <- max(n_markers)
alleles <- letters[1:n_alleles]
all_markers <- paste0("m", 1:max_n_markers) # Marker names
rorder <- sample(all_markers, size = max_n_markers) # Many markers, random order
marker_subsets <- lapply(1:max_n_markers, function(m) rorder[1:m]) # Subsets
# smallest subset over which clones are disallowed (ensures the set of RGs is
# the same for all n_markers):
no_clone_subset <- marker_subsets[[min_n_markers]]

# Map the markers to chromosomes. Assume equally sized chromosomes — reasonable
# providing we later assume an equal number of crossovers per chromosome
chrs_per_marker <- round(seq(0.51, 14.5, length.out = max_n_markers))

for(case in cases) {

  print(case)

  #===============================================================================
  # Generate data
  # When admixture is TRUE, draw rare alleles with high probability for the
  # parent that parents intra-episode parasite only.
  #===============================================================================
  for(c in c_params) {

    # Sample allele frequencies
    fs <- sapply(all_markers, function(m) {
      if(c > c_cutoff) {
        fs_unnamed <- rep(1/n_alleles, n_alleles)
      } else {
        fs_unnamed <- MCMCpack::rdirichlet(1, rep(c, n_alleles))
      }
      setNames(fs_unnamed, alleles)
    }, USE.NAMES = TRUE, simplify = FALSE)

    # Store allele frequencies
    fs_store[[as.character(c)]] <- fs

    for(admixture in c(TRUE, FALSE)) {
      for(i in 1:n_repeats) {

        if (case == "PCLike") {
          # Sample parental genotypes
          parent_clones <- TRUE
          while (parent_clones) {

            if (admixture) {  # Draw rare alleles with high probability
              parent1 <- sapply(all_markers, function(t) sample(alleles, size = 1, prob = 1-fs[[t]]))
            } else {
              parent1 <- sapply(all_markers, function(t) sample(alleles, size = 1, prob = fs[[t]]))
            }
            parent2 <- sapply(all_markers, function(t) sample(alleles, size = 1, prob = fs[[t]]))
            parent_clones <- identical(parent1[no_clone_subset], parent2[no_clone_subset])
          }

          parents <- cbind(parent1, parent2)

          # Sample children genotypes
          children_clones <- TRUE
          while (children_clones) {

            # Sample parental allocations
            cs <- recombine_parent_ids(chrs_per_marker)[,1]
            names(cs) <- all_markers

            # Create recombinant
            child12 <- sapply(all_markers, function(m) parents[m,cs[m]])
            child1 <- parent1
            child2 <- parent2

            children_clones <- identical(child1[no_clone_subset], child12[no_clone_subset])
          }

          # Make parasite infection and data
          initial <- rbind(child1, child12)
          relapse <- rbind(child2)

        } else if (case == "Half") {

          # Sample parental genotypes
          parent_clones <- TRUE
          while (parent_clones) {
            if (admixture) { # Draw rare alleles with high probability
              parent1 <- sapply(all_markers, function(t) sample(alleles, size = 1, prob = 1-fs[[t]]))
            } else {
              parent1 <- sapply(all_markers, function(t) sample(alleles, size = 1, prob = fs[[t]]))
            }
            parent2 <- sapply(all_markers, function(t) sample(alleles, size = 1, prob = fs[[t]]))
            parent3 <- sapply(all_markers, function(t) sample(alleles, size = 1, prob = fs[[t]]))
            parent_clones <- any(identical(parent1[no_clone_subset], parent2[no_clone_subset]),
                                 identical(parent1[no_clone_subset], parent3[no_clone_subset]),
                                 identical(parent2[no_clone_subset], parent3[no_clone_subset]))
          }

          parents12 <- cbind(parent1, parent2)
          parents13 <- cbind(parent1, parent3)
          parents23 <- cbind(parent2, parent3)

          # Sample children genotypes
          children_clones <- TRUE
          while (children_clones) {

            # Sample parental allocations independently
            cs <- sapply(1:3, function(i) recombine_parent_ids(chrs_per_marker)[,1])
            rownames(cs) <- all_markers

            # Construct children genotypes from parental allocations
            child12 <- sapply(all_markers, function(i) parents12[[i,cs[i,1]]])
            child13 <- sapply(all_markers, function(i) parents13[[i,cs[i,2]]])
            child23 <- sapply(all_markers, function(i) parents23[[i,cs[i,3]]])

            # Check for clones in MOI = 2 infection only
            children_clones <- identical(child12[no_clone_subset], child13[no_clone_subset])
          }

          # Make parasite infection and data
          initial <- rbind(child12, child13)
          relapse <- rbind(child23)

        } else {
          stop('case should either be "PCLike" or "Half"')
        }

        y <- list(initial = apply(initial, 2, unique, simplify = F),
                  relapse = apply(relapse, 2, unique, simplify = F))

        # Store the data
        ys_store[[as.character(c)]][[sprintf("admixture_%s", admixture)]][[i]] <- y

      }
    }
  }

  #=============================================================================
  # Generate results for markers 1:max_n_markers
  #=============================================================================
  c <- 100 # For uniform allele frequencies only
  admixture <- FALSE # For parents from the same population
  fs <- fs_store[[as.character(c)]] # Extract frequencies
  for(i in 1:n_repeats){
    y_all_markers <- ys_store[[as.character(c)]][[sprintf("admixture_%s", admixture)]][[i]]
    for(m in 1:max_n_markers){
      marker_subset <- marker_subsets[[m]]
      y <- sapply(y_all_markers, function(x) x[marker_subset], simplify = FALSE)
      ps <- suppressMessages(compute_posterior(y, fs))
      ps_store_all_ms[[as.character(i)]][[paste0("m",m)]] <- ps$marg
    }
  }

  c <- 0.5 # For admixed case with imbalanced allele frequencies
  admixture <- TRUE # For parents from different populations
  fs <- fs_store[[as.character(c)]] # Extract frequencies
  for(i in 1:n_repeats){
    y_all_markers <- ys_store[[as.character(c)]][[sprintf("admixture_%s", admixture)]][[i]]
    for(m in 1:max_n_markers){
      marker_subset <- marker_subsets[[m]]
      y <- sapply(y_all_markers, function(x) x[marker_subset], simplify = FALSE)
      ps <- suppressMessages(compute_posterior(y, fs))
      ps_store_all_ms_admix_rare[[as.character(i)]][[paste0("m",m)]] <- ps$marg
    }
  }

  #=============================================================================
  # Generate results with return.logp = TRUE
  #=============================================================================
  for(i in 1:n_repeats){
    for(c in c_params) {
      fs <- fs_store[[as.character(c)]]
      for(admixture in c(TRUE, FALSE)) {
        y_all_markers <- ys_store[[as.character(c)]][[sprintf("admixture_%s", admixture)]][[i]]
        for(m in n_markers){
          marker_subset <- marker_subsets[[m]]
          y <- sapply(y_all_markers, function(x) x[marker_subset], simplify = FALSE)
          ps <- suppressMessages(compute_posterior(y, fs, return.RG = TRUE, return.logp = TRUE))
          ps_store[[as.character(c)]][[sprintf("admixture_%s", admixture)]][[as.character(i)]][[as.character(m)]] <- ps
        }
      }
    }
  }

  #=============================================================================
  # Bundle magic numbers, data and results
  #=============================================================================
  output <- list(n_alleles = n_alleles,
                 n_repeats = n_repeats,
                 n_markers = n_markers,
                 c_params = c_params,
                 c_cutoff = c_cutoff,
                 seed = seed,
                 fs_store = fs_store,
                 ys_store = ys_store,
                 ps_store = ps_store,
                 ps_store_all_ms = ps_store_all_ms,
                 ps_store_all_ms_admix_rare = ps_store_all_ms_admix_rare)

  output_HalfSib.PClikeSib[[case]] <- output
}

# Save as exported data
usethis::use_data(output_HalfSib.PClikeSib, overwrite = TRUE)
