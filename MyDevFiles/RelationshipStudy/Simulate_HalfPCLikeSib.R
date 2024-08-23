################################################################################
################################################################################
rm(list = ls())
library(Pv3Rs)
library(MCMCpack) # For rdirichlet
library(tictoc) # For timing

#===============================================================================
# Magic numbers / quantities
#===============================================================================
cases <- c("ParentChildLike", "Half")
n_alleles <- 5 # Number of alleles per marker (marker cardinality)
n_repeats <- 10 # Number of simulations per parameter combination
n_markers <- c(10, 50, 100, 150) # Number of markers for which RG likelihood returned
c_params <- c(0.5, 1, 10, 100) # Dirichlet concentration parameter
c_cutoff <- 99 # Switch from Dirichlet r.v. to 1/n_alleles above c_cutoff
seed <- 1 # For reproducibility

#===============================================================================
# Stores for data, frequencies & results
#===============================================================================
ys_store <- list() # y for data
fs_store <- list() # f for frequency
ps_store <- list() # p for posterior
ps_store_all_ms <- list() # all_ms for all marker counts

#===============================================================================
# Set the seed, name markers and get alleles, maker subsets etc.
#===============================================================================
set.seed(seed) # Set the seed
min_n_markers <- min(n_markers)
max_n_markers <- max(n_markers)
alleles <- letters[1:n_alleles]
all_markers <- paste0("m", 1:max_n_markers) # Marker names

# Create marker subsets, randomly sampled over chromosomes
marker_subsets <- list(sample(all_markers, 1))
rorder <- sample(all_markers, size = max_n_markers)
for(m in 2:max_n_markers) marker_subsets[[m]] <- rorder[1:m]
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
  # When rare_enrich is TRUE, draw rare alleles with high probability for one
  # parent (interpret as an migrant from another population) who parents
  # intra-episode parasite.
  #===============================================================================
  tictoc::tic()
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

    for(rare_enrich in c(TRUE, FALSE)) {
      for(i in 1:n_repeats) {

        if (case == "ParentChildLike") {
          # Sample parental genotypes
          parent_clones <- TRUE
          while (parent_clones) {

            # draw rare alleles with high probability
            if (rare_enrich) {
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
            if (rare_enrich) { # draw rare alleles with high probability
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
          stop('case should either be "ParentChildLike" or "Half"')
        }

        y <- list(initial = apply(initial, 2, unique, simplify = F),
                  relapse = apply(relapse, 2, unique, simplify = F))

        # Store the data
        ys_store[[as.character(c)]][[sprintf("rare_enrich_%s", rare_enrich)]][[i]] <- y

      }
    }
  }
  tictoc::toc()

  #===============================================================================
  # Generate results with return.logp = TRUE
  #===============================================================================
  tictoc::tic()
  for(i in 1:n_repeats){
    print(i)
    for(c in c_params) {
      fs <- fs_store[[as.character(c)]]
      for(rare_enrich in c(TRUE, FALSE)) {
        y_all_markers <- ys_store[[as.character(c)]][[sprintf("rare_enrich_%s", rare_enrich)]][[i]]
        for(m in n_markers){
          marker_subset <- marker_subsets[[m]]
          y <- sapply(y_all_markers, function(x) x[marker_subset], simplify = FALSE)
          ps <- suppressMessages(compute_posterior(y, fs, return.RG = TRUE, return.logp = TRUE))
          ps_store[[as.character(c)]][[sprintf("rare_enrich_%s", rare_enrich)]][[as.character(i)]][[as.character(m)]] <- ps
        }
      }
    }
  }
  tictoc::toc()

  #===============================================================================
  # Generate results for markers 1:max_n_markers
  #===============================================================================
  c <- 100 # For uniform allele frequencies only
  rare_enrich <- FALSE # For parents from the same population
  fs <- fs_store[[as.character(c)]] # Extract frequencies
  tictoc::tic()
  for(i in 1:n_repeats){
    print(i)
    y_all_markers <- ys_store[[as.character(c)]][[sprintf("rare_enrich_%s", rare_enrich)]][[i]]
    for(m in 1:max_n_markers){
      marker_subset <- marker_subsets[[m]]
      y <- sapply(y_all_markers, function(x) x[marker_subset], simplify = FALSE)
      ps <- suppressMessages(compute_posterior(y, fs))
      ps_store_all_ms[[as.character(i)]][[paste0("m",m)]] <- ps$marg
    }
  }
  tictoc::toc()

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
                 ps_store_all_ms = ps_store_all_ms)


save(output, file = sprintf("%s_siblings.rda", case))
}
