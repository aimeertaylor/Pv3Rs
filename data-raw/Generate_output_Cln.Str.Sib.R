################################################################################
# Simulate data and generate results for a initial episode with MOI 2 or 3
# meiotic siblings and a monoclonal recurrence, where the recurrent parasite is
# either a stranger, a clone, a regular sibling, or a meiotic sibling. For each
# case, generate results for all marker counts when alleles are equifrequent,
# and for a subset of marker counts otherwise.
#
# Doesn't take two long to run: minute or two on 32 GB RAM laptop
################################################################################
rm(list = ls())
library(Pv3Rs)
library(MCMCpack) # For rdirichlet

#===============================================================================
# Magic numbers / quantities
#===============================================================================
provide_correct_MOIs <- FALSE # Toggle for inference with correct external MOIs
cases <- c("Stranger", "Clone", "Regular_sibling", "Meiotic_sibling")
MOIs_per_infection <- c("2_1", "3_1") # Change num. of meiotic sibs in 1st epi.
n_alleles <- 5 # Number of alleles per marker (marker cardinality)
n_repeats <- 5 # Number of simulations per parameter combination
n_markers <- c(10, 50, 100) # Number of markers for which RG likelihood returned
c_params <- c(0.5, 1, 100) # Dirichlet concentration parameter
c_cutoff <- 99 # Switch from Dirichlet r.v. to 1/n_alleles above c_cutoff
seed <- 1 # For reproducibility

#===============================================================================
# Stores for data, frequencies & results
#===============================================================================
output_Cln.Str.Sib <- list() # For output
ys_store <- list() # y for data
fs_store <- list() # f for frequency
ps_store <- list() # p for posterior
ps_store_all_ms <- list() # for all marker counts (ms)

#===============================================================================
# Set the seed, name markers and get alleles, maker subsets etc.
#===============================================================================
set.seed(seed) # Set the seed
min_n_markers <- min(n_markers)
max_n_markers <- max(n_markers)
alleles <- letters[1:n_alleles]
all_markers <- paste0("m", 1:max_n_markers) # Marker names
rorder <- sample(all_markers, size = max_n_markers) # Many markers, random order
marker_subsets <- lapply(1:max_n_markers, function(m) rorder[1:m]) # Subset
# Smallest subset for which clones are disallowed (ensures the set of RGs is
# the same for all n_markers):
no_clone_subset <- marker_subsets[[min_n_markers]]

# Map the markers to chromosomes. Assume equally sized chromosomes — reasonable
# providing we later assume an equal number of crossovers per chromosome
chrs_per_marker <- round(seq(0.51, 14.5, length.out = max_n_markers))

for(case in cases){

  print(case)

  #===============================================================================
  # Generate data
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

    for(i in 1:n_repeats) {

      parent_clones <- TRUE
      while (parent_clones) { # Avoid clonal parents
        parent1 <- sapply(all_markers, function(t) sample(alleles, size = 1, prob = fs[[t]]))
        parent2 <- sapply(all_markers, function(t) sample(alleles, size = 1, prob = fs[[t]]))
        parent_clones <- identical(parent1[no_clone_subset], parent2[no_clone_subset])
      }
      parents <- cbind(parent1, parent2)

      children_clones <- TRUE # Avoid clonal children
      while (children_clones) {

        # Sample parental allocations
        if(case == "Regular_sibling") {
          # independently for fourth child
          cs <- cbind(recombine_parent_ids(chrs_per_marker)[,1:3],
                      recombine_parent_ids(chrs_per_marker)[,1]) # Regular sibling
        } else {
          # Sample parental allocations dependently
          cs <- recombine_parent_ids(chrs_per_marker)
        }

        # Construct children from parental allocations
        children <- sapply(1:max_n_markers, function(i) {
          sapply(1:ncol(cs), function(j) parents[i,cs[i,j]])
        })
        colnames(children) <- all_markers

        # Check diversity among children in the minimum marker subset
        children_clones <- nrow(unique(children[,no_clone_subset])) < 4
      }

      # For different numbers of genotypes per infection
      for(MOIs in MOIs_per_infection) {

        # Make initial infection: either two or three of four meiotic siblings
        if (MOIs == "2_1") {
          initial <- children[1:2,] # Two of four meiotic siblings
        } else if (MOIs == "3_1") {
          initial <- children[1:3,] # Three of four meiotic siblings
        }

        # Get recurrent parasite
        if (case == "Clone") {
          relapse <- rbind(children[1,])
        } else if (case == "Stranger") {
          stranger <- sapply(all_markers, function(t) sample(alleles, size = 1, prob = fs[[t]]))
          relapse <- rbind(stranger)
        } else {
          relapse <- rbind(children[4,]) # Sibling
        }

        # Format parasite episode data for compute_posterior
        y <- list(initial = apply(initial, 2, unique, simplify = F),
                  relapse = apply(relapse, 2, unique, simplify = F))

        # Store the data
        ys_store[[as.character(c)]][[MOIs]][[as.character(i)]] <- y
      }
    }
  }


  #=============================================================================
  # Generate results for markers 1:max_n_markers
  #=============================================================================
  c <- 100 # For uniform allele frequencies only
  fs <- fs_store[[as.character(c)]] # Extract frequencies

  # Specify the recurrent state with largest expected posterior
  if (case == "Clone") {
    exp_state <- "C"
  } else if (case == "Stranger") {
    exp_state <- "I"
  } else {
    exp_state <- "L"
  }

  for(i in 1:n_repeats){
    for(MOIs in MOIs_per_infection) {
      y_all_markers <- ys_store[[as.character(c)]][[MOIs]][[as.character(i)]]
      for(m in 1:max_n_markers){
        marker_subset <- marker_subsets[[m]]
        y <- sapply(y_all_markers, function(x) x[marker_subset], simplify = F)
        ps <- suppressMessages(compute_posterior(y, fs))
        ps_store_all_ms[[MOIs]][[as.character(i)]][[paste0("m",m)]] <- ps$marg
      }
    }
  }

  #=============================================================================
  # Generate results with return.logp = TRUE
  #=============================================================================
  for(i in 1:n_repeats){
    for(c in c_params) {
      fs <- fs_store[[as.character(c)]]
      for(MOIs in MOIs_per_infection) {
        y_all_markers <- ys_store[[as.character(c)]][[MOIs]][[as.character(i)]]
        for(m in n_markers){
          marker_subset <- marker_subsets[[m]]
          y <- sapply(y_all_markers, function(x) x[marker_subset], simplify = F)
          if(provide_correct_MOIs) {
            ps <- suppressMessages(compute_posterior(y, fs, MOIs = as.numeric(strsplit(MOIs, "_")[[1]]), return.RG = TRUE, return.logp = TRUE))
          } else {
            ps <- suppressMessages(compute_posterior(y, fs, return.RG = TRUE, return.logp = TRUE))
          }
          ps <- suppressMessages(compute_posterior(y, fs, return.RG = TRUE, return.logp = TRUE))
          ps_store[[as.character(c)]][[MOIs]][[as.character(i)]][[as.character(m)]] <- ps
        }
      }
    }
  }

  # Extract probability of expected state
  post_S <- sapply(ps_store, function(X) {
    sapply(X, function(XX) {
      sapply(XX, function(XXX) {
        sapply(XXX, function(post) post$marg[,exp_state])
      })
    }, simplify = F)
  }, simplify = F)

  # Extract probability of the data given the relationship graph
  llikeRGs <- sapply(ps_store, function(X) {
    sapply(X, function(XX) {
      sapply(XX, function(XXX) {
        sapply(XXX, function(post) sapply(post$RGs, function(RG) RG$logp))
      }, simplify = F)
    }, simplify = F)
  }, simplify = F)

  #=============================================================================
  # Bundle magic numbers, data, and results
  #=============================================================================
  output <- list(MOIs_per_infection = MOIs_per_infection,
                 provide_correct_MOIs = provide_correct_MOIs,
                 n_alleles = n_alleles,
                 n_repeats = n_repeats,
                 n_markers = n_markers,
                 c_params = c_params,
                 c_cutoff = c_cutoff,
                 seed = seed,
                 exp_state = exp_state,
                 fs_store = fs_store,
                 ys_store = ys_store,
                 ps_store = ps_store,
                 ps_store_all_ms = ps_store_all_ms,
                 post_S = post_S,
                 llikeRGs = llikeRGs)

  output_Cln.Str.Sib[[case]] <- output
}

# Save as exported data
usethis::use_data(output_Cln.Str.Sib, overwrite = TRUE)
