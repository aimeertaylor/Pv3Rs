################################################################################
# Bulk data from three of four meiotic siblings

# To do: check simulations for c("Stranger", "Clone", "Full_sibling", "Meiotic_sibling")
# are as expected; delete the data_raw/ scripts for these cases

################################################################################
rm(list = ls())
library(MCMCpack) # For rdirichlet
library(tictoc) # For timing

#===============================================================================
# Magic numbers / quantities
#===============================================================================
relapsing_parasites <- c("Stranger", "Clone", "Full_sibling", "Meiotic_sibling")
MOIs_per_infection_all <- c("2_1", "3_1") # Evaluate both for siblings only
n_alleles <- 5 # Number of alleles per marker (marker cardinality)
n_repeats <- 5 # Number of simulations per parameter combination
n_markers <- c(10, 50, 100) # Number of markers for which RG likelihood returned
c_params <- c(0.5, 1, 100) # Dirichlet concentration parameter
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
markers_per_chr <- table(chrs_per_marker)

#===============================================================================
# Generate data
#===============================================================================
for(relapsing_parasite in relapsing_parasites){

  print(relapsing_parasite)

  if(grepl("sibling", relapsing_parasite)) {
    MOIs_per_infection <- MOIs_per_infection_all
  } else {
    MOIs_per_infection <- MOIs_per_infection_all[1]
  }

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

    for(i in 1:n_repeats) {

      parent_clones <- TRUE
      while (parent_clones) {
        parent1 <- sapply(all_markers, function(t) sample(alleles, size = 1, prob = fs[[t]]))
        parent2 <- sapply(all_markers, function(t) sample(alleles, size = 1, prob = fs[[t]]))
        parent_clones <- identical(parent1[no_clone_subset], parent2[no_clone_subset])
      }
      parents <- cbind(parent1, parent2)

      children_clones <- TRUE
      while (children_clones) {

        # Sample parental allocations
        if(relapsing_parasite == "Full_sibling") {
          # independently for fourth child
          cs <- cbind(recombine_parent_ids(markers_per_chr)[,1:3],
                      recombine_parent_ids(markers_per_chr)[,1]) # Full sibling
        } else {
          # Sample parental allocations dependently
          cs <- recombine_parent_ids(markers_per_chr)
        }

        # Construct children genotypes from parental allocations
        children <- sapply(1:max_n_markers, function(i) {
          sapply(1:ncol(cs), function(j) parents[i,cs[i,j]])
        })
        colnames(children) <- all_markers

        # Check diversity among children in the minimum marker subset
        children_clones <- nrow(unique(children[,no_clone_subset])) < 4
      }

      # For different numbers of genotypes per infection
      for(MOIs in MOIs_per_infection) {

        # Make initial infection
        if (MOIs == "2_1") {
          initial <- children[1:2,]
        } else if (MOIs == "3_1") {
          initial <- children[1:3,] # Three of four meiotic siblings
        }

        if (relapsing_parasite == "Clone") {
          relapse <- rbind(children[1,])
        } else if (relapsing_parasite == "Stranger") {
          stranger <- sapply(all_markers, function(t) sample(alleles, size = 1, prob = fs[[t]]))
          relapse <- rbind(stranger)
        } else {
          relapse <- rbind(children[4,]) # Sibling
        }

        # Format parasite infection data for compute_posterior
        y <- list(initial = apply(initial, 2, unique, simplify = F),
                  relapse = apply(relapse, 2, unique, simplify = F))

        # Store the data
        ys_store[[as.character(c)]][[MOIs]][[as.character(i)]] <- y
      }
    }
  }
  tictoc::toc()



  #=============================================================================
  # Generate results with return.logp = TRUE for c_params and n_markers
  #=============================================================================
  tictoc::tic()
  for(i in 1:n_repeats){
    print(i)
    for(c in c_params) {
      fs <- fs_store[[as.character(c)]]
      for(MOIs in MOIs_per_infection) {
        y_all_markers <- ys_store[[as.character(c)]][[MOIs]][[as.character(i)]]
        for(m in n_markers){
          marker_subset <- marker_subsets[[m]]
          y <- sapply(y_all_markers, function(x) x[marker_subset], simplify = FALSE)
          ps <- suppressMessages(compute_posterior(y, fs, return.RG = TRUE, return.logp = TRUE))
          ps_store[[as.character(c)]][[MOIs]][[as.character(i)]][[as.character(m)]] <- ps
        }
      }
    }
  }
  tictoc::toc()

  #=============================================================================
  # Generate results for equifrequent alleles and markers 1:max_n_markers
  #=============================================================================
  c <- tail(c_params, 1) # For uniform allele frequencies only
  fs <- fs_store[[as.character(c)]] # Extract frequencies

  # Specifiy the recurrent state with largested expected posterior
  if (relapsing_parasite == "Clone") {
    exp_state <- "C"
  } else if (relapsing_parasite == "Stranger") {
    exp_state <- "I"
  } else {
    exp_state <- "L"
  }

  tictoc::tic()
  for(i in 1:n_repeats){
    print(i)
    for(MOIs in MOIs_per_infection) {
      y_all_markers <- ys_store[[as.character(c)]][[MOIs]][[as.character(i)]]
      for(m in 1:max_n_markers){
        marker_subset <- marker_subsets[[m]]
        y <- sapply(y_all_markers, function(x) x[marker_subset], simplify = FALSE)
        ps <- suppressMessages(compute_posterior(y, fs))
        ps_store_all_ms[[MOIs]][[as.character(i)]][[paste0("m",m)]] <- ps$marg[, exp_state]
      }
    }
  }
  tictoc::toc()

  # Extract probability of expected state
  post_S <- sapply(ps_store, function(X) {
    sapply(X, function(XX) {
      sapply(XX, function(XXX) {
        sapply(XXX, function(post) post$marg[,exp_state])
      })
    }, simplify = F)
  }, simplify = F)

  # Aside: check graphs are all ordered the same [make into a unit test?]
  # Extract graph summary (i.e., discard logp)
  justRGs <- sapply(ps_store, function(X) {
    sapply(X, function(XX) {
      sapply(XX, function(XXX) {
        sapply(XXX, function(post) {
          sapply(post$RGs, function(RG) c(RG$clone, RG$sib))
        }, simplify = F)
      }, simplify = F)
    }, simplify = F)
  }, simplify = F)

  # Check all the graphs are returned in the same order
  justRG <- justRGs[[1]][[1]][[1]][[1]]
  RGcheck <- sapply(c_params, function(c) {
    sapply(MOIs_per_infection, function(MOIs) {
      sapply(1:n_repeats, function(i) {
        sapply(n_markers, function(m) {
          identical(justRG, justRGs[[as.character(c)]][[MOIs]][[i]][[as.character(m)]])
        })
      })
    })
  })

  if (!all(RGcheck)) stop("graphs not returned in the same order")

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
                 n_alleles = n_alleles,
                 n_repeats = n_repeats,
                 n_markers = n_markers,
                 c_params = c_params,
                 c_cutoff = c_cutoff,
                 seed = seed,
                 fs_store = fs_store,
                 ys_store = ys_store,
                 ps_store = ps_store,
                 ps_store_all_ms = ps_store_all_ms,
                 post_S = post_S,
                 llikeRGs = llikeRGs)

  save(output, file = sprintf("%s.rda", relapsing_parasite))
}
