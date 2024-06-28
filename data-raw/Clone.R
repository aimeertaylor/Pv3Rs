################################################################################
# To do: check simulations for c("Stranger", "Clone", "Full_sibling", "Meiotic_sibling")
# are as expected; delete the data_raw/ scripts for these cases

################################################################################
rm(list = ls())
library(MCMCpack) # For rdirichlet
library(tictoc) # For timing

#===============================================================================
# Magic numbers / quantities
#===============================================================================
seed <- 2 # For reproducibility
n_alleles <- 5 # Number of alleles per marker (marker cardinality)
n_repeats <- 2 # Number of simulations per parameter combination
n_markers <- c(10, 50, 100) # Number of markers
c_params <- c(0.5, 1, 100) # Dirichlet concentration parameter
c_cutoff <- 99 # Above c_cutoff, switch from Dirichlet to 1/n_alleles
relapsing_parasites <- c("Stranger", "Clone", "Full_sibling", "Meiotic_sibling")
MOIs_per_infection_all <- c("2_1", "3_1")

#===============================================================================
# Stores for data, frequencies & results
#===============================================================================
ys_store <- list() # y for data
fs_store <- list() # f for frequency
ps_store <- list() # p for posterior
ps_store_all_ms <- list() # all_ms for all marker counts

#===============================================================================
# Set the seed, name markers and get alleles
#===============================================================================
set.seed(seed) # Set the seed
max_n_markers <- max(n_markers)
min_n_markers <- min(n_markers) # Size of marker subset for which clones are disallowed
alleles <- letters[1:n_alleles] # Alleles
all_markers <- paste0("m", 1:max_n_markers) # Marker names

# Create progressive marker subsets, randomly sampled over chromosomes
marker_subsets <- list()
min_marker_subset <- sample(all_markers, min_n_markers, replace = FALSE)
marker_subsets[[1]] <- sample(min_marker_subset, 1)
for(m in 2:min_n_markers){
  markers_thusfar <- marker_subsets[[m-1]]
  marker_to_add <- sample(setdiff(min_marker_subset, markers_thusfar), 1)
  marker_subsets[[m]] <- c(markers_thusfar, marker_to_add)
}
for(m in min_n_markers:(max_n_markers-1)) {
  markers_thusfar <- marker_subsets[[m]]
  marker_to_add <- sample(setdiff(all_markers, markers_thusfar), 1)
  marker_subsets[[m+1]] <- c(markers_thusfar, marker_to_add)
}

# Map the markers to chromosomes. Assume equal sized chromosomes — reasonable
# providing we later assume an equal number of crossovers per chromosome
chrs_per_marker <- round(seq(0.51, 14.5, length.out = max_n_markers))
markers_per_chr <- table(chrs_per_marker)

# Function to sample allele frequencies
fs_sample <- function(ms, c, c_cutoff){
  x <- sapply(ms, function(m) {
    if(c > c_cutoff) {
      fs_unnamed <- rep(1/n_alleles, n_alleles)
    } else {
      fs_unnamed <- MCMCpack::rdirichlet(1, alpha = rep(c, n_alleles))
    }
    setNames(fs_unnamed, alleles)
  }, USE.NAMES = TRUE, simplify = FALSE)
  return(x)
}



#===============================================================================
# Generate data
#===============================================================================
tictoc::tic()
for(relapsing_parasite in relapsing_parasites){

  if(grepl("sibling", relapsing_parasite)) {
    MOIs_per_infection <- MOIs_per_infection_all
  } else {
    MOIs_per_infection <- MOIs_per_infection_all[1]
  }

  for(c in c_params) {

    fs <- fs_sample(all_markers, c, c_cutoff)
    fs_store[[as.character(c)]] <- fs

    for(i in 1:n_repeats) {

      # (ensure no clones for min_marker_subset and thus always the same no. of
      # vertices in RGs for n_markers)
      parent_clones <- TRUE
      while (parent_clones) {
        parent1 <- sapply(all_markers, function(t) sample(alleles, size = 1, prob = fs[[t]]))
        parent2 <- sapply(all_markers, function(t) sample(alleles, size = 1, prob = fs[[t]]))
        parent_clones <- identical(parent1[min_marker_subset], parent2[min_marker_subset])
      }
      parents <- cbind(parent1, parent2)

      # (ensure no clones for min_marker_subset and thus always the same no. of
      # vertices in RGs for n_markers)
      children_clones <- TRUE
      while (children_clones) {

        if(relapsing_parasite == "Full_sibling") {
          # Sample parental allocations for child 4 independently from 1-3
          cs1 <- recombine_parent_ids(markers_per_chr)[,1:3]
          cs2 <- recombine_parent_ids(markers_per_chr)[,1]
          cs <- cbind(cs1, cs2)
        } else {
          # Sample parental allocations dependently
          cs <- recombine_parent_ids(markers_per_chr)
        }

        # Construct children genotypes from parental allocations
        children <- sapply(1:max_n_markers, function(i) {
          sapply(1:ncol(cs), function(j) {
            parents[i,cs[i,j]]
          })
        })
        colnames(children) <- all_markers

        # Check diversity among children in the minimum marker subset
        children_clones <- nrow(unique(children[,min_marker_subset])) < 4
      }

      # For different numbers of genotypes per infection
      for(MOIs in MOIs_per_infection) {

        # Make initial infection
        if (MOIs == "2_1") {
          initial <- children[1:2,]
        } else if (MOIs == "3_1") {
          initial <- children[1:3,]
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
        ys_store[[as.character(c)]][[MOIs]][[i]] <- y
      }
    }
  }
  tictoc::toc()

  #===============================================================================
  # Generate results for different allele frequency types and for
  # different marker counts
  #===============================================================================
  tictoc::tic()
  for(i in 1:n_repeats){
    print(i)
    for(c in c_params) {
      fs <- fs_store[[as.character(c)]]
      for(MOIs in MOIs_per_infection) {
        y_all_markers <- ys_store[[as.character(c)]][[MOIs]][[i]]
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

  #===============================================================================
  # Generate results for equifrequent alleles and for all markers progressively
  #===============================================================================
  c <- 100 # For uniform allele frequencies only
  fs <- fs_store[[as.character(c)]] # Extract frequencies

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
      y_all_markers <- ys_store[[as.character(c)]][[MOIs]][[i]]
      # compute posterior relapse probabilities for all marker counts from one onwards
      marker_subset <- min_marker_subset
      for(m in 1:max_n_markers){
        marker_subset <- marker_subsets[[m]]
        y <- sapply(y_all_markers, function(x) x[marker_subset], simplify = FALSE)
        ps <- suppressMessages(compute_posterior(y, fs))
        ps_store_all_ms[[MOIs]][[as.character(i)]][[paste0("m",m)]] <- ps$marg[, exp_state]
      }
    }
  }
  tictoc::toc()


  #===============================================================================
  # Bundle data, results and magic numbers
  #===============================================================================
  output <- list(fs_store = fs_store,
                 ys_store = ys_store,
                 ps_store = ps_store,
                 ps_store_all_ms = ps_store_all_ms,
                 c_cutoff = c_cutoff,
                 n_alleles = n_alleles,
                 min_n_markers = min_n_markers,
                 seed = seed)

  save(output, file = sprintf("%s.rda", relapsing_parasite))
}
