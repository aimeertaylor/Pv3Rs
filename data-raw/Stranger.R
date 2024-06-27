################################################################################
################################################################################
rm(list = ls())
library(MCMCpack) # For rdirichlet
library(tictoc) # For timing

#===============================================================================
# Magic numbers
#===============================================================================
seed <- 2 # For reproducibility
c_params <- c(0.5, 1, 100) # Dirichlet concentration parameter
n_markers <- c(10, 50, 100) # Number of makers
n_repeats <- 5 # Number of simulations per c_param, n_marker combination
c_cutoff <- 99 # Above c_cutoff, switch from Dirichlet to 1/n_alleles
n_alleles <- 5 # Number of alleles per marker (marker cardinality)

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
min_marker_subset <- sample(all_markers, min_n_markers, replace = FALSE)

# Create progressive marker subsets, randomly sampled over chromosomes
marker_subsets <- list()
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

# Map the markers to chromosomes. Assume equal sized chromosomes; okay if
# later we assume an equal number of crossovers per chromosome
chrs_per_marker <- round(seq(0.51, 14.5, length.out = max_n_markers))
markers_per_chr <- table(chrs_per_marker)

#===============================================================================
# Generate data
#===============================================================================
tictoc::tic()
for(c in c_params) {

  # Sample allele frequencies
  fs <- sapply(all_markers, function(m) {
    if(c > c_cutoff) {
      fs_unnamed <- rep(1/n_alleles, n_alleles)
    } else {
      fs_unnamed <- MCMCpack::rdirichlet(1, alpha = rep(c, n_alleles))
    }
    setNames(fs_unnamed, alleles)
  }, USE.NAMES = TRUE, simplify = FALSE)
  fs_store[[as.character(c)]] <- fs

  for(i in 1:n_repeats) {

    print(paste(c,i))

    # Sample strangers
    any_clones <- TRUE
    while (any_clones) {
      parent1 <- sapply(all_markers, function(t) sample(alleles, size = 1, prob = fs[[t]]))
      parent2 <- sapply(all_markers, function(t) sample(alleles, size = 1, prob = fs[[t]]))
      stranger <- sapply(all_markers, function(t) sample(alleles, size = 1, prob = fs[[t]]))
      any_clones <- (identical(parent1[min_marker_subset], parent2[min_marker_subset]) |
                       identical(parent1[min_marker_subset], stranger[min_marker_subset]) |
                       identical(parent2[min_marker_subset], stranger[min_marker_subset]))
    }

    parents <- cbind(parent1, parent2)

    # Note condition on sufficient diversity
    children_clones <- TRUE
    while (children_clones) {

      # Sample parental allocations
      cs <- recombine_parent_ids(markers_per_chr)[,1:2]

      # Construct children genotypes from parental allocations
      children <- sapply(1:max_n_markers, function(i) {
        sapply(1:ncol(cs), function(j) {
          parents[i,cs[i,j]]
        })
      })
      colnames(children) <- all_markers

      # Check diversity among children in the first few markers
      children_clones <- nrow(unique(children[,min_marker_subset])) < 2
    }

    # Make parasite infection and data
    initial <- rbind(children[1,], children[2,])
    relapse <- rbind(stranger)

    # Format parasite infection data for compute_posterior
    y <- list(initial = apply(initial, 2, unique, simplify = F),
              relapse = apply(relapse, 2, unique, simplify = F))

    # Store the data
    ys_store[[as.character(c)]][[i]] <- y
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
      y_all_markers <- ys_store[[as.character(c)]][[i]]
      for(m in n_markers){
        marker_subset <- marker_subsets[[m]]
        y <- sapply(y_all_markers, function(x) x[marker_subset], simplify = FALSE)
        ps <- suppressMessages(compute_posterior(y, fs, return.RG = TRUE, return.logp = TRUE))
        ps_store[[as.character(c)]][[as.character(i)]][[as.character(m)]] <- ps
    }
  }
}
tictoc::toc()

#===============================================================================
# Generate results for equifrequent alleles and for all markers progressively
#===============================================================================
c <- 100 # For uniform allele frequencies only
fs <- fs_store[[as.character(c)]] # Extract frequencies
tictoc::tic()
for(i in 1:n_repeats){
  print(i)
    y_all_markers <- ys_store[[as.character(c)]][[i]]
    # compute posterior relapse probabilities for all marker counts from one onwards
    marker_subset <- min_marker_subset
    for(m in 1:max_n_markers){
      marker_subset <- marker_subsets[[m]]
      y <- sapply(y_all_markers, function(x) x[marker_subset], simplify = FALSE)
      ps <- suppressMessages(compute_posterior(y, fs))
      ps_store_all_ms[[as.character(i)]][[paste0("m",m)]] <- ps$marg[,"I"]
  }
}
tictoc::toc()


#===============================================================================
# Bundle data, results and magic numbers
#===============================================================================
Stranger <- list(fs_store = fs_store,
                      ys_store = ys_store,
                      ps_store = ps_store,
                      ps_store_all_ms = ps_store_all_ms,
                      c_cutoff = c_cutoff,
                      n_alleles = n_alleles,
                      min_n_markers = min_n_markers,
                      seed = seed)


#===============================================================================
# Save Half_siblings as exported data
#===============================================================================
usethis::use_data(Stranger, overwrite = TRUE)
