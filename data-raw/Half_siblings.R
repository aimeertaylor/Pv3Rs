################################################################################
# Code to prepare Half_siblings dataset following best practice outlined in
# 7.1.1 of https://r-pkgs.org/data.html. This script only features in the source
# version of the package (it is listed under .Rbuildignore). For c_params in
# 0.1, 1, 10, 100, n_markers in 10, 50, 100 and 150, and 10 repeats, takes a
# couple of minutes to run on a powerful MacBook Pro with 32 GB memory.
################################################################################
rm(list = ls())
library(MCMCpack) # For rdirichlet
library(tictoc) # For timing

#===============================================================================
# Magic numbers
#===============================================================================
seed <- 1 # For reproducibility
c_params <- c(0.1, 1, 10, 100) # Dirichlet concentration parameter
n_markers <- c(10, 50, 100, 150) # Number of makers
n_repeats <- 10 # Number of simulations per c_param, n_marker combination
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
all_markers <- paste0("m", 1:max_n_markers) # Marker names
alleles <- letters[1:n_alleles] # Alleles
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
#
# If rare_enrich = TRUE, draw rare alleles with high probability for one of the
# parents who fathers the intra-episode siblings (interpret as an migrant
# from another population); otherwise, draw alleles proportional to allele
# frequencies for all parents (parents from a single population)
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

  for(rare_enrich in c(TRUE, FALSE)) {
    for(i in 1:n_repeats) {

      print(paste(c,rare_enrich,i))

      # Sample parental genotypes
      # (ensure no clones and thus always the same no. of vertices in RGs)
      anyclones <- TRUE
      while (anyclones) {
        if (rare_enrich) {
          parent1 <- sapply(all_markers, function(t) {
            sample(alleles, size = 1, prob = 1-fs[[t]])}, simplify = F)
        } else {
          parent1 <- sapply(all_markers, function(t) {
            sample(alleles, size = 1, prob = fs[[t]])}, simplify = F)
        }
        parent2 <- sapply(all_markers, function(t) {
          sample(alleles, size = 1, prob = fs[[t]])}, simplify = F)
        parent3 <- sapply(all_markers, function(t) {
          sample(alleles, size = 1, prob = fs[[t]])}, simplify = F)

        anyclones <- any(identical(parent1[min_marker_subset], parent2[min_marker_subset]),
                         identical(parent1[min_marker_subset], parent3[min_marker_subset]),
                         identical(parent2[min_marker_subset], parent3[min_marker_subset]))
      }

      parents12 <- cbind(parent1, parent2)
      parents13 <- cbind(parent1, parent3)
      parents23 <- cbind(parent2, parent3)

      # Sample children genotypes independently (ensure no intra-clones)
      clones <- TRUE
      while (clones) {

        # Sample parental allocations independently
        ps <- sapply(1:3, function(i) recombine_parent_ids(markers_per_chr)[,1])
        rownames(ps) <- all_markers

        # Construct children genotypes from parental allocations
        child12 <- sapply(all_markers, function(i) parents12[[i,ps[i,1]]])
        child13 <- sapply(all_markers, function(i) parents13[[i,ps[i,2]]])
        child23 <- sapply(all_markers, function(i) parents23[[i,ps[i,3]]])

        # Check for clones in MOI = 2 infection only
        clones <- identical(child12[min_marker_subset], child13[min_marker_subset])
      }

      # Make parasite infection and data
      initial <- rbind(unlist(child12), unlist(child13))
      relapse <- rbind(unlist(child23))
      y <- list(initial = apply(initial, 2, unique, simplify = F),
                relapse = apply(relapse, 2, unique, simplify = F))

      # Store the data
      ys_store[[as.character(c)]][[sprintf("rare_enrich_%s", rare_enrich)]][[i]] <- y

    }
  }
}
tictoc::toc()

#===============================================================================
# Generate results for different allele frequency types, for a migrant parent
# (rare_enrich TRUE) as well as parents from a single population, and for
# different marker counts
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
# Generate results for equifrequent alleles, for parents from a single
# population, and for all marker counts from one onwards.
#===============================================================================
c <- 100 # For uniform allele frequencies only
rare_enrich <- FALSE # For parents from the same population
fs <- fs_store[[as.character(c)]] # Extract frequencies
tictoc::tic()
for(i in 1:n_repeats){
  print(i)
  y_all_markers <- ys_store[[as.character(c)]][[sprintf("rare_enrich_%s", rare_enrich)]][[i]]
  # compute posterior relapse probabilities for all marker counts from one onwards
  for(m in 1:max_n_markers){
    marker_subset <- marker_subsets[[m]]
    y <- sapply(y_all_markers, function(x) x[marker_subset], simplify = FALSE)
    ps <- suppressMessages(compute_posterior(y, fs))
    ps_store_all_ms[[as.character(i)]][[paste0("m",m)]] <- ps$marg[,"L"]
  }
}
tictoc::toc()


#===============================================================================
# Bundle data, results and magic numbers
#===============================================================================
Half_siblings <- list(fs_store = fs_store,
                      ys_store = ys_store,
                      ps_store = ps_store,
                      ps_store_all_ms = ps_store_all_ms,
                      c_cutoff = c_cutoff,
                      n_alleles = n_alleles,
                      seed = seed)


#===============================================================================
# Save Half_siblings as exported data
#===============================================================================
usethis::use_data(Half_siblings, overwrite = TRUE)
