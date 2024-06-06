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
max_n_markers <- max(n_markers)
n_repeats <- 10 # Number of simulations per c_param, n_marker combination
c_cutoff <- 99 # Above c_cutoff, switch from Dirichlet to 1/n_alleles
n_alleles <- 3 # Number of alleles per marker (marker cardinality)
m_min_clone <- 4 # Minimum number of markers for which clone is disallowed

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
all_markers <- paste0("m", 1:max_n_markers) # Marker names
alleles <- letters[1:n_alleles] # Alleles

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

      print(paste(c,rare_enrich,i))

      # Sample parental genotypes
      # (ensure no clones and thus always the same no. of vertices in RGs)
      anyclones <- TRUE
      while (anyclones) {
        parent1 <- sapply(all_markers, function(t) {
          sample(alleles, size = 1, prob = fs[[t]])}, simplify = F)
        parent2 <- sapply(all_markers, function(t) {
          sample(alleles, size = 1, prob = fs[[t]])}, simplify = F)
        anyclones <- any(identical(parent1[1:m_min_clone], parent2[1:m_min_clone]))
      }


      # Map the markers to chromosomes. Assume equal sized chromosomes; okay if
      # later we assume an equal number of crossovers per chromosome
      chrs <- round(seq(0.51, 14.5, length.out = max_n_markers))
      marker_per_chr <- table(chrs)

      # Per chromosome parent1 segment length for recombinant chromatids one and two
      c1_p1_segment_length <- sapply(marker_per_chr, sample, size = 1)
      c2_p1_segment_length <- sapply(marker_per_chr, sample, size = 1)

      # Crossover for chromatid one
      c1 <- do.call(c, sapply(1:14, function(i) {
          x <- rep(2, marker_per_chr[i])
          x[1:c1_p1_segment_length[i]] <- 1
          return(x)
      }))

      # Crossover for chromatid two
      c2 <- do.call(c, sapply(1:14, function(i) {
        x <- rep(2, marker_per_chr[i])
        x[1:c2_p1_segment_length[i]] <- 1
        return(x)
      }))

      # Complements of c1 and c2
      c3 <- abs(c1-2) + 1
      c4 <- abs(c2-2) + 1

      # Check complements
      if (!all(c(c1+c3, c2+c4) == 3)) {
        stop ("recombinant chromatids not complementary")
      }

      # Independent orientation +++++++++++ This is where I've got to +++++++++++++
      recomb_chromatid_ids <- sapply(1:14, function(chr) sample(x = 1:4, size = 4))

      # Sample children genotypes dependently
      clones <- TRUE
      while (clones) {

        # Crossing over
        # Assume both pairs of homologous chromosomes crosses over
        A1 <- sample(0:1, replace = T, size = 14) + 1
        B1 <- abs(A1-1) + 1
        A2 <- sample(0:1, replace = T, size = 14) + 1
        B2 <- abs(A2-1) + 1

        # Independent orientation of chromosomes
        # ++++++++++++++++++++
        # This is where I'm at: need to map tp 150 markers
        sapply(1:14, function(i) sample(1:4, size = 4))


        if (!all(rowSums(cbind(x, x_c)) == 3)) stop("Not complementary")

        child1 <- sapply(all_markers, function(t) sample(c(parent1[[t]], parent2[[t]]), 1), simplify = F)
        child1c <-  sapply(all_markers, function(t) sample(c(parent1[[t]], parent3[[t]]), 1), simplify = F)
        child2 <-  sapply(all_markers, function(t) sample(c(parent2[[t]], parent3[[t]]), 1), simplify = F)
        child2c

        clones <- identical(child13[1:m_min_clone], child23[1:m_min_clone])
      }

      # Make parasite infection and data
      initial <- rbind(unlist(child12))
      relapse <- rbind(unlist(child13), unlist(child23)) # MOI recurrence > initial
      y <- list(initial = apply(initial, 2, unique, simplify = F),
                relapse = apply(relapse, 2, unique, simplify = F))

      # Store the data
      ys_store[[as.character(c)]][[i]] <- y


  }
}
tictoc::toc()

#===============================================================================
# Generate results for different allele frequency types, for a migrant parent
# (rare_enrich TRUE) as well as parents from a single population, and for
# different marker counts
#===============================================================================
tictoc::tic()
for(c in c_params) {
  fs <- fs_store[[as.character(c)]]
    for(i in 1:n_repeats){
      y_all_markers <- ys_store[[as.character(c)]][[sprintf("rare_enrich_%s", rare_enrich)]][[i]]
      for(m in n_markers){
        y <- sapply(y_all_markers, function(x) x[1:m], simplify = FALSE)
        ps <- compute_posterior(y, fs, return.RG = TRUE)
        ps_store[[as.character(c)]][[sprintf("rare_enrich_%s", rare_enrich)]][[as.character(i)]][[as.character(m)]] <- ps
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
  y_all_markers <- ys_store[[as.character(c)]][[sprintf("rare_enrich_%s", rare_enrich)]][[i]]
  # compute posterior relapse probabilities for all marker counts from one onwards
  for(m in 1:max(n_markers)){
    y <- sapply(y_all_markers, function(x) x[1:m], simplify = FALSE)
    ps <- compute_posterior(y, fs)
    ps_store_all_ms[[as.character(i)]][[paste0("m",m)]] <- ps$marg[,"L"]
  }
}
tictoc::toc()


#===============================================================================
# Bundle data, results and magic numbers
#===============================================================================
Mieotic_siblings <- list(fs_store = fs_store,
                      ys_store = ys_store,
                      ps_store = ps_store,
                      ps_store_all_ms = ps_store_all_ms,
                      c_cutoff = c_cutoff,
                      n_alleles = n_alleles,
                      m_min_clone = m_min_clone,
                      seed = seed)


#===============================================================================
# Save Half_siblings as exported data
#===============================================================================
usethis::use_data(Mieotic_siblings, overwrite = TRUE)
