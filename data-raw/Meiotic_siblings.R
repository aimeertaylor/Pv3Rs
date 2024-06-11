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
seed <- 2 # For reproducibility
c_params <- c(0.5, 1, 100) # Dirichlet concentration parameter
n_markers <- c(10, 50, 100) # Number of makers
n_repeats <- 5 # Number of simulations per c_param, n_marker combination
c_cutoff <- 99 # Above c_cutoff, switch from Dirichlet to 1/n_alleles
n_alleles <- 4 # Number of alleles per marker (marker cardinality)
MOIs_per_infection <- c("2_1", "3_1")

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
marker_subsets[[as.character(min_n_markers)]] <- min_marker_subset
for(m in min_n_markers:(max_n_markers-1)) {
  markers_thusfar <- marker_subsets[[as.character(m)]]
  marker_to_add <- sample(setdiff(all_markers, markers_thusfar), 1)
  marker_subsets[[as.character(m+1)]] <- c(markers_thusfar, marker_to_add)
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

    # Sample parental genotypes
    # (ensure no clones and thus always the same no. of vertices in RGs)
    parental_clones <- TRUE
    while (parental_clones) {
      parent1 <- sapply(all_markers, function(t) {
        sample(alleles, size = 1, prob = fs[[t]])}, simplify = F)
      parent2 <- sapply(all_markers, function(t) {
        sample(alleles, size = 1, prob = fs[[t]])}, simplify = F)
      parental_clones <- identical(parent1[min_marker_subset], parent2[min_marker_subset])
    }
    parents <- cbind(parent1, parent2)

    # Note condition on sufficient diversity
    children_clones <- TRUE
    while (children_clones) {
      # Per chromosome parent1 segment length for recombinant chromatids one and two
      c1_p1_segment_length <- sapply(markers_per_chr, sample, size = 1)
      c2_p1_segment_length <- sapply(markers_per_chr, sample, size = 1)

      # Parent allocation for chromatid one after crossover
      c1 <- do.call("c", sapply(1:14, function(i) {
        x <- rep(2, markers_per_chr[i])
        x[1:c1_p1_segment_length[i]] <- 1
        return(x)
      }))

      # Parent allocation for chromatid two after crossover
      c2 <- do.call("c", sapply(1:14, function(i) {
        x <- rep(2, markers_per_chr[i])
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

      # Parent allocation for chromatid bundle pre independent orientation
      cs_pre <- cbind(c1,c2,c3,c4)

      # Independent orientation
      recomb_chromatid_ids <- sapply(1:14, function(chr) sample(x = 1:4, size = 4))

      # Parent allocation for chromatid bundle post independent orientation
      cs <- sapply(1:4, function(i) {
        c_ind <- recomb_chromatid_ids[i,]
        do.call("c", sapply(1:14, function(j) {
          m_ind <- which(chrs_per_marker == j)
          cs_pre[m_ind,c_ind[j]]
        }))
      })

      # Child genotypes
      children <- sapply(1:max_n_markers, function(i) {
        sapply(1:4, function(j) {
          parents[i,cs[i,j]]
        })
      })
      colnames(children) <- all_markers

      # Check diversity among children in the first few markers
      children_clones <- nrow(unique(children[,min_marker_subset])) < 4
    }

    # For different numbers of genotypes per infection
    for(MOIs in MOIs_per_infection) {

      # Make parasite infection and data
      # MOI initial > relapse (recrudescence plausible)
      if (MOIs == "2_1") {
        initial <- rbind(unlist(children[1,]),unlist(children[2,]))
        relapse <- rbind(unlist(children[3,]))
      } else if (MOIs == "3_1") {
        initial <- rbind(unlist(children[1,]),unlist(children[2,]),unlist(children[3,]))
        relapse <- rbind(unlist(children[4,]))
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
for(c in c_params) {
  fs <- fs_store[[as.character(c)]]
  for(i in 1:n_repeats){
    for(MOIs in MOIs_per_infection) {
      y_all_markers <- ys_store[[as.character(c)]][[MOIs]][[i]]
      for(m in n_markers){
        marker_subset <- marker_subsets[[as.character(m)]]
        y <- sapply(y_all_markers, function(x) x[marker_subset], simplify = FALSE)
        ps <- compute_posterior(y, fs, return.RG = TRUE)
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
tictoc::tic()
for(i in 1:n_repeats){
  for(MOIs in MOIs_per_infection) {
    y_all_markers <- ys_store[[as.character(c)]][[MOIs]][[i]]
    # compute posterior relapse probabilities for all marker counts from one onwards
    for(m in min_n_markers:max_n_markers){
      marker_subset <- marker_subsets[[as.character(m)]]
      y <- sapply(y_all_markers, function(x) x[1:m], simplify = FALSE)
      ps <- compute_posterior(y, fs)
      ps_store_all_ms[[MOIs]][[as.character(i)]][[paste0("m",m)]] <- ps$marg[,"L"]
    }
  }
}
tictoc::toc()


#===============================================================================
# Bundle data, results and magic numbers
#===============================================================================
Meiotic_siblings <- list(fs_store = fs_store,
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
usethis::use_data(Meiotic_siblings, overwrite = TRUE)
