################################################################################
# Code to prepare Half_siblings dataset following best practice outlined in
# 7.1.1 of https://r-pkgs.org/data.html. This script only features in the source
# version of the package (it is listed under .Rbuildignore). For XX c_params
# and XX n_markers it takes XXX seconds to run on MacBook Pro with 32 GB memory
################################################################################
rm(list = ls())
library(MCMCpack) # For rdirichlet
library(tictoc) # For timing
set.seed(2)

c_params <- c(0.1, 1, 10, 100) # Dirichlet concentration parameter
c_cutoff <- 99 # Above c_cutoff, switch from Dirichlet to 1/n_alleles
n_markers <- c(10, 50, 100, 150) # Number of makers
n_alleles <- 5 # Number of alleles per marker (marker cardinality)
n_repeats <- 16 # Number of simulations per c_param, n_marker combination

tictoc::tic()
# If rare_enrich = TRUE, draw rare alleles with high probability for one of the parent (could
# interpret as an migrant from another population); otherwise, draw alleles
# proportional to allele frequencies for all parents
for(rare_enrich in c(TRUE, FALSE)) {

  ys <- list()
  ps <- list()
  ys_store <- list() # y for data
  ps_store <- list() # p for posterior
  fs_store <- list() # f for frequency
  rs_store <- list() # r for rare enrichment

  for(c in c_params) {
    for(m in n_markers) {

      markers <- paste0("m", 1:m) # Marker names
      alleles <- letters[1:n_alleles] # Alleles

      # Sample allele frequencies
      fs <- sapply(markers, function(m) {
        if(c > c_cutoff) {
          fs_unnamed <- rep(1/n_alleles, n_alleles)
        } else {
          fs_unnamed <- MCMCpack::rdirichlet(1, alpha = rep(c, n_alleles))
        }
        setNames(fs_unnamed, alleles)
      }, USE.NAMES = TRUE, simplify = FALSE)

      for(i in 1:n_repeats) {
        print(paste(rare_enrich,c,m,i))

        # Sample parental genotypes (ensure no clones)
        anyclones <- TRUE
        while (anyclones) {
          parent1 <- sapply(markers, function(t) {
            sample(alleles, size = 1, prob = fs[[t]])}, simplify = F)
          parent2 <- sapply(markers, function(t) {
            sample(alleles, size = 1, prob = fs[[t]])}, simplify = F)

          if (rare_enrich) {
            parent3 <- sapply(markers, function(t) {
              sample(alleles, size = 1, prob = 1-fs[[t]])}, simplify = F)
          } else {
            parent3 <- sapply(markers, function(t) {
              sample(alleles, size = 1, prob = fs[[t]])}, simplify = F)
          }
          anyclones <- any(identical(parent1, parent2),
                           identical(parent2, parent3),
                           identical(parent1, parent3))
        }


        # Sample children genotypes (ensure no clones)
        anyclones <- TRUE
        counter <- 0
        while (anyclones) {
          child12 <-  sapply(markers, function(t) sample(c(parent1[[t]], parent2[[t]]), 1), simplify = F)
          child13 <-  sapply(markers, function(t) sample(c(parent1[[t]], parent3[[t]]), 1), simplify = F)
          child23 <-  sapply(markers, function(t) sample(c(parent2[[t]], parent3[[t]]), 1), simplify = F)
          anyclones <- any(identical(child12, child13),
                           identical(child12, child23),
                           identical(child13, child23))
          counter <- counter + 1
          print(counter)
        }

        # Make parasite infection (incompatible with recrudescence) and data
        initial <- rbind(unlist(child12))
        relapse <- rbind(unlist(child13), unlist(child23)) # MOI recurrence > initial
        y <- list(initial = apply(initial, 2, unique, simplify = F),
                  relapse = apply(relapse, 2, unique, simplify = F))

        # Compute posterior
        ps[[i]] <- compute_posterior(y, fs, return.RG = TRUE)

        # Store the data
        ys[[i]] <- y

      }
      ys_store[[as.character(c)]][[as.character(m)]] <- ys
      ps_store[[as.character(c)]][[as.character(m)]] <- ps
      fs_store[[as.character(c)]][[as.character(m)]] <- fs
    }
  }

  rs_store[[sprintf("rare_enrich_%s", rare_enrich)]] <- list(ys_store = ys_store,
                                                            ps_store = ps_store,
                                                            fs_store = fs_store)
}
tictoc::toc()

Half_siblings <- rs_store

# Save Half_siblings as exported data
usethis::use_data(Half_siblings, overwrite = TRUE)
