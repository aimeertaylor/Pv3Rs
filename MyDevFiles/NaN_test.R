library(MCMCpack) # For rdirichlet
set.seed(1)

c_params <- c(1) # rep(1/n_alleles, n_alleles)
n_markers <- c(1000) # 10000 is also fine, but takes longer
n_repeats <- 10

ys <- list()
ys_store <- list()
posteriors <- list()
posteriors_store <- list()
fs_store <- list()

for(c in c_params) {
  for(m in n_markers) {

    markers <- paste0("m", 1:m) # Marker names
    alleles <- letters # Alleles
    n_alleles <- length(alleles) # Number of alleles per marker

    # Sample allele frequencies
    fs <- sapply(markers, function(m) {
      if(c > 999) {
        fs_unnamed <- rep(1/n_alleles, n_alleles)
      } else {
        fs_unnamed <- rdirichlet(1, alpha = rep(c, n_alleles))
      }
      setNames(fs_unnamed, alleles)
    }, USE.NAMES = TRUE, simplify = FALSE)

    for(i in 1:n_repeats) {

      # Sample parental genotypes
      parent1 <- sapply(markers, function(t) {
        sample(alleles, size = 1, prob = 1-fs[[t]])}, simplify = F)
      parent2 <- sapply(markers, function(t) {
        sample(alleles, size = 1, prob = 1-fs[[t]])}, simplify = F)

      # Sample children genotypes (ensure all different, s.t. we can focus on a subset of RGs)
      anyclones <- TRUE
      while (anyclones) {
        child1 <-  sapply(markers, function(t) sample(c(parent1[[t]], parent2[[t]]), 1), simplify = F)
        child2 <-  sapply(markers, function(t) sample(c(parent1[[t]], parent2[[t]]), 1), simplify = F)
        child3 <-  sapply(markers, function(t) sample(c(parent1[[t]], parent2[[t]]), 1), simplify = F)
        anyclones <- any(identical(child1, child2),
                         identical(child2, child3),
                         identical(child1, child3))
      }

      # Make parasite infection (incompatible with recrudescence)
      initial <- rbind(unlist(child1))
      relapse <- rbind(unlist(child2), unlist(child3))
      y <- list(initial = apply(initial, 2, unique, simplify = F),
                relapse = apply(relapse, 2, unique, simplify = F))

      # Store the data
      ys[[i]] <- y

      # Compute posterior
      posteriors[[i]] <- compute_posterior(y, fs, return.RG = TRUE)
    }
    posteriors_store[[as.character(c)]][[as.character(m)]] <- posteriors
    ys_store[[as.character(c)]][[as.character(m)]] <- ys
    fs_store[[as.character(c)]][[as.character(m)]] <- fs
  }
}

for (m in n_markers) {
  print(m)
  for (p in posteriors_store[["1"]][[as.character(m)]]) {
    print(p$marg)
  }
}

post <- posteriors_store[["1"]][["92"]][[1]]

logps <- sapply(post$RGs, "[[", "logp")
