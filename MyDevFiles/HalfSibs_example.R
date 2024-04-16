# ==============================================================================
# aside: underflow / overflow problems with large marker counts
# aside: how to plot a graph for a real sample where gs and ts_per_gs are unknown? - record
# Hypothesis: something to do with where the rare alleles
# Edge case: conc_param = 1, 50 markers
# Correct answer with the wrong reasoning.


# For different numbers of markers, different concentration parameters, over X repeats
# To-do: plot RGs
# Process llikeRGs
# ==============================================================================
rm(list = ls())
library(MCMCpack) # For rdirichlet
set.seed(1)

c_params <- c(1,10,1000) # rep(1/n_alleles, n_alleles)
n_markers <- c(10,30,50,70)
n_repeats <- 15

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
      parent3 <- sapply(markers, function(t) {
        sample(alleles, size = 1, prob = 1-fs[[t]])}, simplify = F)

      # Sample children genotypes (ensure all different, s.t. we can focus on a subset of RGs)
      anyclones <- TRUE
      while (anyclones) {
        child12 <-  sapply(markers, function(t) sample(c(parent1[[t]], parent2[[t]]), 1), simplify = F)
        child13 <-  sapply(markers, function(t) sample(c(parent1[[t]], parent3[[t]]), 1), simplify = F)
        child23 <-  sapply(markers, function(t) sample(c(parent2[[t]], parent3[[t]]), 1), simplify = F)
        anyclones <- any(identical(child12, child13),
                         identical(child12, child23),
                         identical(child13, child23))
      }

      # Make parasite infection (incompatible with recrudescence)
      initial <- rbind(unlist(child12))
      relapse <- rbind(unlist(child13), unlist(child23))
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


# Extract probability of relapse
post_L <- sapply(posteriors_store, function(X) {
  sapply(X, function(XX) {
    sapply(XX, function(XXX) XXX$marg[,"L"])
  }, simplify = F)
}, simplify = F)

# Plots posterior relapse probabilities
par(mfrow = c(length(c_params),1))
for(c in c_params){
  plot(NULL, xlim = range(n_markers)+c(-10,10), ylim = c(0,1),
       xaxt = "n", bty = "n", panel.first = grid(nx = NA, ny = NULL),
       ylab = "Posterior relapse probability",
       xlab = "Number of markers (with added jitter)",
       main = sprintf("Concentration parameter: %s", c))
  axis(side = 1, at = n_markers)
  for(m in n_markers) {
    points(y = post_L[[as.character(c)]][[as.character(m)]],
           x = jitter(rep(m, n_repeats), amount = 5), pch = 4)
  }
}



#===============================================================================
# Extract graphs (i.e., discard logp)
justRGs <- sapply(posteriors_store, function(X) {
  sapply(X, function(XX) {
    sapply(XX, function(post) {
      sapply(post$RGs, function(RG) c(RG$clone, RG$sib))
    }, simplify = F)
  }, simplify = F)
}, simplify = F)


# Check all the graphs are returned in the same order
justRG <- justRGs[[1]][[1]][[1]]
RGcheck <- sapply(n_markers, function(m) {
  sapply(c_params, function(c) {
    sapply(2:n_repeats, function(i) {
      identical(justRG, justRGs[[as.character(c)]][[as.character(m)]][[i]])
    })
  })
})

if (!all(RGcheck)) stop("graphs not returned in the same order")

# Extract probability of the data given the graph (check)
llikeRGs <- sapply(posteriors_store, function(X) {
  sapply(X, function(XX) {
    sapply(XX, function(post) {
      sapply(post$RGs, function(RG) RG$logp)
    }, simplify = F)
  }, simplify = F)
}, simplify = F)


# Plot data
for(c in c_params){
  for(m in 10){ #n_markers){
    ys <- ys_store[[as.character(c)]][[as.character(m)]]
    names(ys) <- 1:n_repeats
    plot_data(ys, fs = fs_store[[as.character(c)]][[as.character(m)]], marker_annotate = F)
  }
}


