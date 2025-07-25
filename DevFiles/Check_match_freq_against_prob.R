# NOTE: from deprecated model

# Check match frequency given sampled lineages is equal to match
# probability given lineage frequencies

# Generate
n <- 40
lineage_probs <- as.vector(gtools::rdirichlet(1, alpha = rep(0.1, n)))
names(lineage_probs) <- generate_lineages(n)
f_comp <- sum(lineage_probs^2)
g_comp <- 0.5*(f_comp + 1)


# Function to sample lineages given sibling pairs
sim_sibling_cluster_counts <- function (n_sim, pop_lineages, rg_size) {
  sim_res <- sapply(1:n_sim, function(i){
    parents <- sample(pop_lineages$lineages,
                      size = 2,
                      replace = T,
                      prob = pop_lineages$probs)
    rg_lineages <- sample(parents, rg_size, replace = T)
    length(unique(rg_lineages))})
}


# Simulate lots of stranger match frequencies
f_sims <- sapply(1:100, function(i) {
  str_match <- sapply(1:100, function(j){
    lineage_pair <- sample(names(lineage_probs),
                           size = 2,
                           replace = T,
                           prob = lineage_probs)
    length(unique(lineage_pair)) == 1
  })
  mean(str_match)
})

# Simulate lots of sibling match frequencies
g_sims <- sapply(1:100, function(i) {
  sib_match <- sapply(1:100, function(j){
    parents <- sample(names(lineage_probs),
                      size = 2,
                      replace = T,
                      prob = lineage_probs)
    lineage_pair <- sample(parents, 2, replace = T)
    length(unique(lineage_pair)) == 1
  })
  mean(sib_match)
})


# Plot
par(mfrow = c(2,1))
H <- hist(f_sims, main = "", xlab = "Stranger match frequency")
abline(v = f_comp, col = "red")
text(x = f_comp, y = max(H$counts), pos = 4, labels = "Stranger match probability",
     col = "red", cex = 0.5)
H <- hist(g_sims, main = "", xlab = "Sibling match frequency")
abline(v = g_comp, col = "red")
text(x = g_comp, y = max(H$counts), pos = 4, labels = "Sibling match probability",
     col = "red", cex = 0.5)

