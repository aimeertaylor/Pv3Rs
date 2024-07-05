################################################################################
# Script illustrating that the bulk data from a pool of three recombinant filial
# chromatids is indistinguishable from the bulk data from their two parents
#
# See following for examples of complex chiasmata
# https://biology.stackexchange.com/questions/42288/do-only-one-or-both-pairs-of-homologous-chromatids-exchange-genetic-material-dur
################################################################################
rm(list = ls())
library(Pv3Rs)
library(MCMCpack)
max_n_markers <- 150
all_markers <- paste0("m", 1:max_n_markers)
n_alleles <- 5
alleles <- letters[1:n_alleles]
fs_param <- setNames(rep(1, n_alleles), alleles)
fs <- sapply(all_markers, function(i) setNames(rdirichlet(1, fs_param), alleles), simplify = F)
parent1 <- sapply(all_markers, function(i) sample(alleles, size = 1, prob = fs[[i]]))
parent2 <- sapply(all_markers, function(i) sample(alleles, size = 1, prob = fs[[i]]))
parents <- cbind(parent1, parent2)
chrs_per_marker <- round(seq(0.51, 14.5, length.out = max_n_markers))
markers_per_chr <- table(chrs_per_marker)
cs_meiotic <- recombine_parent_ids(markers_per_chr)
cs_full <- sapply(1:4, function(i) recombine_parent_ids(markers_per_chr)[,1])

# Construct children genotypes from parental allocations
children_meiotic <- sapply(1:max_n_markers, function(i) {
  sapply(1:ncol(cs_meiotic), function(j) parents[i,cs_meiotic[i,j]])
})
colnames(children_meiotic) <- all_markers

# Construct children genotypes from parental allocations
children_full <- sapply(1:max_n_markers, function(i) {
  sapply(1:ncol(cs_full), function(j) parents[i,cs_full[i,j]])
})
colnames(children_full) <- all_markers

# Format parasite infection data for compute_posterior
y_full <- list(initial = apply(children_full[1:3,], 2, unique, simplify = F),
                   relapse = apply(children_full[4,,drop=FALSE], 2, unique, simplify = F))

y_meitotic <- list(initial = apply(children_meiotic[1:3,], 2, unique, simplify = F),
                   relapse = apply(children_full[4,,drop=FALSE], 2, unique, simplify = F))

y_parents <- list(initial = apply(t(parents), 2, unique, simplify = F),
                  relapse = apply(children_full[4,,drop=FALSE], 2, unique, simplify = F))

# Check they give the same posterior state probabilities:
suppressMessages(compute_posterior(y = y_meitotic, fs)$marg)
suppressMessages(compute_posterior(y = y_parents, fs)$marg)

# And that the problem does not occur given full siblings:
suppressMessages(compute_posterior(y = y_full, fs)$marg)

# Plot data
plot_data(list(meiotic = y_meitotic,
               parents = y_parents,
               full = y_full), marker_annotate = F)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# To-finish: compute IBD-based relatedness between meiotic sibs
all_pairs <- gtools::combinations(4,2)
apply(all_pairs, 1, function(i) {
  pair <- cs_meiotic[,all_pairs[i,]]
  mean(pair[,1] == pair[,2])
})

