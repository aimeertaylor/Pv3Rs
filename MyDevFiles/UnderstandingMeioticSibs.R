################################################################################
# Script illustrating that perfect bulk data from a pool of three meiotic
# siblings is indistinguishable from perfect bulk data from their two parents;
# that a single error is enough to recover full-sibling like behaviour; and that
# meiotic siblings are 33% related.
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
fs <- sapply(all_markers, function(i) {
  setNames(rdirichlet(1, fs_param), alleles)}, simplify = F)
parent1 <- sapply(all_markers, function(i) sample(alleles, size = 1, prob = fs[[i]]))
parent2 <- sapply(all_markers, function(i) sample(alleles, size = 1, prob = fs[[i]]))
parents <- cbind(parent1, parent2)
chrs_per_marker <- round(seq(0.51, 14.5, length.out = max_n_markers))
cs_meiotic <- recombine_parent_ids(chrs_per_marker)
cs_full <- sapply(1:4, function(i) recombine_parent_ids(chrs_per_marker)[,1])

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

# Meiotic sibs and their parents indistinguishable:
suppressMessages(compute_posterior(y = y_parents, fs)$marg)
suppressMessages(compute_posterior(y = y_meitotic, fs)$marg)

# Providing the true but unknowable MOI makes little difference (unknowable
# without single-cell sequences)
suppressMessages(compute_posterior(y = y_meitotic, fs, MOIs = c(3,1))$marg)

# Because we assume no error, a single error is enough to solve the problem:
y_meitotic$initial$m150 <- "a"
suppressMessages(compute_posterior(y = y_meitotic, fs)$marg)

# No problem with full siblings:
suppressMessages(compute_posterior(y = y_full, fs)$marg)

# Plot data:
plot_data(list(meiotic = y_meitotic,
               parents = y_parents,
               full = y_full), marker_annotate = F)

# Check siblings are related as expected:
all_pairs <- gtools::combinations(4,2)

average_meiotic_relatedness <- mean(apply(all_pairs, 1, function(ind){
  pair <- cs_meiotic[,ind]
  mean(pair[,1] == pair[,2])
}))

average_full_relatedness <- mean(apply(all_pairs, 1, function(ind){
  pair <- cs_full[,ind]
  mean(pair[,1] == pair[,2])
}))

average_meiotic_relatedness
average_full_relatedness
