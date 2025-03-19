################################################################################
# Investigate equivalence classes for example in `demonstrate-usage.Rmd`
################################################################################

library(Pv3Rs)
library(dplyr) # for `near`

y <- list("Enrollment" = list(m1 = c('C','G','T'),
                              m2 = c('A','C'),
                              m3 = c('C','G','T')),
          "Recurrence 1" = list(m1 = c('C','T'),
                                m2 = c('A'),
                                m3 = c('A','C')),
          "Recurrence 2" = list(m1 = c('T'),
                                m2 = c('A'),
                                m3 = c('A')))

fs <- list(m1 = c(A = 0.27, C = 0.35, G = 0.18, T = 0.20),
           m2 = c(A = 0.78, C = 0.14, G = 0.07, T = 0.01),
           m3 = c(A = 0.21, C = 0.45, G = 0.26, T = 0.08))

post <- compute_posterior(y, fs, return.RG = TRUE, return.logp = TRUE)
lliks <- sapply(post$RGs, function(RG) RG$logp)
gs <- paste0("g", 1:6)
ts_per_gs <- rep(1:length(y), determine_MOIs(y))

## How many graphs have the same logl?
sorted_lliks <- sort(lliks, decreasing = T)
plot(sorted_lliks[1:50])
# are the (1st, 2nd, 3rd, ...) logls are equal to (2nd, 3rd, 4th, ...) logls?
# adj_equal is FALSE at the indices where the sorted log-likelihoods decrease
adj_equal <- near(head(sorted_lliks, -1), tail(sorted_lliks, -1))
decr_idxs <- which(adj_equal == FALSE) # 2, 8, 14, 20, 32, ...
# first 2 have same logl, next 6 have the same logl etc.
class_sizes <- c(decr_idxs[1], diff(decr_idxs))

# Warning: In the code, RGs with the same logl are assumed to be in the same
# equivalent class. Technically, it is possible for two non-equivalent RGs (in
# the genotype  permutation symmetry sense) to have the same logl.

## Which EC has the highest probability?

# logl of a representative from each 'equivalence class' (EC)
lliks_unique <- sorted_lliks[decr_idxs]
# prob of EC
class_ps <- exp(lliks_unique)*class_sizes
# which EC has highest probability
max_class_p <- which(class_ps == max(class_ps)) # 5
max_idx <- decr_idxs[max_class_p] # 32
max_size <- class_sizes[max_class_p] # 12
# this EC consists of the RGs with logl rank 21-32
RG_order <- order(lliks, decreasing = T) # order RGs by logl
for(i in (max_idx-max_size+1):max_idx) {
  RG <- post$RGs[[RG_order[i]]]
  seqs_comp_MLE_RG <- compatible_rstrs(RG, split(gs, ts_per_gs))
  par(mar = rep(0.1,4))
  plot_RG(RG_to_igraph(RG, gs, ts_per_gs), edge.curved=0.25, vertex.size=20)
}
# Compared to the RG with highest logl, these RGs have more sibling edges,
# which intuitively seems to support the data better given the genetic
# similarity between the first two episodes. Interestingly, this RG is not
# supported by the MLE recurrence sequence, IC.
