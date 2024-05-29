################################################################################
# Script illustrating that the bulk data from a pool of four recombinant filial
# chromatids is indistinguishable from the bulk data from their two parents
#
# See following for examples of complex chiasmata
# https://biology.stackexchange.com/questions/42288/do-only-one-or-both-pairs-of-homologous-chromatids-exchange-genetic-material-dur
################################################################################

library(Pv3Rs) # For plot_data
rm(list = ls()) # For cleanliness
set.seed(1) # For reproducibility

chr_length <- 20 # Length of the chromosome
n_recom_events <- 2 # Number of independent recombination events

# Sample parents
parentA <- sample(c(0,1), replace = T, size = chr_length)
parentB <- sample(c(0,1), replace = T, size = chr_length)

# Record sister chromatids pre crossover
pre_CO <- cbind(parentA, parentA, parentB, parentB)

# For i in 1 to some number of recombination events
post_COs <- lapply(1:n_recom_events, function(i) {

  # Sample number of crossovers
  # Although there could be more than two crossovers, let's stick to at most two
  n_cross_over <- sample(c(1,2), 1)

  # Sample chiasmata position
  # Although a cross-over could technically occur at the same position on
  # complementary chromatids, let's avoid sampling chiasmata positions with
  # replacement
  if (n_cross_over > chr_length) stop # check we can sample without replacement
  chiasmata_pos <- sample(1:(chr_length-1), replace = F, n_cross_over)

  # Sample sister chromatids per chiasmata position
  # Although two cross-overs could involve the same sister chromatid, let's assume
  # at most one crossover per chromatid
  if(length(chiasmata_pos) == 2) {

    post_CO <- cbind(c(parentA[1:chiasmata_pos[1]], parentB[(chiasmata_pos[1]+1):chr_length]),
                     c(parentA[1:chiasmata_pos[2]], parentB[(chiasmata_pos[2]+1):chr_length]),
                     c(parentB[1:chiasmata_pos[1]], parentA[(chiasmata_pos[1]+1):chr_length]),
                     c(parentB[1:chiasmata_pos[2]], parentA[(chiasmata_pos[2]+1):chr_length]))

  } else if(length(chiasmata_pos) == 1) {

    post_CO <- cbind(c(parentA[1:chiasmata_pos[1]], parentB[(chiasmata_pos[1]+1):chr_length]),
                     parentA,
                     c(parentB[1:chiasmata_pos[1]], parentA[(chiasmata_pos[1]+1):chr_length]),
                     parentB)
  }

  return(post_CO)
})

# Check the post_COs are different
identity_matrix <- sapply(1:(length(post_COs)-1), function(i) {
  sapply((i+1):length(post_COs), function(j) {
    identical(post_COs[[i]], post_COs[[j]])
  })
})

# Some of the recombination events are the same
# (not so surprising given the chr_length is small)
sapply(identity_matrix, function(x) sum(as.numeric(x)))

# Reduce pool of recombinant chromatids to bulk data
ys_post <- lapply(post_COs, function(post_CO) {
  x <- apply(post_CO, 1, unique, simplify = F)
  names(x) <- paste0("m", 1:chr_length)
  return(x)
})
names(ys_post) <- paste0("recombin", 1:100)

# Reduce chromatids pre recombination to bulk data
ys_pre <- apply(pre_CO, 1, unique, simplify = F)
names(ys_pre) <- paste0("m", 1:chr_length)

# Format for plot
pre <- list(ys_pre = ys_pre)
ys <- c(pre, ys_post)
x <- list(ys = ys)

# Plot data
plot_data(x, marker_annotate = F)





