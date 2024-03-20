# NOTE: from deprecated model

# Script written to check if expected relatedness is 1/3 between meiotic products
set.seed(1)

sims <- sapply(1:100, function(i) {

  # Make homologous chromosomes: one from mum, one from dad
  h1 <- rep("m", 14)
  h3 <- rep("d", 14);

  # Make Sister chromatids (S-phase of interphase)
  h2 <- h1
  h4 <- h3

  # On average we expect one cross over per sister chromatid
  # Sample occurrence before some marker of interest
  crossed_over <- sample(c(T,F), 14, replace=T)
  h2[crossed_over] <- "d"
  h4[crossed_over] <- "m"

  # Align
  tetra <- cbind(h1,h2,h3,h4)

  # Independent orientation
  tetra_p <- t(apply(tetra, 1, function(x) return(x[sample(4,4)])))

  # Compute relatedness between each pair of haploid products
  all_hpairs <- gtools::combinations(4,2)
  rs <- apply(all_hpairs, 1, function(pair) {
    mean(tetra_p[,pair[1]] == tetra_p[,pair[2]])})

  return(rs)
})

hist(as.vector(sims))
abline(v = 0.33)
