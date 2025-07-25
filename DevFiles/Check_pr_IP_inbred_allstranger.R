# NOTE: from deprecated model

###############################################
# Model checking
# Conclusion: a single inbreeding value can be used to compute the probability
# of an IBD graph among stranger parasites if their lineages are equifrequent;
# otherwise,
# if lineages are not equifrequent, f under-epresents the number that can be sampled
# if relationships are not stranger, f over-represents the number that can be sampled
# Consider writing up as an example in Pv3Rs::compute_pr_partition
###############################################

# Set up
rm(list = ls())
set.seed(1)
n <- 5
rep_num <- 10 # Number of repeats for pvalues and histograms
sim_num <- 100 # Number of simulations from which to compute frequencies
par(mfrow = c(2,2))

# Copied from https://stackoverflow.com/questions/57233914/how-to-code-elementary-symmetric-polynomials-in-r
# Required for non-uniform frequency Birthday problem formula from Munford 1977
sympoly <- function(k, x) sum(combn(x, k, prod))

# If FALSE, compute_pr_partitions does not work; otherwise, works.
Equifrequent_lineages <- TRUE

# Generate n founder lineages for simulation
if(Equifrequent_lineages) {
  lineage_probs <- rep(1/n, n)
} else {
  lineage_probs <- as.vector(gtools::rdirichlet(1, alpha = rep(0.1, n)))
}

names(lineage_probs) <- generate_lineages(n)
f <- sum(lineage_probs^2) # Compute probability of being IBD

for(g_count in 2:6){

  # Generate relationship graph for sample_IPs_given_RG
  g_names <- paste0("g", 1:g_count)
  RM <- array(0, dim = rep(g_count, 2), dimnames = list(g_names, g_names))
  RG <- igraph::graph_from_adjacency_matrix(adjmatrix = RM, mode='undirected', diag = F)

  # Compute Pr(IG|RG)
  Pr_IG_RG <- compute_pr_partitions(g_count, f)

  sapply(partitions::listParts(g_count), compute_pr_IP_RG, RG)

  # Remain Pr(IG|RG) to compare with output of sample_IPs_given_RG
  names(Pr_IG_RG) <- sapply(names(Pr_IG_RG), function(x) {
    x <- as.numeric(strsplit(x, split = "")[[1]])
    names(x) <- paste0("g", 1:length(x))
    convert_IP_to_sting(convert_lineages_to_IP(x))
  })

  # Generate store for pvalues
  Fr_IG_RG_reps = array(0, dim = c(rep_num, length(Pr_IG_RG)), dimnames = list(NULL, names(Pr_IG_RG)))
  pvalue_store <- array(NA, dim = length(Pr_IG_RG), dimnames = list(names(Pr_IG_RG)))


  # Compute frequency multiple times to plot histogram and compute pvalue
  pb <- txtProgressBar(min = 0, max = rep_num, style = 3)
  for(sim_i in 1:rep_num){
    IPs <- sapply(1:100, function(i) {
      convert_IP_to_sting(sample_IPs_given_RG(RG, outbred = FALSE, lineage_probs = lineage_probs))
    })

    # Compute frequency of IP given RG and store
    Fr_IG_RG <- table(IPs)/100
    Fr_IG_RG_reps[sim_i, names(Fr_IG_RG)] <- Fr_IG_RG
    setTxtProgressBar(pb, sim_i)
  }


  # Plot agreement between frequencies and probabilities
  for(IG_name in names(Pr_IG_RG)){

    if((Pr_IG_RG[IG_name] == 0) & all(Fr_IG_RG_reps[,IG_name] == 0)) next()

    # Compute and store pvalue
    pvalue <- mean(Fr_IG_RG_reps[,IG_name] < Pr_IG_RG[IG_name])
    pvalue_store[IG_name] <- pvalue

    # Generate barplot using hist with breaks = rep_num-1
    H <- hist(Fr_IG_RG_reps[,IG_name], yaxt = "n", xlab = "",
              main = paste(IG_name,signif(Pr_IG_RG[IG_name],3)),
              las = 2,
              breaks = rep_num-1, xlim = c(0,1)) # bar per simulation => barplot
    abline(v = Pr_IG_RG[IG_name],col = "red")
    text(x = Pr_IG_RG[IG_name], y = max(H$counts),
         labels = paste('p =', round(pvalue, 3)), pos = 4, col = "red", cex = 1)
  }

  # Eyeball the difference between computed probabilities and frequencies
  print(rbind(colMeans(Fr_IG_RG_reps), Pr_IG_RG))

  # Compare with different calculations
  if(g_count < 5) { # Otherwise too combinatorially complicated
    print(c(tail(Pr_IG_RG, 1),
            factorial(g_count)*sympoly(g_count, rep(1/n,n)),  # Formula from Munford 1977
            factorial(n)/((n^g_count)*factorial(n-g_count)))) # Assuming equal lineage frequencies
  }
}
















