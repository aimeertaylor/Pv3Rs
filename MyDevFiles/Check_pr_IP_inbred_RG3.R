###############################################
# Model checking
###############################################

# Set up
rm(list = ls())
set.seed(1)
n <- 5
rep_num <- 100 # Number of repeats for pvalues and histograms
sim_num <- 1000 # Number of simulations from which to compute frequencies
par(mfrow = c(2,2))

# If FALSE, compute_pr_partitions does not work; otherwise, works.
Equifrequent_lineages <- TRUE

# Generate n founder lineages for simulation
if(Equifrequent_lineages) {
  lineage_probs <- rep(1/n, n)
} else {
  lineage_probs <- as.vector(MCMCpack::rdirichlet(1, alpha = rep(0.1, n)))
}

names(lineage_probs) <- generate_lineages(n)
f <- sum(lineage_probs^2) # Compute probability of being IBD

g_count <- 3

# Generate all RGs compatible with three vertices across infections
RGs <- enumerate_RGs(c(1,1,1))
IPs <- enumerate_IPs(c(1,1,1))

# Plot RGs
par(mfrow = c(4,3), mar = c(0,0,0,0))
sapply(RGs, plot_RG, edge.curved = 0.5)

# By-hand computed probabilities of IP given RGs
# (I don't yet have a function to compute probabilities given inbreeding)
Pr_IPs_RGs <- cbind(RG1 = c("(g1,g2,g3)" = f^2,
                            "(g1,g3)(g2)" = f*(1-f),
                            "(g1,g2)(g3)" = f*(1-f),
                            "(g2,g3)(g1)" = f*(1-f),
                            "(g1)(g2)(g3)" = (1-f)*(1-2*f)),

                    RG2 = c("(g1,g2,g3)" = f^2 + 0.5*f*(1-f),
                            "(g1,g3)(g2)" = 0.5*f*(1-f),
                            "(g1,g2)(g3)" = (3/2)*f*(1-f) + 0.5*(1-f)*(1-2*f),
                            "(g2,g3)(g1)" = 0.5*f*(1-f),
                            "(g1)(g2)(g3)" =  0.5*(1-f)*(1-2*f)),

                    RG3 = c("(g1,g2,g3)" = f,
                            "(g1,g3)(g2)" = 0,
                            "(g1,g2)(g3)" = 1-f,
                            "(g2,g3)(g1)" = 0,
                            "(g1)(g2)(g3)" = 0),

                    RG4 = c("(g1,g2,g3)" = f^2 + 0.5*f*(1-f),
                            "(g1,g3)(g2)" = (3/2)*f*(1-f) + 0.5*(1-f)*(1-2*f),
                            "(g1,g2)(g3)" = 0.5*f*(1-f),
                            "(g2,g3)(g1)" = 0.5*f*(1-f),
                            "(g1)(g2)(g3)" =  0.5*(1-f)*(1-2*f)),

                    RG5 = c("(g1,g2,g3)" = f,
                            "(g1,g3)(g2)" = 1-f,
                            "(g1,g2)(g3)" = 0,
                            "(g2,g3)(g1)" = 0,
                            "(g1)(g2)(g3)" = 0),

                    RG6 = c("(g1,g2,g3)" = f^2 + 0.5*f*(1-f),
                            "(g1,g3)(g2)" = 0.5*f*(1-f),
                            "(g1,g2)(g3)" = 0.5*f*(1-f),
                            "(g2,g3)(g1)" = (3/2)*f*(1-f) + 0.5*(1-f)*(1-2*f),
                            "(g1)(g2)(g3)" =  0.5*(1-f)*(1-2*f)),

                    RG7 = c("(g1,g2,g3)" = 0.25*(1-f) + f,
                            "(g1,g3)(g2)" = 0.25*(1-f),
                            "(g1,g2)(g3)" = 0.25*(1-f),
                            "(g2,g3)(g1)" = 0.25*(1-f),
                            "(g1)(g2)(g3)" = 0),

                    RG8 = c("(g1,g2,g3)" = f + 0.5*(1-f),
                            "(g1,g3)(g2)" = 0,
                            "(g1,g2)(g3)" = 0.5*(1-f),
                            "(g2,g3)(g1)" = 0,
                            "(g1)(g2)(g3)" = 0),

                    RG9 = c("(g1,g2,g3)" = f + 0.5*(1-f),
                            "(g1,g3)(g2)" = 0.5*(1-f),
                            "(g1,g2)(g3)" = 0,
                            "(g2,g3)(g1)" = 0,
                            "(g1)(g2)(g3)" = 0),

                    RG10 = c("(g1,g2,g3)" = f,
                             "(g1,g3)(g2)" = 0,
                             "(g1,g2)(g3)" = 0,
                             "(g2,g3)(g1)" = 1-f,
                             "(g1)(g2)(g3)" = 0),

                    RG11 = c("(g1,g2,g3)" = f + 0.5*(1-f),
                             "(g1,g3)(g2)" = 0,
                             "(g1,g2)(g3)" = 0,
                             "(g2,g3)(g1)" = 0.5*(1-f),
                             "(g1)(g2)(g3)" = 0),

                    RG11 = c("(g1,g2,g3)" = 1,
                             "(g1,g3)(g2)" = 0,
                             "(g1,g2)(g3)" = 0,
                             "(g2,g3)(g1)" = 0,
                             "(g1)(g2)(g3)" = 0))

colSums(Pr_IPs_RGs) # check sum to one

# Generate store for pvalues
Fr_IP_RG_reps = array(0, dim = c(nrow(Pr_IPs_RGs), ncol(Pr_IPs_RGs), rep_num),
                      dimnames = list(rownames(Pr_IPs_RGs), colnames(Pr_IPs_RGs), NULL))
pvalue_store <- array(NA, dim = c(nrow(Pr_IPs_RGs), ncol(Pr_IPs_RGs)),
                      dimnames = list(rownames(Pr_IPs_RGs), colnames(Pr_IPs_RGs)))

# Compute frequency multiple times to plot histogram and compute pvalue
for(RG_i in 1:length(RGs)) {

  RG <- RGs[[RG_i]]
  pb <- txtProgressBar(min = 0, max = rep_num, style = 3)
  for(sim_i in 1:rep_num){
    IPs <- sapply(1:100, function(i) {
      IP_to_character(sample_IP_given_RG(RG, outbred = FALSE, lineage_probs = lineage_probs))
    })

    # Compute frequency of IP given RG and store
    Fr_IP_RG <- table(IPs)/100
    Fr_IP_RG_reps[names(Fr_IP_RG), RG_i, sim_i] <- Fr_IP_RG
    setTxtProgressBar(pb, sim_i)
  }
}

# Plot agreement between frequencies and probabilities
par(mfrow = c(3,2), mar = rep(3,4))
for(RG_i in 1:length(RGs)) {

  plot_RG(RGs[[RG_i]], edge.curved = 0.5)

  for(IP_name in rownames(Pr_IPs_RGs)){

    if((Pr_IPs_RGs[IP_name,RG_i] == 0) & all(Fr_IP_RG_reps[IP_name,RG_i,] == 0)) {
      plot.new()
      next()
    }

    # Compute and store pvalue
    pvalue <- mean(Fr_IP_RG_reps[IP_name,RG_i,] < Pr_IPs_RGs[IP_name,RG_i])
    pvalue_store[IP_name, RG_i] <- pvalue

    # Generate barplot using hist with breaks = rep_num-1
    H <- hist(Fr_IP_RG_reps[IP_name,RG_i,], yaxt = "n", xlab = "",
              main = paste(IP_name, signif(Pr_IPs_RGs[IP_name,RG_i],3)),
              las = 2,
              breaks = rep_num-1, xlim = c(0,1)) # bar per simulation => barplot
    abline(v = Pr_IPs_RGs[IP_name, RG_i],col = "red")
    text(x = Pr_IPs_RGs[IP_name, RG_i], y = max(H$counts),
         labels = paste('p =', round(pvalue, 3)), pos = 4, col = "red", cex = 1)
  }

  # Eyeball the difference between computed probabilities and frequencies
  print(cbind("freq" = rowMeans(Fr_IP_RG_reps[,RG_i,]),
        "prob" = Pr_IPs_RGs[,RG_i]))
}

