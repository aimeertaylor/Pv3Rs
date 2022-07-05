# Compare the probabilities of IBD partitions given RG and an outbred population
# and the frequencies of IBD partitions sampled given RG and an outbred population
# outbred population
sim_count <- 100
par(mfrow = c(2,2))

# First try some large random RG
set.seed(3)
RG <- sample_single_recurrence_RG(MOI_init = 5, cause = "reinfection")
IPs <- enumerate_IPs(igraph::vcount(RG))
pr_IP_RG <- sapply(IPs, compute_pr_IP_RG, RG)
names(pr_IP_RG) <- sapply(IPs, convert_IP_to_string)
# Third simulate the frequency of the IPs under the outbred simulation model
simulated_IPs <- sapply(1:100, function(j) convert_IP_to_string(sample_IPs_given_RG(RG)))
simulated_IP_fr <- table(simulated_IPs)/100

# Plot agreement and relationship graph
par(mfrow = c(1,2), pty = "s")
plot(x = pr_IP_RG[names(simulated_IP_fr)],
     y = as.vector(simulated_IP_fr),
     xlab = "Probability", ylab = "Frequency",
     xlim = c(0,1), ylim = c(0,1),
     pch = 4)
abline(a = 0, b = 1)
plot_RG(RG)

# Construct various MOIs summing to at most five (six takes ages and beyond six
# there are two many graphs to enumerate) over at most three infections
MOIs <- c(do.call(c, apply(partitions::parts(2), 2, list)),
          do.call(c, apply(partitions::parts(3), 2, list)),
          do.call(c, apply(partitions::restrictedparts(4,3), 2, list)),
          do.call(c, apply(partitions::restrictedparts(5,3), 2, list)))

# Remove zero-valued MOIs
MOIs <- lapply(MOIs, function(x) x[x>0])

for(m in MOIs){

  RGs <- enumerate_RGs(m)
  IPs <- enumerate_IPs(m)
  pr_IP_RG <- array(0, dim = c(length(RGs), length(IPs)),
                    dimnames = list(NULL, sapply(IPs, convert_IP_to_string)))
  fr_IP_RG <- pr_IP_RG

  for(i in 1:length(RGs)){

    # Extract RG and convert to a sparse matrix
    RG <- RGs[[i]]
    pr_IP_RG[i,] <- sapply(IPs, compute_pr_IP_RG, RG)

    # Third simulate the frequency of the IPs under the outbred simulation model
    simulated_IPs <- sapply(1:sim_count, function(j) convert_IP_to_string(sample_IPs_given_RG(RG)))
    simulated_IP_fr <- table(simulated_IPs)/sim_count
    fr_IP_RG[i, names(simulated_IP_fr)] <- simulated_IP_fr
  }

  all(rowSums(pr_IP_RG) == 1) # Check sum to one
  all(rowSums(fr_IP_RG) == 1) # Check sum to one

  Rmsq <- vector(length = nrow(pr_IP_RG))
  for(i in 1:nrow(pr_IP_RG)) {
    Rmsq[i] <- sqrt(mean((pr_IP_RG[i,] - fr_IP_RG[i,])^2))
  }

  plot(x = as.vector(pr_IP_RG), y = as.vector(fr_IP_RG),
       xlab = "Probability computed", ylab = "Frequency simulated",
       pch = 4, main = paste(m, collapse = ""))
  abline(a = 0, b = 1)
}


