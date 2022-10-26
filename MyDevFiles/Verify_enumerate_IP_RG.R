# original approach enumerates all possible IPs and calculates the probability of each given RG
# new approach takes an RG, and enumerates all IPs consistent with the RG
# claim: all consistent IPs have equal probability

# note enumerate_IPs_RG only works with RGs obtained from enumerate_RGs_alt,
# does not work with RGs obtained from enumerate_RGs

library(Pv3Rs)

MOIs <- c(2,2)

# igraph must be true for compute_pr_IP_RG to work
RGs <- enumerate_RGs_alt(MOIs, igraph=T)

# only needed for original approach
IPs <- enumerate_IPs(sum(MOIs))
IPs <- lapply(IPs, eq_canonical) # makes all IPs sorted for easier comparison

for(RG in RGs) {
  # get all probabilities using original approach
  IP_prs <- sapply(IPs, compute_pr_IP_RG, RG)

  # check all positive probabilities are equal
  stopifnot(length(unique(IP_prs[IP_prs>0])) == 1)

  # new approach enumerates IPs consistent with given RG
  IPs_RG <- enumerate_IPs_RG(RG)
  # makes all IPs sorted for easier comparison
  IPs_RG <- lapply(IPs_RG, eq_canonical)

  # check number of IPs enumerated is correct
  stopifnot(sum(IP_prs>0) == length(IPs_RG))

  # check that IPs enumerated are correct
  IP_withpos <- IPs[IP_prs>0]
  for(IP in IPs_RG) {
    # this line relies on eq_canonical being run first
    stopifnot(sum(sapply(IP_withpos, function(x) identical(x, IP))) == 1)
  }
}
writeLines("\nTests for correctness passed")

# The following compares enumeration of RGs
library(tictoc)

# time comparison on slightly larger example, 250 RGs
MOIs <- c(2,1,2)

# original approach of enumerating RGs
tic("\nOriginal enumeration of RGs with MOIs = (2, 1, 2)")
RGs_ori <- enumerate_RGs(MOIs)
toc()

# new approach of enumerating RGs
tic("\nNew enumeration of RGs with MOIs = (2, 1, 2)")
RGs <- enumerate_RGs_alt(MOIs, igraph=F)
toc()

# original approach
tic("\nCalculating p(IBD|RG) for every IBD with original approach")
for(RG in RGs_ori) {
  IP_prs <- sapply(IPs, compute_pr_IP_RG, RG)
}
toc()

# new approach
tic("\nFinding IBDs such that p(IBD|RG) > 0 with new approach")
for(RG in RGs) {
  enumerate_IPs_RG(RG)
}
toc()
