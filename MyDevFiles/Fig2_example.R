# ==============================================================================
# Made up example, based on Figure in manuscript
# ==============================================================================

# Make parasite genotypes, where green-ish genotypes are yellow-blue mosaics;
# Alleles are named by their ancestor.
Yel_child <- c(m1 = "yellow", m2 = "yellow", m3 = "yellow", m4 = "yellow")
YelBlu_child_1 <- c(m1 = "yellow", m2 = "yellow", m3 = "blue", m4 = "yellow")
YelBlu_child_2 <- c(m1 = "yellow", m2 = "blue", m3 = "blue", m4 = "yellow")
YelBlu_child_3 <- c(m1 = "yellow", m2 = "blue", m3 = "yellow", m4 = "blue")
YelBlu_child_4 <- c(m1 = "blue", m2 = "yellow", m3 = "blue", m4 = "blue")
Pink <- c(m1 = "hotpink", m2 = "hotpink", m3 = "hotpink", m4 = "hotpink")
Red <- c(m1 = "red", m2 = "red", m3= "red", m4 = "red")

# Make parasite infections
initial <- rbind(Yel_child, YelBlu_child_1, YelBlu_child_2, YelBlu_child_3)
relapse <- rbind(YelBlu_child_4, Pink)
recrud <- rbind(YelBlu_child_4)
reinf <- rbind(Red)

# Make list of per-marker list of observed alleles
# The model will mistakenly think there are two not four genotypes in the initial
# infection
y <- list(initial = apply(initial, 2, unique, simplify = F),
          relapse = apply(relapse, 2, unique, simplify = F),
          recrud = apply(recrud, 2, unique, simplify = F),
          reinf = apply(reinf, 2, unique, simplify = F))


# Function to draw allele frequencies that sum to one
rdirichlet(1, alpha = rep(0.1,4))

set.seed(1)
fs <- list(
  m1 = setNames(rdirichlet(1, alpha = rep(0.1,4)), c("yellow", "hotpink", "red", "blue")),
  m2 = setNames(rdirichlet(1, alpha = rep(0.1,4)), c("yellow", "hotpink", "red", "blue")),
  m3 = setNames(rdirichlet(1, alpha = rep(0.1,4)), c("yellow", "hotpink", "red", "blue")),
  m4 = setNames(rdirichlet(1, alpha = rep(0.1,4)), c("yellow", "hotpink", "red", "blue"))
)

# compute
post <- compute_posterior(y, fs, return.RG=TRUE, return.logp=TRUE)

post$marg # Makes sense
sort(post$joint, decreasing = T) # Makes sense

# How to plot a graph for a real sample where gs and ts_per_gs are unknown? - record
# Most likely has a clonal edge
RGlogp <- sapply(post$RGs, function(RG) RG$logp) # log likelihoods of graphs
RG <- post$RGs[[which.max(RGlogp)]]
gs <- paste0("g", 1:6)
ts_per_gs <- c(1,1,2,2,3,4)
par(mar = c(0.5, 0.5, 0.5, 0.5))
plot_RG(RG_to_igraph(RG, gs, ts_per_gs), edge.curved=0.25, vertex.size=20, vertex_palette = "Blues")
seqs_comp_MLE_RG <- compatible_rstrs(RG, split(gs, ts_per_gs))
