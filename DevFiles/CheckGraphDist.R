# ==============================================================================
# P(G | S) is uniform
# P(G) = sum over states P(G | S)* P(S) is not uniform
# ==============================================================================
library(Pv3Rs)
rm(list = ls())

# Specify MOIs and prior state relative frequencies
MOIs <- c(2,1) # Script only works for MOIs c(1,1) and c(2,1)
freq <- c(L = 10, I = 10, C = 1)

prob <- freq / sum(freq)
RGs <- enumerate_RGs(MOIs)
n_RGs <- length(RGs)

# Compute plotting parameters
mat <- matrix(data = 1, nrow = n_RGs, ncol = n_RGs)
mat[n_RGs, 1:n_RGs] <- 2:(n_RGs+1)
layout(mat)
par(mar = c(1,1,1,1))

# Generate a list of genotypes per episodes
gs <- paste0("g", 1:sum(MOIs)) # genotype names (graph vertices)
ts <- 1:length(MOIs) # episode indices
ts_per_gs <- rep(ts, MOIs) # episode index of each genotype
gs_per_ts <- split(gs, ts_per_gs) # genotypes grouped by episode

# For each RG, extract compatible states
compatible_states <- sapply(RGs, Pv3Rs:::compatible_rstrs, gs_per_ts)

# Compute conditional prior probabilities
C_comp <- sapply(compatible_states, function(s) "C" %in% s) # zero or not
I_comp <- sapply(compatible_states, function(s) "I" %in% s) # zero or not
prob_G_C <- C_comp * 1/(sum(C_comp))
prob_G_I <- I_comp * 1/(sum(I_comp))
prob_G_L <- rep(1/n_RGs, n_RGs)

# Compute unconditional prior probabilities
prob_G <- (prob_G_L * prob["L"]) + (prob_G_I * prob["I"]) + (prob_G_C * prob["C"])

# Order graphs by compatibility
RGs_sorted <- sort.int(sapply(compatible_states, paste, collapse = ""), index.return = T)

# Plot distribution and graphs
barplot(prob_G[RGs_sorted$ix], axes = FALSE)
sapply(RGs_sorted$ix, function(i){
  plot_RG(RGs[[i]], vertex.label = NA)
  box(col = "grey")
})

# Checks
if (any((c(sum(prob_G_C),
           sum(prob_G_I),
           sum(prob_G_L),
           prob_C + prob_I + prob_L,
           sum(prob_G)) - 1) > .Machine$double.eps)) stop("Problem with probabilities summing to one")


# ==============================================================================
# More generally, distribution over transitive graphs in Fig B2 of Mehra G3 2025
# ==============================================================================
m <- 4 # Cluster size used
RGs <- enumerate_RGs(MOIs = rep(1, m), igraph = F)

# Plot all graphs
sapply(RGs, function(g){
  plot_RG(RG_to_igraph(g), vertex.label = NA, edge.curved = 0.5)
})

# Ignoring node identities (this is impefect: graphs 10 and 11 are conflated)
remove_node_ids <- function(g){
  paste(c(sort(sapply(g$clone, length)), sort(sapply(g$sib, length))),
        collapse = "")
}

x <- table(sapply(RGs, remove_node_ids))
x/sum(x) # n.b. graphs of type 10 and 11 are conflated

