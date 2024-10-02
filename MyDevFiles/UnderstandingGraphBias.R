################################################################################
# To update and make into a markdown potentially
#
# The output of compute posterior converges approximately to the prior-weighted
# relative fraction of relationship graphs compatible with the observed data.
################################################################################
rm(list = ls())

# Simplex plot parameters
pardefault <- par()
par(mar = c(0,0,0,0))

# Define function to get gs, RGs, states given RGs, and thus avoid repetition
get_graph_dist <- function(x) {
  gs <- paste0("g", 1:sum(x)) # name genotypes
  ts <- 1:length(x) # episode indices
  ts_per_gs <- rep(ts, x) # episode index of each genotype
  gs_per_ts <- split(gs, ts_per_gs) # genotypes grouped by episode
  RGs <- enumerate_RGs(x) # generate all relationship graphs
  CIL_gvn_RGs <- sapply(RGs, compatible_rstrs, gs_per_ts) # list compatible states
  return(list(gs = gs, RGs = RGs, CIL_gvn_RGs = CIL_gvn_RGs))
}




#===============================================================================
# Homoallelic example without recurrent data using MOIs to change graphs.
#
# Homoallelic call limits summation of RGs to RGs with intra-episode
# relatedness.
# ===============================================================================
fs = list(m1 = c("A" = 0.001, "B" = 0.999)) # Allele frequencies
y <- list(enroll = list(m1 = c('A')), recur = list(m1 = NA)) # Data
MOIs <- list(c(2,1), c(3,1), c(2,2), c(3,2), c(3,3)) # MOIs
results <- sapply(MOIs, function(x) compute_posterior(y, fs, MOIs = x)$marg)
plot_simplex(c("Recrudescence", "Relapse", "Reinfection"), 0.5) # Plot simplex
xy <- apply(results, 2, project2D) # Project probabilities
points(x = xy["x", ], y = xy["y", ], pch = 20, cex = 1, col = 1:length(MOIs))
legend("left", col = 1:length(MOIs), pch = 20, title = "MOIs", inset = 0,
       legend = sapply(MOIs, paste, collapse = " & "), bty = "n")

# ------------------------------------------------------------------------------
# Posterior approximation based on graphs only
# ------------------------------------------------------------------------------
graph_frac_hom <- sapply(1:length(MOIs), FUN = function(i) {

  stuff <- get_graph_dist(MOIs[[i]])
  RGs_C <- sapply(stuff$CIL_gvn_RGs, function(RG) "C" %in% RG)
  RGs_I <- sapply(stuff$CIL_gvn_RGs, function(RG) "I" %in% RG)

  # Identify RGs with intra-episode relatedness for the initial episode
  if(identical(MOIs[[i]], c(1,1))) {
    RGs_edges_first <- rep(TRUE, length(stuff$RGs)) # monoclonal case is trivial
  } else { # Make adjacency matrix for the first episode
    MOI1 <- MOIs[[i]][1]
    gs1 <- stuff$gs[1:MOI1]
    mat_first <- array(1, dim = c(MOI1, MOI1), dimnames = list(gs1, gs1))
    edges_first <- igraph::as_ids(igraph::E(igraph::graph_from_adjacency_matrix(mat_first, mode = "undirected", diag = F)))
    RGs_edges <- lapply(stuff$RGs, function(RG) igraph::as_ids(igraph::E(RG))) # Get edges for every RG
    RGs_edges_first <- sapply(RGs_edges, function(RG_E) all(edges_first %in% RG_E))
  }

  graph_frac_unnormalised <- c(frac_C = sum(RGs_edges_first*RGs_C)/sum(RGs_C),
                               frac_L = sum(RGs_edges_first)/length(stuff$RGs),
                               frac_I = sum(RGs_edges_first*RGs_I)/sum(RGs_I))
  graph_frac <- graph_frac_unnormalised/sum(graph_frac_unnormalised)
  return(graph_frac)
})
xy_graph_frac_hom <- apply(graph_frac_hom, 2, project2D)
points(x = xy_graph_frac_hom["x", ], y = xy_graph_frac_hom["y", ],
       col = 1:length(MOIs), cex = 2)
legend("bottom", pch = c(20, 1), pt.cex = c(1, 2), # Add legends
       legend = c("Posterior probability", "Relative weighted graphs"), bty = "n")

# ------------------------------------------------------------------------------
# Understanding approximation based on graphs only
# ------------------------------------------------------------------------------
MOIs <- c(3,3) # Illustrative example
stuff <- get_graph_dist(MOIs)
result33  <- compute_posterior(y, fs, MOIs = MOIs, return.logp = TRUE)
logps <- sapply(result33$RGs, function(RG) RG$logp) # extract logp

# Make a vector of intra-infection edges by first creating a block diag. matrix
mat_first <- array(1, dim = c(3, 3), dimnames = list(paste0("g", 1:3), paste0("g", 1:3)))
edges_first <- igraph::as_ids(igraph::E(igraph::graph_from_adjacency_matrix(mat_first, mode = "undirected", diag = F)))
RGs_edges <- lapply(stuff$RGs, function(RG) igraph::as_ids(igraph::E(RG))) # Get edges for every RG
RGs_edges_first <- sapply(RGs_edges, function(RG_E) all(edges_first %in% RG_E))

# Inspect likelihoods:
if (max(logps[!(RGs_edges_first)]) < min(logps[RGs_edges_first])) {
  writeLines("Expected-most-likely graphs are most likely ")
  #for(i in which(RGs_edges_first)) plot_RG(RGs[[i]])
  high_logps <- table(logps[RGs_edges_first])
  writeLines(sprintf("%s likelihood values among %s graphs", length(high_logps), sum(high_logps)))
  print(high_logps)
}























#===============================================================================
# Recrudescence vs relapse using MOIs to change graphs
# ===============================================================================
# ------------------------------------------------------------------------------
# Graphs with notable likelihood in the MOI c(3,3) case are limited to those
# with intra-infection relatedness compatible with recrudescence or relapse but
# not reinfection; there are multiple likelihood values among them.
# ------------------------------------------------------------------------------
MOIs <- c(3,3) # Illustrative example
gs <- paste0("g", 1:sum(MOIs))
ts <- 1:length(MOIs)
ts_per_gs <- rep(ts, MOIs)
gs_per_ts <- split(gs, ts_per_gs)
result33  <- compute_posterior(y, fs, MOIs = MOIs, return.logp = TRUE)
RGs <- sapply(result33$RGs, function(RG) RG_to_igraph(RG, gs, ts_per_gs))
CIL_gvn_RGs <- sapply(RGs, compatible_rstrs, gs_per_ts = gs_per_ts) # states
RGs_C_or_L <- sapply(CIL_gvn_RGs, function(CIL_gvn_RG) !("I" %in% CIL_gvn_RG))
RGs_C <- sapply(CIL_gvn_RGs, function(CIL_gvn_RG) "C" %in% CIL_gvn_RG)
logps <- sapply(result33$RGs, function(RG) RG$logp) # extract logp

# Make a vector of intra-infection edges by first creating a block diag. matrix
mat_within <- Matrix::bdiag(lapply(MOIs,function(x) matrix(1, ncol=x, nrow=x)))
colnames(mat_within) <- gs; rownames(mat_within) <- gs
edges_within <- igraph::as_ids(igraph::E(igraph::graph_from_adjacency_matrix(mat_within, mode = "undirected", diag = F)))
RGs_edges <- lapply(RGs, function(RG) igraph::as_ids(igraph::E(RG)))
RGs_edges_within <- sapply(RGs_edges, function(RG_E) all(edges_within %in% RG_E))


# Inspect likelihoods:
if (max(logps[!(RGs_C_or_L & RGs_edges_within)]) < min(logps[RGs_C_or_L & RGs_edges_within])) {
  writeLines("Expected-most-likely graphs are most likely ")
  for(i in which(RGs_C_or_L & RGs_edges_within)) plot_RG(RGs[[i]])
  high_logps <- table(logps[RGs_C & RGs_edges_within])
  writeLines(sprintf("%s likelihood values among %s graphs", length(high_logps), sum(high_logps)))
  print(high_logps)
}











