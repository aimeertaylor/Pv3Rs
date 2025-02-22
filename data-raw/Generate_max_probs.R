################################################################################
# Script to generate maximum probabilities for all MOI vectors summing to
# at-most eight. Takes 5 hours to run.
################################################################################
rm(list = ls())
set.seed(1)
n <- 10000
many_MOIs <- list(c(1,1), c(2,1), c(2,2), c(1,1,1), # Should all be indifferent
                  c(3,1), c(3,2), c(3,3), c(3,3,2),
                  c(4,1), c(4,2), c(4,3), c(4,3,1))

all_MOIs <- c(do.call(c, sapply(2:7, function(x) { # 254 different vectors of MOIs
  y <- gtools::permutations(n = 7, r = x, repeats.allowed = T)
  z <- y[rowSums(y) < 9, ]
  lapply(1:nrow(z), function(i) z[i,])
})), list(rep(1,8)))

# Function to get size of largest sibling clique
get_max_clique_size <- function(RG, intra_edges) {
  all_edges <- igraph::as_ids(igraph::E(RG))
  inter_edges <- setdiff(all_edges, intra_edges)
  if (length(inter_edges) > 0) {
    RG_intra <- igraph::delete_edges(RG, inter_edges)
  } else {
    RG_intra <- RG
  }
  return(igraph::clique_num(RG_intra))
}

# Function to sample a graph given a recurrence
# as.numeric(as.character()) needed to override convenience behaviour; see ?sample
sample_graph <- function(r, RGs_ind, RGs_C_ind, RGs_I_ind) {
  if (r == "L") {
    sample(RGs_ind, 1)
  } else if (r == "C") {
    as.numeric(sample(as.character(RGs_C_ind), 1))
  } else if (r == "I") {
    as.numeric(sample(as.character(RGs_I_ind), 1))
  }
}

tictoc::tic()
max_probs <- sapply(all_MOIs, function(MOIs){

  #=============================================================================
  # Get requisite information for maximum probability computation
  #=============================================================================
  gs <- paste0("g", 1:sum(MOIs)) # genotype names (graph vertices)
  ts <- 1:length(MOIs) # episode indices
  ts_per_gs <- rep(ts, MOIs) # episode index of each genotype
  gs_per_ts <- split(gs, ts_per_gs) # genotypes grouped by episode

  # Get all intra_edges by first creating block diag. matrix
  mat <- Matrix::bdiag(lapply(MOIs, function(x) matrix(1, ncol=x, nrow=x)))
  colnames(mat) <- gs; rownames(mat) <- gs
  graph <- igraph::graph_from_adjacency_matrix(mat, diag = F, mode = "undirected")
  intra_edges <- igraph::as_ids(igraph::E(graph))

  RGs <- enumerate_RGs(MOIs) # Get graphs
  CIL_gvn_RGs <- sapply(RGs, compatible_rstrs, gs_per_ts) # compatible states

  if (all(ts < 3)) {
    RGs_C_log <- sapply(CIL_gvn_RGs, function(x) "C" %in% x) # C compatible
    RGs_I_log <- sapply(CIL_gvn_RGs, function(x) "I" %in% x) # I compatible
  } else {
    First_CIL_gvn_RGs <- sapply(CIL_gvn_RGs, function(x) {
      do.call(rbind, strsplit(x, split = ""))[,1]
    })
    RGs_C_log <- sapply(First_CIL_gvn_RGs, function(x) "C" %in% x) # C compatible
    RGs_I_log <- sapply(First_CIL_gvn_RGs, function(x) "I" %in% x) # I compatible
  }

  # Get the size of the largest intra-episode sibling clique per graph
  max_clique_sizes <- sapply(RGs, get_max_clique_size, intra_edges = intra_edges)
  keep_log <- max_clique_sizes < 3
  nRGs_with <- length(RGs) # No. of relapse-compatible graphs summing over all
  nRGs_wout <- sum(keep_log) # Not summing over all
  RGs_with <- 1:nRGs_with
  RGs_wout <- which(keep_log)
  RGs_C_with <- which(RGs_C_log)
  RGs_I_with <- which(RGs_I_log)
  RGs_C_wout <- which(RGs_C_log & keep_log)
  RGs_I_wout <- which(RGs_I_log & keep_log)

  #=============================================================================
  # Compute theoretical maximum probabilities with and w/out summation over all
  #=============================================================================
  if(sum(RGs_C_log) == 0){
    theory_C_with <- 0
    theory_C_wout <- 0
  } else {
    theory_C_with <- nRGs_with/(sum(RGs_C_log) + nRGs_with)
    theory_C_wout <- nRGs_wout/(sum(RGs_C_log & keep_log) + nRGs_wout)
  }

  theory_I_with <- nRGs_with/(sum(RGs_I_log) + nRGs_with)
  theory_I_wout <- nRGs_wout/(sum(RGs_I_log & keep_log) + nRGs_wout)

  #=============================================================================
  # Compute numerical maximum probabilities with and w/out summation over all
  #=============================================================================
  if(sum(RGs_C_log) == 0) {
    CIL <- sample(c("I", "L"), n, replace = T) # sample recurrence
  } else {
    CIL <- sample(c("C", "I", "L"), n, replace = T) # sample recurrence
  }

  graphs_with <- sapply(CIL, sample_graph, RGs_with, RGs_C_with, RGs_I_with)
  graphs_wout <- sapply(CIL, sample_graph, RGs_wout, RGs_C_wout, RGs_I_wout)

  if(sum(RGs_C_log) == 0) {
    sim_C_with <- 0
    sim_C_wout <- 0
  } else {
    # simulation C with
    ind <- graphs_with %in% RGs_C_with # Among graphs compatible with C
    sim_C_with <- (table(names(graphs_with[ind]))/sum(ind))["C"]

    # simulation C with
    ind <- graphs_wout %in% RGs_C_wout # Among graphs compatible with C
    sim_C_wout <- (table(names(graphs_wout[ind]))/sum(ind))["C"]
  }

  # simulation I with
  ind <- graphs_with %in% RGs_I_with # Among graphs compatible with I
  sim_I_with <- (table(names(graphs_with[ind]))/sum(ind))["I"]

  # simulation I with
  ind <- graphs_wout %in% RGs_I_wout # Among graphs compatible with I
  sim_I_wout <- (table(names(graphs_wout[ind]))/sum(ind))["I"]

  return(c(theory_C_with = theory_C_with,
           sim_C_with = as.numeric(sim_C_with),
           theory_C_wout = theory_C_wout,
           sim_C_wout = as.numeric(sim_C_wout),
           theory_I_with = theory_I_with,
           sim_I_with = as.numeric(sim_I_with),
           theory_I_wout = theory_I_wout,
           sim_I_wout = as.numeric(sim_I_wout)))
})
tictoc::toc()

# Name by MOIs
colnames(max_probs) <- sapply(all_MOIs, function(x) paste(x, collapse = ""))

# Save as exported data
usethis::use_data(max_probs, overwrite = TRUE)
