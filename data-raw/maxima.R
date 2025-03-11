################################################################################
# Script to generate maximum probabilities for all MOI vectors summing to
# at-most eight assuming recurrent states are equally likely a priori. Takes 5
# hours to run.
################################################################################
rm(list = ls())
set.seed(1)
n <- 10000

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

tictoc::tic()
maxima <- sapply(all_MOIs[1:10], function(MOIs){

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
  states <- unique(unlist(CIL_gvn_RGs))
  states_C_1st <- states[do.call(rbind, strsplit(states, split = ""))[,1] == "C"]
  states_I_1st <- states[do.call(rbind, strsplit(states, split = ""))[,1] == "I"]

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
  RGs_with <- 1:length(RGs)
  RGs_wout <- which(keep_log)
  RGs_C_with <- which(RGs_C_log)
  RGs_I_with <- which(RGs_I_log)
  RGs_C_wout <- which(RGs_C_log & keep_log)
  RGs_I_wout <- which(RGs_I_log & keep_log)

  RGs_gvn_CIL <- sapply(states, function(s) which(sapply(CIL_gvn_RGs, function(x) s%in%x)), simplify = F)
  RGsC_gvn_CIL <- sapply(RGs_gvn_CIL, function(x) intersect(x, RGs_C_with), simplify = F)
  RGsI_gvn_CIL <- sapply(RGs_gvn_CIL, function(x) intersect(x, RGs_I_with), simplify = F)

  RGs_gvn_CIL_wout <- sapply(RGs_gvn_CIL, function(x) intersect(x, RGs_wout), simplify = F)
  RGsC_gvn_CIL_wout <- sapply(RGs_gvn_CIL, function(x) intersect(x, RGs_C_with), simplify = F)
  RGsI_gvn_CIL_wout <- sapply(RGs_gvn_CIL, function(x) intersect(x, RGs_I_with), simplify = F)


  #=============================================================================
  # Compute theoretical maximum probabilities with and w/out summation over all
  #=============================================================================
  if(sum(RGs_C_log) == 0){
    theory_C_with <- 0
    theory_C_wout <- 0
  } else {
    theory_s_with_un <- sapply(RGsC_gvn_CIL, length)/
      sapply(RGs_gvn_CIL, length)
    theory_s_wout_un <- sapply(RGsC_gvn_CIL_wout, length)/
      sapply(RGs_gvn_CIL_wout, length)
    theory_C_with <- sum((theory_s_with_un/sum(theory_s_with_un))[states_C_1st])
    theory_C_wout <- sum((theory_s_wout_un/sum(theory_s_wout_un))[states_C_1st])
  }

  theory_s_with_un <- sapply(RGsI_gvn_CIL, length)/
    sapply(RGs_gvn_CIL, length)
  theory_s_wout_un <- sapply(RGsI_gvn_CIL_wout, length)/
    sapply(RGs_gvn_CIL_wout, length)
  theory_I_with <- sum((theory_s_with_un/sum(theory_s_with_un))[states_I_1st])
  theory_I_wout <- sum((theory_s_wout_un/sum(theory_s_wout_un))[states_I_1st])

  #=============================================================================
  # Compute numerical maximum probabilities with and w/out summation over all
  #=============================================================================
  ss <- sample(states, n, replace = T) # sample recurrence states uniformly at random

  # as.numeric(as.character()) needed to override convenience behaviour; see ?sample
  graphs_with <- sapply(ss, function(s) {
    as.numeric(sample(as.character(RGs_gvn_CIL[[s]]), 1))})
  graphs_wout <- sapply(ss, function(s) {
    as.numeric(sample(as.character(intersect(RGs_gvn_CIL[[s]], RGs_wout)), 1))})

  if(sum(RGs_C_log) == 0) {
    sim_C_with <- 0
    sim_C_wout <- 0
  } else {
    # Of graphs compatible with C 1st, proportion generated by sequences with C 1st
    ind <- graphs_with %in% RGs_C_with
    sim_C_with <- sum((table(names(graphs_with[ind]))/sum(ind))[states_C_1st])

    ind <- graphs_wout %in% RGs_C_wout # Among graphs compatible with C
    sim_C_wout <- sum((table(names(graphs_wout[ind]))/sum(ind))[states_C_1st])
  }

  # Of graphs compatible with I 1st, proportion generated by sequences with I 1st
  ind <- graphs_with %in% RGs_I_with
  sim_I_with <- sum((table(names(graphs_with[ind]))/sum(ind))[states_I_1st])

  ind <- graphs_wout %in% RGs_I_wout
  sim_I_wout <- sum((table(names(graphs_wout[ind]))/sum(ind))[states_I_1st])

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
colnames(maxima) <- sapply(all_MOIs, function(x) paste(x, collapse = ""))

# Save as exported data
#usethis::use_data(maxima, overwrite = TRUE)

