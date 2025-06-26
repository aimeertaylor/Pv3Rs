################################################################################
# Script to generate maximum probabilities for MOI vectors summing to at-most
# eight assuming recurrence states are equally likely a priori. Run time ~ 20 min
################################################################################
rm(list = ls())

# Generate all MOI combinations given a single recurrence
y <- gtools::permutations(n = 7, r = 2, repeats.allowed = T)
z <- y[rowSums(y) < 9, ] # Total genotype count of 8
all_MOIs <- c(lapply(1:nrow(z), function(i) z[i,]), # Make into a list
              sapply(3:8, function(x) rep(1,x), simplify = F)) # All monoclonal given two or more recurrences

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
maxima <- sapply(all_MOIs, function(MOIs){

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
      do.call(rbind, strsplit(x, split = ""))[,1] # Get state of 1st recurrence
    })
    RGs_C_log <- sapply(First_CIL_gvn_RGs, function(x) "C" %in% x) # C compatible
    RGs_I_log <- sapply(First_CIL_gvn_RGs, function(x) "I" %in% x) # I compatible
  }

  # Get the size of the largest intra-episode sibling clique per graph
  max_clique_sizes <- sapply(RGs, get_max_clique_size, intra_edges = intra_edges)
  keep_log <- max_clique_sizes < 3 # Keep only those with at most sibling pairs
  RGs_with <- 1:length(RGs) # Count all RGs
  RGs_wout <- which(keep_log) # Count all RGs with at most sibling pairs
  RGs_C_with <- which(RGs_C_log) # Repeat for RGs compatible with C 1st
  RGs_I_with <- which(RGs_I_log)
  RGs_C_wout <- which(RGs_C_log & keep_log) # Repeat for RGs compatible with I 1st
  RGs_I_wout <- which(RGs_I_log & keep_log)

  # Count RGs for different states
  RGs_gvn_CIL <- sapply(states, function(s) which(sapply(CIL_gvn_RGs, function(x) s%in%x)), simplify = F)
  RGsC_gvn_CIL <- sapply(RGs_gvn_CIL, function(x) intersect(x, RGs_C_with), simplify = F)
  RGsI_gvn_CIL <- sapply(RGs_gvn_CIL, function(x) intersect(x, RGs_I_with), simplify = F)

  RGs_gvn_CIL_wout <- sapply(RGs_gvn_CIL, function(x) intersect(x, RGs_wout), simplify = F)
  RGsC_gvn_CIL_wout <- sapply(RGs_gvn_CIL, function(x) intersect(x, RGs_C_with), simplify = F)
  RGsI_gvn_CIL_wout <- sapply(RGs_gvn_CIL, function(x) intersect(x, RGs_I_with), simplify = F)

  # Compute theoretical maximum probabilities with and w/out summation over all
  if(sum(RGs_C_log) == 0){
    C_with <- 0
    C_wout <- 0
  } else {

    # Compute vectors of un-normalised graph space ratios
    s_with_un <- sapply(RGsC_gvn_CIL, length)/sapply(RGs_gvn_CIL, length)
    s_wout_un <- sapply(RGsC_gvn_CIL_wout, length)/sapply(RGs_gvn_CIL_wout, length)

    # Example of above for C, L, I given single recurrence:
    # (|Gc|/|Gc|, |Gc|/|GL|, 0/|GI|) = (1, |Gc|/|GL|, 0)

    # Example of above for CC, IC, LC... LL given two recurrences:
    # (|Gcc|/|Gcc|, 0/|GIc|, |Gcc|/|GLC|, ..., |GcL|/|GLL|) = (1, 0, |Gcc|/|GLC|, ..., |GcL|/|GLL|)

    # Sum over vectors of normalised graph space ratios (taking marginal)
    C_with <- sum((s_with_un/sum(s_with_un))[states_C_1st])
    C_wout <- sum((s_wout_un/sum(s_wout_un))[states_C_1st])

    # Example of above for C, L, I given single recurrence:
    # (1, |Gc|/|GL|, 0) / (1 + (|Gc|/|GL|)) = (|GL|/(|GC + GL|), |GC|/(|GC + GL|), 0)
  }

  # Compute vectors of un-normalised graph space ratios
  s_with_un <- sapply(RGsI_gvn_CIL, length)/sapply(RGs_gvn_CIL, length)
  s_wout_un <- sapply(RGsI_gvn_CIL_wout, length)/sapply(RGs_gvn_CIL_wout, length)

  # Sum over vectors of normalised graph space ratios (taking marginal)
  I_with <- sum((s_with_un/sum(s_with_un))[states_I_1st])
  I_wout <- sum((s_wout_un/sum(s_wout_un))[states_I_1st])


  return(c(C_with = C_with,
           C_wout = C_wout,
           I_with = I_with,
           I_wout = I_wout))
})
tictoc::toc()

# Name by MOIs
colnames(maxima) <- sapply(all_MOIs, function(x) paste(x, collapse = ""))

# Save as exported data
usethis::use_data(maxima, overwrite = TRUE)

