################################################################################
# Generate distribution over first edge type given 2 to 4 monoclonal episodes
################################################################################

# Number of clonal, sibling and stranger edges given L, LL, and LLL
n_recur <- 2:4
names(n_recur) <- paste0(n_recur, " recurrence")
counts <- sapply(n_recur, function(k){

  MOIs <- rep(1,k+1)

  # Get requisite information for maximum probability computation
  gs <- paste0("g", 1:sum(MOIs)) # genotype names (graph vertices)
  ts <- 1:length(MOIs) # episode indices
  ts_per_gs <- rep(ts, MOIs) # episode index of each genotype
  gs_per_ts <- split(gs, ts_per_gs) # genotypes grouped by episode
  RGs <- enumerate_RGs(MOIs, igraph = F) # Get graphs
  CIL_gvn_RGs <- sapply(RGs, compatible_rstrs, gs_per_ts) # compatible states

  # C / I compatibility with first recur (g1 to g2 edge clone / stranger)
  if (all(ts < 3)) {
    RGs_C_log <- sapply(CIL_gvn_RGs, function(x) "C" %in% x) # g1--g2 C
    RGs_I_log <- sapply(CIL_gvn_RGs, function(x) "I" %in% x) # g1--g2 I
  } else {
    First_CIL_gvn_RGs <- sapply(CIL_gvn_RGs, function(x) {
      do.call(rbind, strsplit(x, split = ""))[,1]
    })
    RGs_C_log <- sapply(First_CIL_gvn_RGs, function(x) "C" %in% x) # g1--g2 C
    RGs_I_log <- sapply(First_CIL_gvn_RGs, function(x) "I" %in% x) # g1--g2 I
  }

  # Edge count
  return(c("Clone" = sum(RGs_C_log),
           "Sibling" = (length(RGs)-sum(RGs_C_log)-sum(RGs_I_log)),
           "Stranger" = (sum(RGs_I_log))))

})

rbind(counts, total = colSums(counts)) # As fraction
t(t(counts)/colSums(counts)) # As probability
