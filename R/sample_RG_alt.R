#' Sample a transitive relationship graph, alternate version
#'
#' Use the techniques in enumerate_RGs_alt instead
#' @export
sample_RG_alt <- function(MOIs){
  # Check MOIs are positive whole numbers
  if (!all(is.wholenumber(MOIs)) | any(MOIs < 1)) stop("MOIs should be positive integers")

  infection_count <- length(MOIs) # Number of time points
  gs_count <- sum(MOIs) # Number of genotypes

  # remove restriction for now as this method should scale better
  # if (gs_count > 6 | infection_count > 3) stop("Sorry, too many RGs")
  # if (gs_count <= 1) stop("Sorry, no RGs")

  gs <- paste0("g", 1:gs_count) # Genotype names
  ts_per_gs <- rep(1:infection_count, MOIs)

  CP_list <- enumerate_CPs(MOIs)
  sizes <- sapply(CP_list, function(x) ncol(part.list[[max(x)]]))
  CP <- unlist(sample(CP_list, 1, prob=sizes))

  n.units <- max(CP)
  units <- split(1:gs_count, CP)

  # given clonal relationships, generate all compatible sibling relationships
  if (n.units > 10) stop("Currently supports up to 10 genotypes")

  sib.parts <- part.list[[n.units]]
  j <- sample(1:ncol(sib.parts), 1)

  adj_all <- array(0, dim = rep(gs_count, 2), dimnames = list(gs,gs))

  n.sib.units <- max(sib.parts[,j])
  sib.units <- split(1:n.units, sib.parts[,j])
  for(k in 1:n.sib.units) {
    selection <- unlist(units[sib.units[[k]]])
    adj_all[selection, selection] <- 0.5;
  }
  for(u in 1:n.units) {
    adj_all[units[[u]], units[[u]]] <- 1;
  }
  diag(adj_all) <- 0

  # Convert adjacency matrix into an igraph item
  RG <- igraph::graph_from_adjacency_matrix(adjmatrix = adj_all,
                                            mode = "lower",
                                            diag = F,
                                            weighted = T)
  # Add a time-point attribute
  RG <- igraph::set_vertex_attr(RG, "group", value = ts_per_gs)

  return(RG)
}
