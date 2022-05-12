#' Sample a transitive relationship graph (RG) for a single recurrence
#'
#' Samples one transitive graph of stranger and sibling relationships between
#' distinct parasite genotypes within the initial and recurrent infections and,
#' depending on the specified cause, stranger, sibling and/or clonal
#' relationships across parasite genotypes in the initial and recurrent
#' infections. The function requires the user to specify the multiplicity of the
#' initial infection, \code{MOI_init}. If \code{MOI_init} is greater than 4,
#' \code{sample_recurrent_RG} is liable to take a long time because it is
#' testing for transitivity while sampling among not necessarily transitive
#' graphs, of which there will be many; see the message that
#' \code{sample_recurrent RG} prints to the screen.
#'
#' @param MOI_init The number of distinct parasite genotypes (a.k.a the
#'   multiplicity of infection, MOI) of the initial infection.
#' @param cause A character string specifying the cause of the recurrent
#'   infection: \code{"recrudescence"}, \code{"reinfection"} or
#'   \code{"relapse"}.
#' @param prob_recrudesce A probability that each genotype survives in a
#'   recrudescent infection.
#' @param beta_relapse A multiplicative factor that modifies the average MOI of
#'   a reinfection to obtain the average MOI of a relapse.
#'
#' @return Returns a transitive relationship graph (RG); see \code{\link{enumerate_RGs}} for details.
#'
#' @section To-do:
#' Better to test transitivity on intra before whole?
#'
#' @examples
#'
#'for (cause in c("recrudescence", "reinfection", "relapse")) {
#'   RG <- sample_recurrent_RG(MOI_init = 3, cause)
#'   plot_RG(RG)
#'   print(RG)
#' }
#'
#'
#' @export
sample_recurrent_RG <- function(MOI_init, cause,
                               prob_recrudesce = 0.5, beta_relapse = 1){

  # Check MOI is a whole number
  if (!is.wholenumber(MOI_init)) stop("MOI should be a positive integer")
  cause %in% c("recrudescence", "reinfection", "relapse")

  # Sample the MOI of the recurrent infection
  MOI_recur <- 0
  if (cause == "reinfection") {
    while(MOI_recur < 1) MOI_recur <- rpois(n = 1, lambda = MOI_init)
  } else if (cause == "relapse") {
    while(MOI_recur < 1) MOI_recur <- rpois(n = 1, lambda = MOI_init * beta_relapse)
  } else if (cause == "recrudescence") {
    while(MOI_recur < 1) MOI_recur <- rbinom(n = 1, size = MOI_init, prob = prob_recrudesce)
  }

  MOIs <- c(MOI_init, MOI_recur)
  intra_edge_counts <- sapply(MOIs, choose, k = 2)
  inter_edge_count <- prod(MOIs)
  gs_count <- sum(MOIs) # Number of genotypes
  gs <- paste0("g", 1:gs_count) # Genotype names
  ts_per_gs <- rep(1:2, MOIs) # Infection indices per genotype
  gs_per_ts <- split(gs, ts_per_gs) # Genotypes per infection

  # Hard code relationship types to satisfy the test_transitive function
  relationship_types <- c(stranger = 0, sibling = 0.5, clone = 1)
  intra_relationship_types <- relationship_types[setdiff(names(relationship_types), "clone")]

  # Compute the size of the not necessarily transitive sample space
  if (cause == "relapse") {
    RGs_to_eval_count <- prod(3^inter_edge_count, 2^intra_edge_counts)
    if(RGs_to_eval_count > 1000){
      writeLines(sprintf("\nThis might take a while: sampling among %s not necessarily transitive RGs",
                         RGs_to_eval_count))
    }
  } else {
    RGs_to_eval_count <- prod(2^intra_edge_counts)
    if(RGs_to_eval_count > 1000){
      writeLines(sprintf("\nThis might take a while: sampling among %s not necessarily transitive RGs",
                         RGs_to_eval_count))
    }
  }

  # Allocate adjacency matrices
  adj_all <- array(NA, dim = rep(gs_count, 2), dimnames = list(gs,gs))
  intra_adjs <- lapply(MOIs, function(MOI) array(NA, dim = rep(MOI, 2)))

  transitive <- FALSE
  while (!transitive) {

    # Sample the within infection relationships and populate adj_all
    if (cause == "reinfection" | cause == "relapse") {
      for(t in 1:2) {
        intra_relationships <- sample(intra_relationship_types, choose(MOIs[t],2), replace = T)
        intra_adjs[[t]][lower.tri(intra_adjs[[t]])] <- intra_relationships
        adj_all[gs_per_ts[[t]], gs_per_ts[[t]]] <- intra_adjs[[t]]
      }
    } else if (cause == "recrudescence") {
      relationships_init <- sample(intra_relationship_types, choose(MOI_init,2), replace = T)
      intra_adjs[[1]][lower.tri(intra_adjs[[1]])] <- relationships_init
      intra_adjs[[2]] <- intra_adjs[[1]][1:MOI_recur, 1:MOI_recur] # Copy first MOI_recur
      for(t in 1:2) adj_all[gs_per_ts[[t]], gs_per_ts[[t]]] <- intra_adjs[[t]] # Populate adj_all
    }

    # Specify/sample and populate the between time-point relationship
    if (cause == "reinfection") {
      inter_relationships <- rep(relationship_types["stranger"], inter_edge_count)
      adj_all[gs_per_ts[[2]],gs_per_ts[[1]]] <- inter_relationships
    } else if (cause == "relapse") {
      inter_relationships <- sample(relationship_types, inter_edge_count, replace = T)
      adj_all[gs_per_ts[[2]],gs_per_ts[[1]]] <- inter_relationships
    } else if (cause == "recrudescence") {
      # Extract intra recrudescent
      intra_recur <- intra_adjs[[1]][1:MOI_recur, 1:MOI_recur]
      # Make symmetrical
      intra_recur[upper.tri(intra_recur)] <- 0; intra_recur <- t(intra_recur) + intra_recur
      # Add clonal relationships
      diag(intra_recur) <- 1
      # Give recrudescent new names
      rownames(intra_recur) <- gs_per_ts[[2]]

      # Extract inter not recrudescent and recrudescent
      if(MOI_init > MOI_recur) {
        inter <- intra_adjs[[1]][(MOI_recur+1):MOI_init, 1:MOI_recur, drop = F]
        rownames(inter) <- gs_per_ts[[1]][(MOI_recur+1):MOI_init]
        colnames(inter) <- gs_per_ts[[2]]
        adj_all[colnames(inter), rownames(inter)] <- inter
      }

      # Populate adj_all
      adj_all[gs_per_ts[[2]],gs_per_ts[[1]][1:MOI_recur]] <- intra_recur
    }

    # Convert adjacency matrix into an igraph item
    RG <- igraph::graph_from_adjacency_matrix(adjmatrix = adj_all,
                                              mode = "lower",
                                              diag = F,
                                              weighted = T)
    # Add a time-point attribute
    RG <- igraph::set_vertex_attr(RG, "group", value = ts_per_gs)

    # Test transitivity and store if transitive
    transitive <- test_transitive(RG)
  }

  return(RG)
}



