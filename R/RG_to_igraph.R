#' Converts relationship graph to an \code{igraph} object
#'
#' Makes a relationship graph object output by \code{\link{enumerate_RGs_alt}}
#' compatible with the \code{igraph} output of \code{\link{enumerate_RGs}}.
#'
#' @param RG Relationship graph output by \code{\link{enumerate_RGs_alt}} with
#'   \code{igraph=FALSE}.
#' @param gs Vector of genotype names.
#' @param ts_per_gs Vector of infection numbers for each genotype.
#'
#' @return An \code{igraph} object along with the original variables in
#'   \code{RG}.
#'
#' @export
RG_to_igraph <- function(RG, gs, ts_per_gs) {
  gs_count <- length(gs)

  n.sib.clones <- max(RG$sib.vec) # number of sibling cells
  n.clones <- max(RG$clone.vec) # number of clonal cells
  clones <- RG$clone # clonal partition as list

  # make adjacent matrix
  adj_all <- array(0, dim = rep(gs_count, 2), dimnames = list(gs, gs))
  for (k in 1:n.sib.clones) {
    selection <- unlist(clones[RG$sib[[k]]])
    adj_all[selection, selection] <- 0.5 # sibling and clonal edges
  }
  for (u in 1:n.clones) {
    adj_all[clones[[u]], clones[[u]]] <- 1 # overwrite clonal edges
  }
  diag(adj_all) <- 0 # remove self loops

  # Convert adjacency matrix into an igraph item
  RG_igraph <- igraph::graph_from_adjacency_matrix(
    adjmatrix = adj_all,
    mode = "lower",
    diag = F,
    weighted = T
  )
  # Add a time-point attribute
  RG_igraph <- igraph::set_vertex_attr(RG_igraph, "group", value = ts_per_gs)

  RG_igraph$clone <- RG$clone
  RG_igraph$clone.vec <- RG$clone.vec
  RG_igraph$sib <- RG$sib
  RG_igraph$sib.vec <- RG$sib.vec

  RG_igraph
}
