#' Converts a relationship graph (RG) encoded as a list to an \code{igraph} object
#'
#' Converts an RG encoded as a list to an \code{igraph} object, which requires
#' more memory allocation but can be plotted using \code{\link{Plot_RG}}.
#'
#' @param RG List encoding an RG; see **Value** of
#'   \code{\link{enumerate_RGs}} when \code{igraph = FALSE}.
#' @param gs Character vector of genotype names.
#' @param ts_per_gs Numeric vector containing episode number for each genotype,
#'   e.g., `rep(1:length(y), determine_MOIs(y))`.
#'
#' @return A weighted graph whose edge weights 1 and 0.5 encode clonal and
#'   sibling relationships, respectively.
#'
#' @examples
#' set.seed(5)
#' RG_as_list <- sample_RG(c(3, 2), igraph = FALSE)
#' utils::object.size(RG_as_list)
#' RG_as_igraph <- RG_to_igraph(RG_as_list, paste0("g", 1:5), c(1, 1, 1, 2, 2))
#' utils::object.size(RG_as_igraph)
#' plot_RG(RG_as_igraph)
#' igraph::E(RG_as_igraph)$weight
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
