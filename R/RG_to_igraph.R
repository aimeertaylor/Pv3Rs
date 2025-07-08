#' Converts a relationship graph (RG) encoded as a list to an \code{igraph} object
#'
#' Converts an RG encoded as a list to an \code{igraph} object, which requires
#' more memory allocation but can be plotted using \code{\link{plot_RG}}.
#'
#' @param RG List encoding an RG; see **Value** of
#'   \code{\link{enumerate_RGs}} when \code{igraph = FALSE}.
#' @param MOIs Vector of per-episode multiplicities of infection (MOIs), i.e.,
#'   numbers of vertices per episode; adds to the graph an attribute
#'   that is used by \code{\link{plot_RG}} to group genotypes by episodes.
#'
#' @return A weighted graph whose edge weights 1 and 0.5 encode clonal and
#'   sibling relationships, respectively.
#'
#' @examples
#' MOIs <- c(3,2)
#' set.seed(6)
#' RG_as_list <- sample_RG(MOIs, igraph = FALSE)
#' RG_as_igraph <- RG_to_igraph(RG_as_list,  MOIs)
#'
#' # RG encoded as a list requires less memory allocation
#' utils::object.size(RG_as_list)
#' utils::object.size(RG_as_igraph)
#'
#' # RG encoded as an igraph object can be plotted using plot_RG() and
#' # manipulated using igraph functions
#' plot_RG(RG_as_igraph, margin = rep(0,4), vertex.label = NA)
#'
#' # Edge weights 1 and 0.5 encode clonal and sibling relationships
#' igraph::E(RG_as_igraph)$weight
#'
#' # Vertex attribute group encodes episode membership
#' igraph::V(RG_as_igraph)$group
#'
#' @export
RG_to_igraph <- function(RG, MOIs) {

  ts_per_gs <- rep(1:length(MOIs), MOIs)
  gs_count <- sum(MOIs)
  gs <- unlist(RG$clone)

  if (length(gs) != gs_count) stop("sum(MOIs) should equal number of genotypes in RG")

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
