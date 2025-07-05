#' Plot a relationship graph (RG)
#'
#' This function is a wrapper around \code{\link[igraph]{plot.igraph}}, written
#' to group parasite genotypes by infection, both spatially and using vertex
#' colour. Specifically, parasite genotypes within infections are vertically
#' distributed with some horizontal jitter when \code{layout_by_group} is TRUE
#' (default), and coloured the same. It also makes sure clonal and sibling edges
#' are plotted differently using different line types.
#'
#' @param RG A relationship graph, which is an \code{igraph} graph; see return
#'   value of \code{\link{RG_to_igraph}}.
#' @param edge_col A vector of edge colours corresponding to different
#'   relationships, where 0.5 represents a sibling and 1 represents a clone.
#' @param edge_lty A vector of edge line types corresponding to different
#'   relationships, where 0.5 represents a sibling and 1 represents a clone.
#' @param vertex_palette A character string specifying an RColorBrewer palette.
#'   Overrides the default \code{palette} of \code{\link[igraph]{plot.igraph}}.
#' @param layout_by_group A logical argument which if TRUE (default) overrides
#'   the default layout of \code{\link[igraph]{plot.igraph}} so that vertices
#'   that represent parasite genotypes from different infections are distributed
#'   horizontally and vertices that represent genotypes within infections are
#'   distributed vertically.
#' @param edge_width Overrides the default \code{edge.width} of
#'   \code{\link[igraph]{plot.igraph}}.
#' @param ... Additional arguments to pass to \code{\link[igraph]{plot.igraph}}, e.g.
#'   \code{edge.curved}.
#'
#'
#' @section Provenance: This function was adapted from \code{plot_Vivax_model} at
#' \url{https://github.com/jwatowatson/RecurrentVivax/blob/master/Genetic_Model/iGraph_functions.R}.
#'
#' @examples
#' RGs <- enumerate_RGs(c(2, 1, 1))
#' cpar <- par() # record current par before changing
#' par(mfrow = c(3, 4), mar = c(0.1, 0.1, 0.1, 0.1))
#' for (i in 12:23) {
#'   plot_RG(RGs[[i]], edge.curved = 0.1)
#'   box()
#' }
#' par(cpar) # reset par
#' @export

plot_RG <- function(RG,
                    layout_by_group = TRUE,
                    vertex_palette = "Set2",
                    edge_lty = c("0.5" = "dashed", "1" = "solid"),
                    edge_col = c("0.5" = "black", "1" = "black"),
                    edge_width = 1.5,
                    ...) {

  ts_per_gs <- igraph::vertex_attr(RG)$group
  if (is.null(ts_per_gs)) stop("RG vertices need a group attribute")
  MOIs <- as.vector(table(ts_per_gs)) # extract MOIs
  infection_count <- length(MOIs)

  if (layout_by_group) { # Compute graph layout, grouping genotypes per episode
    gs_count <- sum(MOIs)
    RG_layout <- array(dim = c(gs_count, 2))
    X <- rbind(seq(0, 1, length.out = infection_count), MOIs)
    z <- as.numeric(rep(MOIs, MOIs) > 2)
    my_offset <- (-1)^(1:length(z)) * (.2 / infection_count) # offset within time point
    RG_layout[, 1] <- unlist(apply(X, 2, function(x) rep(x[1], x[2]))) + my_offset * z
    RG_layout[, 2] <- unlist(sapply(MOIs, function(x) {
      if (x == 1) {
        return(0.5)
      } else {
        return(seq(0, 1, length.out = x))
      }
    }))
  }


  # Create infection colours for vertices
  infection_colours <- RColorBrewer::brewer.pal(n = 8, vertex_palette)
  if (length(infection_colours) < infection_count) {
    infection_colour_fun <- grDevices::colorRampPalette(infection_colours)
    infection_colours <- infection_colour_fun(infection_count)
  }

  # Plot the graph
  igraph::plot.igraph(RG,
    layout = if (layout_by_group) RG_layout else igraph::layout_nicely,
    vertex.color = infection_colours[igraph::vertex_attr(RG)$group],
    edge.color = edge_col[as.character(igraph::edge_attr(RG)$weight)],
    edge.lty = edge_lty[as.character(igraph::edge_attr(RG)$weight)],
    edge.width = edge_width,
    ...
  )
}
