#' Plot a relationship graph (RG)
#'
#' This function is a wrapper around \code{\link[igraph]{plot.igraph}}, written
#' to group parasite genotypes by episode both spatially and using vertex
#' colour (specifically, parasite genotypes within episodes are vertically
#' distributed with some horizontal jitter when \code{layout.by.group = TRUE}
#' (default), and equicolored), and to ensure clone and sibling edges
#' are plotted using different line types.
#'
#' To see how to plot relationship graphs outputted by
#' \code{\link{compute_posterior}}, please refer to **Exploration of relationship graphs** in
#' \HTMLVignette{demonstrate-usage}{exploration-of-relationship-graphs}{Demonstrate Pv3Rs usage}.
#'
#' @param RG \code{igraph} object encoding an RG; see
#'   \code{\link{RG_to_igraph}}.
#' @param edge.col Named vector of edge colours corresponding to different
#'   relationships.
#' @param edge.lty Named vector of edge line types corresponding to different
#'   relationships.
#' @param vertex.palette A character string specifying an RColorBrewer palette.
#'   Overrides the default \code{palette} of \code{\link[igraph]{plot.igraph}}.
#' @param layout.by.group Logical; if TRUE (default) overrides
#'   the default layout of \code{\link[igraph]{plot.igraph}} so that vertices
#'   that represent parasite genotypes from different episodes are distributed
#'   horizontally and vertices that represent genotypes within episodes are
#'   distributed vertically.
#' @param edge.width Overrides the default \code{edge.width} of
#'   \code{\link[igraph]{plot.igraph}}.
#' @param ... Additional arguments to pass to \code{\link[igraph]{plot.igraph}}, e.g.,
#'   \code{edge.curved}.
#'
#'
#' @section Provenance: This function was adapted from \code{plot_Vivax_model} at
#' \url{https://github.com/jwatowatson/RecurrentVivax/blob/master/Genetic_Model/iGraph_functions.R}.
#'
#' @return None
#'
#' @examples
#' RGs <- enumerate_RGs(c(2, 1, 1), progress.bar = FALSE)
#' oldpar <- par(no.readonly = TRUE) # record user's options
#' par(mfrow = c(3, 4), mar = c(0.1, 0.1, 0.1, 0.1))
#' for (i in 12:23) {
#'   plot_RG(RGs[[i]],
#'   edge.col = c(sibling = "gray", clone = "black"),
#'   edge.lty = c(sibling = "dotted", clone = "solid"),
#'   edge.curved = 0.1)
#'   box()
#' }
#' par(oldpar) # restore user's options
#' @export

plot_RG <- function(RG,
                    layout.by.group = TRUE,
                    vertex.palette = "Set2",
                    edge.lty = c(sibling = "dashed", clone = "solid"),
                    edge.col = c(sibling = "black", clone = "black"),
                    edge.width = 1.5,
                    ...) {

  # Convert relationship names to weights
  widths_by_relationship <- c(sibling = 0.5, clone = 1)
  names(edge.lty) <- as.character(widths_by_relationship[names(edge.lty)])
  names(edge.col) <- as.character(widths_by_relationship[names(edge.col)])

  ts_per_gs <- igraph::vertex_attr(RG)$group
  if (is.null(ts_per_gs)) stop("RG vertices need a group attribute")
  MOIs <- as.vector(table(ts_per_gs)) # extract MOIs
  infection_count <- length(MOIs)

  if (layout.by.group) { # Compute graph layout, grouping genotypes per episode
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
  infection_colours <- RColorBrewer::brewer.pal(n = 8,vertex.palette)
  if (length(infection_colours) < infection_count) {
    infection_colour_fun <- grDevices::colorRampPalette(infection_colours)
    infection_colours <- infection_colour_fun(infection_count)
  }

  # Plot the graph
  igraph::plot.igraph(RG,
    layout = if (layout.by.group) RG_layout else igraph::layout_nicely,
    vertex.color = infection_colours[igraph::vertex_attr(RG)$group],
    edge.color = edge.col[as.character(igraph::edge_attr(RG)$weight)],
    edge.lty = edge.lty[as.character(igraph::edge_attr(RG)$weight)],
    edge.width = edge.width,
    ...
  )
}
