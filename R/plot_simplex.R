#' Plots a 2D simplex
#'
#' Plots a 2D simplex (a triangle with unit sides centered at the origin) onto
#' which per-recurrence posterior probabilities of recrudescence, relapse,
#' reinfection (or any other probability triplet summing to
#' one) can be projected; see [project2D()] and **Examples** below.
#'
#' @param v.labels Vertex labels anticlockwise from top (default:
#'   "Recrudescence", "Relapse", "Reinfection"). If NULL, vertices are not
#'   labelled.
#'
#' @param v.cutoff Number between 0.5 and 1 that separates lower vs higher
#'   probability regions. Use with caution for recrudescence and reinfection
#'   classification; see
#'   ["Understand posterior estimates"](https://aimeertaylor.github.io/Pv3Rs/articles/understand-posterior.html).
#'
#' @param v.colours Vertex colours anticlockwise from top.
#'
#' @param plot.tri Logical; draws the triangular boundary if `TRUE` (default).
#'
#' @param p.coords Matrix of 3D simplex coordinates (e.g., per-recurrence
#'   probabilities of recrudescence, relapse and reinfection), one vector of 3D
#'   coordinates per row, each row is projected onto 2D coordinates using
#'   [project2D()] and plotted as a single simplex point using
#'   [graphics::points()]. If the user provides a vector encoding a probability
#'   triplet summing to one, it is converted to a matrix with one row.
#'
#' @param p.labels Labels of points in \code{p.coords} (default row names of
#'   \code{p.coords}) No labels if `NA`.
#'
#' @param p.labels.pos Position of \code{p.labels}: \code{1} = below, \code{2} =
#'   left, \code{3} = above (default) and \code{4} = right. Can be a single
#'   value or a vector.
#'
#' @param ... Additional parameters passed to [graphics::points()].
#'
#' @examples
#' # Plot 2D simplex
#' plot_simplex(p.coords = diag(3),
#'              p.labels = c("(1,0,0)", "(0,1,0)", "(0,0,1)"),
#'              p.labels.pos = c(1,3,3))
#'
#' # ==============================================================================
#' # Given data on an enrollment episode and a recurrence,
#' # compute the posterior probabilities of the 3Rs and plot the deviation of the
#' # posterior from the prior
#' # ==============================================================================
#'
#' # Some data:
#' y <- list(list(m1 = c('a', 'b'), m2 = c('c', 'd')), # Enrollment episode
#'           list(m1 = c('a'), m2 = c('c'))) # Recurrent episode
#'
#' # Some allele frequencies:
#' fs <- list(m1 = setNames(c(0.4, 0.6), c('a', 'b')),
#'            m2 = setNames(c(0.2, 0.8), c('c', 'd')))
#'
#' # A vector of prior probabilities:
#' prior <- array(c(0.2, 0.3, 0.5), dim = c(1,3),
#'                dimnames = list(NULL, c("C", "L", "I")))
#'
#' # Compute posterior probabilities
#' post <- compute_posterior(y, fs, prior, progress.bar = FALSE)
#'
#' # Plot simplex with the prior and posterior
#' plot_simplex(p.coords = rbind(prior, post$marg),
#'              p.labels = c("Prior", "Posterior"),
#'              pch = 20)
#'
#' # Add the deviation between the prior and posterior: requires obtaining 2D
#' # coordinates manually
#' xy_prior <- project2D(as.vector(prior))
#' xy_post <- project2D(as.vector(post$marg))
#' arrows(x0 = xy_prior["x"], x1 = xy_post["x"],
#'        y0 = xy_prior["y"], y1 = xy_post["y"], length = 0.1)
#' @export
plot_simplex <- function(v.labels = c("Recrudescence", "Relapse", "Reinfection"),
                         v.cutoff = 0.5,
                         v.colours = c("yellow","purple","red"),
                         plot.tri = TRUE,
                         p.coords = NULL,
                         p.labels = rownames(p.coords),
                         p.labels.pos = 3,
                         ...) {

  # Define some constants:
  h <- sqrt(3)/2 # Height of equilateral triangle with unit sides
  r <- 1/sqrt(3) # Radius of circle encompassing triangle with unit sides
  k <- h-r # Distance from (0,0) to bottom of triangle with unit sides

  # Null plot
  plot(NULL, xlim = c(-0.6, 0.6), ylim = c(-(k + 0.1), r + 0.1), asp = 1,
       xaxt = "n", yaxt = "n", bty = "n",
       ylab = "", xlab = "")

  # Plot equilateral triangle:
  if(plot.tri) graphics::polygon(x = c(-0.5, 0.5, 0), y = c(-k, -k, r))

  # Annotate vertices:
  if (!is.null(v.labels)) {
    graphics::text(x = c(0, -0.5, 0.5),
                   y = c(r, -k, -k),
                   labels = v.labels,
                   pos = c(3,1,1))
  }

  # Project critical points to colour vertices (L for left, T for top, R for
  # right)
  xyT <- project2D(c(1,0,0))
  xyL <- project2D(c(0,1,0))
  xyR <- project2D(c(0,0,1))

  xyTLR <- project2D(c(1/3,1/3,1/3))
  xyTL <- project2D(c(0.5,0.5,0))
  xyTR <- project2D(c(0.5,0,0.5))
  xyLR <- project2D(c(0,0.5,0.5))

  # Delineate weak classification
  p1 <- rbind(xyT, xyTR, xyTLR, xyTL)
  p2 <- rbind(xyL, xyTL, xyTLR, xyLR)
  p3 <- rbind(xyR, xyLR, xyTLR, xyTR)

  graphics::polygon(x = p1[,"x"], y = p1[,"y"], border = NA, col = grDevices::adjustcolor(v.colours[1], alpha.f = 0.35))
  graphics::polygon(x = p2[,"x"], y = p2[,"y"], border = NA, col = grDevices::adjustcolor(v.colours[2], alpha.f = 0.35))
  graphics::polygon(x = p3[,"x"], y = p3[,"y"], border = NA, col = grDevices::adjustcolor(v.colours[3], alpha.f = 0.35))


  # Delineate strong classification
  if(!is.null(v.cutoff)) {

    if(v.cutoff < 0.5 | v.cutoff > 1) stop("The v.cutoff should be a number between 0.5 and 1")

    # Project critical points to colour high probability regions (L for left, T
    # for top, R for right)
    xyTL_T <- project2D(c(v.cutoff,1-v.cutoff,0))
    xyTL_L <- project2D(c(1-v.cutoff,v.cutoff,0))
    xyTR_T <- project2D(c(v.cutoff,0,1-v.cutoff))
    xyTR_R <- project2D(c(1-v.cutoff,0,v.cutoff))
    xyLR_L <- project2D(c(0,v.cutoff,1-v.cutoff))
    xyLR_R <- project2D(c(0,1-v.cutoff,v.cutoff))

    # Regions with more than v.cutoff probability
    p1 <- rbind(xyT, xyTR_T, xyTL_T)
    p2 <- rbind(xyL, xyTL_L, xyLR_L)
    p3 <- rbind(xyR, xyLR_R, xyTR_R)

    # Triangles plotted on top of quadrilaterals (0.35 + 0.45 = 0.8)
    graphics::polygon(x = p1[,"x"], y = p1[,"y"], border = NA, col = grDevices::adjustcolor(v.colours[1], alpha.f = 0.45))
    graphics::polygon(x = p2[,"x"], y = p2[,"y"], border = NA, col = grDevices::adjustcolor(v.colours[2], alpha.f = 0.45))
    graphics::polygon(x = p3[,"x"], y = p3[,"y"], border = NA, col = grDevices::adjustcolor(v.colours[3], alpha.f = 0.45))
  }

  # Plot points if given
  if(is.null(p.coords)) return()
  # if p.coords is a single vector, make it a matrix of one row
  if(is.vector(p.coords)) p.coords <- t(p.coords)
  # note that p.coords has one point per row, p_2Dcoords has one point per column
  xy <- apply(p.coords, 1, project2D)
  graphics::points(x = xy["x",], y = xy["y",], ...)
  graphics::text(x = xy["x",], y = xy["y",], pos = p.labels.pos, labels = p.labels)
}


#' Project 3D probability coordinates onto 2D simplex coordinates
#'
#' Project three probabilities that sum to one (e.g., per-recurrence
#' probabilities of recrudescence, relapse and reinfection) onto the coordinates
#' of a 2D simplex centred at the origin (i.e., a triangle centred at (0,0) with
#' unit-length sides).
#'
#' The top, left, and right vertices of the 2D simplex correspond with the
#' first, second and third entries of `v`, respectively. Each probability is
#' proportional to the distance from the point on the simplex to the side
#' opposite the corresponding probability; see **Examples** below and
#' [plot_simplex()] for more details.
#'
#' @param v A numeric vector of three numbers in zero to one that sum to one.
#'
#' @return A numeric vector of two coordinates that can be used to plot the
#'   probability vector `v` on the origin-centred 2D simplex.
#'
#' @examples
#' probabilities_of_v1_v2_v3 <- c(0.75,0.20,0.05)
#' coordinates <- project2D(v = probabilities_of_v1_v2_v3)
#'
#' # Plot probability vector on 2D simplex
#' plot_simplex(v.labels = c("v1", "v2", "v3"))
#' points(x = coordinates[1], y = coordinates[2], pch = 20)
#'
#' # Plot the distances that represent probabilities
#' # get vertices, get points on edges by orthogonal projection, plot arrows
#' v <- apply(matrix(c(1,0,0,0,1,0,0,0,1), nrow = 3), 1, project2D)
#' p3 <- v[,1] + sum((coordinates - v[,1]) * (v[,2] - v[,1])) * (v[,2] - v[,1])
#' p1 <- v[,2] + sum((coordinates - v[,2]) * (v[,3] - v[,2])) * (v[,3] - v[,2])
#' p2 <- v[,3] + sum((coordinates - v[,3]) * (v[,1] - v[,3])) * (v[,1] - v[,3])
#' arrows(x0 = coordinates[1], y0 = coordinates[2], x1 = p1[1], y1 = p1[2], length = 0.1)
#' arrows(x0 = coordinates[1], y0 = coordinates[2], x1 = p2[1], y1 = p2[2], length = 0.1)
#' arrows(x0 = coordinates[1], y0 = coordinates[2], x1 = p3[1], y1 = p3[2], length = 0.1)
#'
#' @export
project2D <- function(v){

  # Define some constants:
  h <- sqrt(3)/2 # Height of equilateral triangle with unit sides
  r <- 1/sqrt(3) # Radius of circle encompassing triangle with unit sides

  if(!abs(1-round(sum(v))) < .Machine$double.eps^0.5) stop('Input does not sum to one')
  A <- matrix(c(0.5,0,1,h,0,0), byrow = T, nrow = 2) # Projection matrix
  xy <- as.vector(A%*%v) - c(0.5, h-r) # Project and translate
  names(xy) <- c('x','y')

  return(xy)
}







