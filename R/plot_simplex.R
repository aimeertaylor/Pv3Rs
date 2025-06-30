#' Plots a 2D simplex
#'
#' Plots a 2D simplex, a triangle with unit sides centered at the origin, onto
#' which marginal posterior probabilities of relapse, reinfection and
#' recrudescence (or any other vector of three numbers in zero to one summing to
#' one) can be projected; see [project2D()] and examples below.
#'
#' @param v.labels A vector of labels that annotate vertices anticlockwise from
#'   top (default: "Recrudescence", "Relapse", "Reinfection"). If NULL, vertices
#'   are not annotated.
#'
#' @param v.cutoff An arbitrary number between 0.5 and 1 that separates regions
#'   of lower and higher probability. Beware the use of cut-offs for probable
#'   recrudescence classification and probable reinfection classification; see
#'   ["Understand posterior estimates"](https://aimeertaylor.github.io/Pv3Rs/articles/understand-posterior.html).
#'
#' @param v.colours A vector of colours associated with the vertices
#'   anticlockwise from top; see example below.
#'
#' @param plot.tri Whether to plot the triangular boundary (default true).
#'
#' @param p.coords Matrix of simplex coordinates (3D) to plot with \code{points},
#'   one row per point. If a vector is given, this is converted to a matrix with
#'   a single row.
#' @param p.labels Labels of the points \code{p.coords}. If labels are undesired,
#'   set this to \code{NA}. The default is to use the row names of \code{p.coords}.
#'
#' @param p.labels.pos Which side to plot the \code{p.labels}. Values of \code{1},
#' \code{2}, \code{3} and \code{4}, respectively indicate positions below, to
#' the left of, above and to the right of the points. Can be either a single
#' integer (default 3) or a vector of integers.
#'
#' @param ... Further graphical parameters passed to points.
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
#' post <- compute_posterior(y, fs, prior)
#'
#' p.coords <- rbind(prior, post$marg)
#'
#' # Plot simplex with probability greater than 0.8 considered relatively
#' # certain
#' plot_simplex(v.cutoff = 0.8, p.coords = p.coords,
#'              p.labels = c("Prior", "Posterior"))
#'
#' # Plot the deviation of the posterior from the prior, this requires manually
#' # obtaining 2D coordinates
#' xy_prior <- project2D(as.vector(prior))
#' xy_post <- project2D(as.vector(post$marg))
#' arrows(x0 = xy_prior["x"], x1 = xy_post["x"],
#'        y0 = xy_prior["y"], y1 = xy_post["y"], length = 0.1)
#' @export
plot_simplex <- function(v.labels = c("Recrudescence", "Relapse", "Reinfection"),
                         v.cutoff = 0.5,
                         v.colours = c("yellow","purple","red"),
                         plot.tri = T,
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
#' Project three probabilities that sum to one (e.g., the marginal probabilities
#' of the 3Rs) onto the coordinates of a 2D simplex centred at the origin
#' (i.e., a triangle centred at (0,0) with unit sides). Each probability is
#' proportional to the distance between the point and the side opposite to the
#' corner corresponding to the probability.
#'
#' @param v A numeric vector of three numbers in zero to one that sum to one.
#'
#' @return A numeric vector of two coordinates that can be used to plot the
#'   probability vector `v` on the origin-centred 2D simplex (see
#'   [plot_simplex()]), where the top, left, and right vertices of the simplex
#'   correspond with the first, second and third entries of `v`
#'   respectively.
#'
#' @examples
#' project2D(v = c(0.75,0.20,0.05))
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







