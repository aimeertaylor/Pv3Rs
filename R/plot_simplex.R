#' Plots a 2D simplex
#'
#' Plots a 2D simplex, a triangle with unit sides centred at the origin, onto
#' which marginal posterior probabilities of relapse, reinfection and
#' recrudescence (or any other vector of three numbers in zero to one summing to
#' one) can be projected; see [project2D()] and examples below.
#'
#' @param v_labels A vector of labels that annotate vertices anticlockwise from
#'   top. If NULL (default), vertices are not annotated.
#'
#' @param v_cutoff An arbitrary number between 0.5 and 1 that separates regions
#'   of lower and higher probability. Beware the use of fixed cut-offs for
#'   probable recrudescence classification and probable reinfection
#'   classification; see [vignette("model-output", package = "Pv3Rs")].
#'
#' @param v_colours A vector of colours associated with the vertices
#'   anticlockwise from top; see example below.
#'
#' @examples
#' # Plot 2D simplex
#' plot_simplex()
#'
#' xy <- project2D(v = c("C" = 1, "L" = 0, "I" = 0))
#' points(x = xy["x"], xy["y"], pch = "C")
#' graphics::text(x = xy["x"], xy["y"], labels = "(1,0,0)", pos = 3)
#'
#' xy <- project2D(v = c("C" = 0, "L" = 1, "I" = 0))
#' points(x = xy["x"], xy["y"], pch = "L")
#' graphics::text(x = xy["x"], xy["y"], labels = "(0,1,0)", pos = 3)
#'
#' xy <- project2D(v = c("C" = 0, "L" = 0, "I" = 1))
#' points(x = xy["x"], xy["y"], pch = "I")
#' graphics::text(x = xy["x"], xy["y"], labels = "(0,0,1)", pos = 3)
#'
#'
#' # ==============================================================================
#' # Given data on an enrollment episode and a recurrence,
#' # compute the posterior probabilities of the 3Rs and plot the deviation of the
#' # posterior from the prior
#' # ==============================================================================
#'
#' # Some data:
#' y <- list(list(m1 = c("A", "C"), m2 = c("G", "T")), # Enrollment episode
#'           list(m1 = c("A"), m2 = c("G"))) # Recurrent episode
#'
#' # Some allele frequencies:
#' fs <- list(m1 = setNames(c(0.4, 0.6), c("A", "C")),
#'            m2 = setNames(c(0.2, 0.8), c("G", "T")))
#'
#' # A vector of prior probabilities:
#' prior <- array(c(0.2, 0.3, 0.5), dim = c(1,3),
#'                dimnames = list(NULL, c("C", "L", "I")))
#'
#' # Compute posterior probabilities
#' post <- compute_posterior(y, fs, prior)
#'
#' # Projev_cutoff marginal prior probabilities onto x and y coordinates:
#' xy_prior <- project2D(as.vector(prior))
#'
#' # Projev_cutoff marginal posterior probabilities onto x and y coordinates:
#' xy_post <- project2D(as.vector(post$marg))
#'
#' # Plot simplex with probability greater than 0.8 considered relatively
#' # certain
#' plot_simplex(colnames(post$marg), 0.8)
#'
#' # Plot the deviation of the posterior from the prior
#' arrows(x0 = xy_prior["x"], x1 = xy_post["x"],
#'        y0 = xy_prior["y"], y1 = xy_post["y"], length = 0.1)
#' @export
plot_simplex <- function(v_labels = NULL,
                         v_cutoff = 0.5,
                         v_colours = c("yellow","purple","red")) {

  # Define some constants:
  h <- sqrt(3)/2 # Height of equilateral triangle with unit sides
  r <- 1/sqrt(3) # Radius of circle encompassing triangle with unit sides
  k <- h-r # Distance from (0,0) to bottom of triangle with unit sides

  # Null plot
  plot(NULL, xlim = c(-0.6, 0.6), ylim = c(-(k + 0.1), r + 0.1), asp = 1,
       xaxt = "n", yaxt = "n", bty = "n",
       ylab = "", xlab = "")

  # Plot equilateral triangle:
  graphics::polygon(x = c(-0.5, 0.5, 0), y = c(-k, -k, r))

  # Annotate vertices:
  if (!is.null(v_labels)) {
    graphics::text(x = c(0, -0.5, 0.5),
                   y = c(r, -k, -k),
                   labels = v_labels,
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

  graphics::polygon(x = p1[,"x"], y = p1[,"y"], border = NA, col = grDevices::adjustcolor(v_colours[1], alpha.f = 0.35))
  graphics::polygon(x = p2[,"x"], y = p2[,"y"], border = NA, col = grDevices::adjustcolor(v_colours[2], alpha.f = 0.35))
  graphics::polygon(x = p3[,"x"], y = p3[,"y"], border = NA, col = grDevices::adjustcolor(v_colours[3], alpha.f = 0.35))


  # Delineate strong classification
  if(!is.null(v_cutoff)) {

    if(v_cutoff < 0.5 | v_cutoff > 1) stop("The v_cutoff should be a number between 0.5 and 1")

    # Project critical points to colour high probability regions (L for left, T
    # for top, R for right)
    xyTL_T <- project2D(c(v_cutoff,1-v_cutoff,0))
    xyTL_L <- project2D(c(1-v_cutoff,v_cutoff,0))
    xyTR_T <- project2D(c(v_cutoff,0,1-v_cutoff))
    xyTR_R <- project2D(c(1-v_cutoff,0,v_cutoff))
    xyLR_L <- project2D(c(0,v_cutoff,1-v_cutoff))
    xyLR_R <- project2D(c(0,1-v_cutoff,v_cutoff))

    # Regions with more than v_cutoff probability
    p1 <- rbind(xyT, xyTR_T, xyTL_T)
    p2 <- rbind(xyL, xyTL_L, xyLR_L)
    p3 <- rbind(xyR, xyLR_R, xyTR_R)

    graphics::polygon(x = p1[,"x"], y = p1[,"y"], border = NA, col = v_colours[1])
    graphics::polygon(x = p2[,"x"], y = p2[,"y"], border = NA, col = v_colours[2])
    graphics::polygon(x = p3[,"x"], y = p3[,"y"], border = NA, col = v_colours[3])
  }

}


#' Project 3D probability coordinates onto 2D simplex coordinates
#'
#' Project three probabilities that sum to one (e.g., the marginal probabilities
#' of the 3Rs) onto the coordinates of a 2D simplex centred at the origin
#' (i.e., a triangle centred at (0,0) with unit sides).
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







