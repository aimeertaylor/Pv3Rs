#' Plots a 2D simplex
#'
#' Plots a 2D simplex, a triangle with unit sides centered at the origin, to
#' which marginal posterior probabilities of the 3Rs can be added.
#'
#' @param v_labels A vector of labels with which to annotate vertices. If NULL,
#'   vertices are not annotated.
#'
#' @examples
#' # Plot 2D simplex
#' plot_simplex(v_labels = c("(1,0,0)", "(0,0,1)", "(0,1,0)"))
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
#' # Project marginal prior probabilities onto x and y coordinates:
#' xy_prior <- project2D(as.vector(prior))
#'
#' # Project marginal posterior probabilities onto x and y coordinates:
#' xy_post <- project2D(as.vector(post$marg))
#'
#' # Plot the deviation of the posterior from the prior
#' arrows(x0 = xy_prior["x"], x1 = xy_post["x"],
#'        y0 = xy_prior["y"], y1 = xy_post["y"], length = 0.1)
#' @export
plot_simplex <- function(v_labels) {

  # Define some constants:
  h <- sqrt(3)/2 # Height of equilateral triangle with unit sides
  r <- 1/sqrt(3) # Radius of circle encompassing triangle with unit sides
  k <- h-r # Distance from (0,0) to bottom of triangle with unit sides

  # Null plot
  plot(NULL, xlim = c(-0.6, 0.6), ylim = c(-(k + 0.1), r + 0.1), pty = 's',
       xaxt = "n", yaxt = "n", bty = "n",
       ylab = "", xlab = "")

  # Plot equilateral triangle:
  polygon(x = c(-0.5, 0.5, 0), y = c(-k, -k, r))

  # Annotate vertices:
  if (!is.null(v_labels)) {
    text(x = c(-0.5, 0.5, 0), y = c(-k, -k, r), labels = v_labels, pos = c(1,1,3))
  }

}


#' Project 3D probability coordinates onto 2D simplex coordinates
#'
#' Project three probabilities that sum to one (e.g., the marginal probabilities
#' of the 3Rs) onto the coordinates of a 2D simplex centered at the origin
#' (i.e., a triangle centered at (0,0) with unit sides).
#'
#' @param v A numeric vector of three probabilities that sum to one.
#' @return A numeric vector of two coordinates that can be used to plot the
#'   probability vector on the origin-centered 2D simplex.
#' @examples
#' project2D(v = c(0.75,0.20,0.05))
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







