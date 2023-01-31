#' Compute the probabilities of partitions
#'
#' Computes the probabilities of all possible partitions of k balls given each
#' pair of balls has probability f of being in the same box (i.e. there are 1/f
#' equally probable boxes). Can be used to compute the probability of IBD
#' partitions among parasites who are all equally related to one another with
#' relatedness f. Can also be used to inefficiently solve the birthday problem
#' provided that there are not too many balls.
#'
#' @param k Number of balls to partition
#' @param f Probability that two balls are in the same box
#'
#' @examples
#'
#' compute_pr_partitions(3, f = 0.5)
#' compute_pr_partitions(3, f = 1)
#' compute_pr_partitions(3, f = 0)
#'
#' # Using compute_pr_partition to solve the birthday problem with...
#' k <- 5 # people
#' days <- 365 # days in the year
#'
#' # birthday problem:
#' stats::pbirthday(n = k, classes = days, coincident = 2)
#'
#' # birthday problem computed using compute_pr_partition:
#' pr_parts <- compute_pr_partitions(k, f = 1 / days)
#' 1 - tail(pr_parts, 1)
compute_pr_partitions <- function(k, f) {
  # Enumerate partitions
  Ps <- partitions::setparts(x = k)
  colnames(Ps) <- apply(Ps, 2, paste, collapse = "")

  box_count <- apply(Ps, 2, function(x) length(unique(x))) # number of filled boxes
  max_box_count <- k # maximum number of filled boxes possible
  new_box_count <- box_count - 1 # number of newly filled boxes conditional on first being filled

  # if all relationships have the same f...
  # list of vectors of (1-f) (1-2f)...
  join_new_box_probs <- sapply(new_box_count, function(x, y) {
    z <- 1 - (0:x) * y # z is a vector of (1-f) (1-2f) etc. some of whose elements may be negative
    z_non_neg <- sapply(z, function(x) max(0, x)) # replace any negative element of z with zero
    prod(z_non_neg) # Compute (1-f)*(1-2f)*(1-3f) etc.
  }, f)

  join_old_box_probs <- f^(max_box_count - box_count)
  overall_partition_probs <- join_new_box_probs * join_old_box_probs
  return(overall_partition_probs)
}
