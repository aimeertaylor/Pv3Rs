#' Enumerate IBD partitions
#'
#' Wrapper function around \code{\link[partitions]{listParts}}. Modifies the
#' output of \code{\link[partitions]{listParts}} to make genotype names
#' explicit.
#'
#' @param MOIs Numeric vector of multiplicities of infections (MOIs) that gets
#'   reduced to its sum, which must be less than or equal to 11.
#'
#' @examples
#' enumerate_IPs(3)
#'
#' @section To-do:
#' Convert to use \code{\link[partitions]{setparts}} whose output format takes
#' less memory, such that larger \code{MOI}s can be dealt with.
#'
#' @export
enumerate_IPs <- function(MOIs) {
  genotype_count <- sum(MOIs)
  # previously hard coded a hard limit, now removed
  # if (genotype_count > 11) stop("Too many IBD partitions to enumerate.")
  ps <- partitions::listParts(genotype_count)

  IPs <- lapply(ps, function(p) {
    out <- lapply(p, function(x) paste0("g", x))
    class(out) <- c(class(out), "equivalence")
    return(out)
  })

  return(IPs)
}
