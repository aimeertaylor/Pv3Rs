#' Generates \code{n} population lineages
#'
#'   Letters represent lineages with repeats if more than 26 lineages are
#'   required.
#'
#' @param n Number of lineages to generate.
#'
#' @examples
#' rev(generate_lineages(28))
#'
#' @export
generate_lineages <- function(n) {

  max_reps <- ceiling(n/length(letters))
  lineages <- c()
  for (reps in max_reps:1) {
    new_lineages <- apply(matrix(rep(letters, each = reps), nrow = reps), 2, paste0, collapse = "")
    lineages <- c(new_lineages, lineages)
  }
  return(lineages[1:n])
}
