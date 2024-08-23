#' List all allelic draws for three half siblings
#'
#' Given a specified number of alleles for a single marker,
#' `generate_halfsib_alleles()` enumerates all the ways three half
#' siblings can draw alleles from their respective parents by firstly
#' enumerating all allelic combinations for the parents, and by secondly
#' enumerating all the inheritable combinations for the children.
#'
#' @param n_alleles Positive whole number specifying a per-marker number of
#'   alleles, otherwise known as marker cardinality.
#'
#' @returns A character matrix. Each column is an individual. Each row is a
#'   possible allelic draw. Alleles are represented by the first `n_alleles`
#'   letters of the latin alphabet.
#' @examples
#' enumerate_halfsib_alleles(3)
#'
#' @export
enumerate_halfsib_alleles <- function(n_alleles){


  # Check positive integer
  if(!(n_alleles > 0 & is.wholenumber(n_alleles))) {
    stop ("n_alleles should be a positive wholenumber")
  }

  # Per-marker alphabet
  alleles <- LETTERS[1:n_alleles]

  # Generate all permutations of parental alleles
  parental_alleles <- as.matrix(expand.grid(p1 = alleles,
                                            p2 = alleles,
                                            p3 = alleles,
                                            stringsAsFactors = FALSE))

  # Generate all permutations of inheritance
  child_parents <- as.matrix(expand.grid(child12 = c("p1", "p2"),
                                         child13 = c("p1", "p3"),
                                         child23 = c("p2", "p3"),
                                         stringsAsFactors = FALSE))

  # Extract inherited alleles
  child_alleles <- do.call(rbind, apply(child_parents, 1, function(ps) {
    parental_alleles[,ps]
  }, simplify = F))

  # Rename
  colnames(child_alleles) <- paste("Half sibling", 1:3)

  # End of function
  return(child_alleles)
}
