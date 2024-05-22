# ==============================================================================
# Function to generate alleles for three half sibs by enumeration (naive)
# ==============================================================================
generate_halfsib_alleles <- function(n_alleles){

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

  # End of function
  return(child_alleles)
}
