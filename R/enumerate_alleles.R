#' Find all allele assignments for genotypes within the same infection
#'
#' Note that this function is not tested for input with alleles mixed with NA.
#'
#' @param y.inf List of alleles observed across markers for genotypes within one
#'   infection. Alleles must not be repeated for the same marker.
#' @param gs.inf Vector of genotype names for genotypes within one infection.
#' @param use.sym Boolean for permutation symmetry is exploited as a
#'   computational shortcut. Due to permutation symmetry of intra-infection
#'   genotypes, we can fix a single assignment for one of the markers whose
#'   number of alleles observed is equal to the MOI (consider it the anchor)
#'   and permute the rest, discarding combinations that under-represent the
#'   observed marker diversity. The default behaviour is to use this symmetry
#'   such that less assignments have to be considered.
#'
#' @return List of dataframes, one for each marker. The columns correspond to
#'   genotypes, while the rows correspond to allele assignments for a marker.
#'
#' @examples
#' # 3 markers
#' y.inf <- list(m1 = c("A", "B"), m2 = c("B", "C", "D"), m3 = c("C"))
#' enumerate_alleles(y.inf, c("g1", "g2", "g3"))
#' # 6 assignments for m1 (BAA, ABA, BBA, AAB, BAB, ABB)
#' # 1 assignment for m2 (accounting for permutation symmetry)
#' # 1 assignment for m3 (CCC)
#'
#' @noRd
enumerate_alleles <- function(y.inf, gs.inf, use.sym = TRUE) {
  # edge case: only one genotype
  if (length(gs.inf) == 1) {
    return(lapply(y.inf, function(x) {
      x <- as.data.frame(x)
      colnames(x) <- gs.inf
      x
    }))
  }

  for (y.m in y.inf) {
    stopifnot("Allele repeats are not allowed"=length(unique(y.m))==length(y.m))
  }

  marker_names <- names(y.inf)
  # find single marker for fixed allele assignment
  # must be a marker whose number of observed alleles = MOI
  # arbitrarily take the first such marker to be the 'anchor' marker
  y.lens <- sapply(y.inf, length)
  MOI <- length(gs.inf)

  m.fix <- NA
  for (m in marker_names) {
    if (y.lens[m] == MOI) {
      m.fix <- m
      break
    }
  }
  if(is.na(m.fix)) use.sym <- F

  comb_per_m <- list() # stores allele assignments for each marker
  for (m in marker_names) {
    if (m == m.fix && use.sym) { # only one fixed assignment for 'anchor' marker
      comb_per_m[[m]] <- as.data.frame(t(y.inf[[m]]))
      colnames(comb_per_m[[m]]) <- gs.inf
    } else {
      # get all combinations of allele assignments
      n.alleles <- length(y.inf[[m]])
      combs <- expand.grid(rep(list(y.inf[[m]]), MOI), stringsAsFactors = F)
      # remove assignments that under-represent the observed marker diversity
      res <- apply(
        combs, 1,
        function(row) {
          length(unique(row)) == n.alleles
        }
      )
      comb_per_m[[m]] <- dplyr::distinct(combs[res,])
      colnames(comb_per_m[[m]]) <- gs.inf
    }
  }
  comb_per_m
}
