#' Determine MOIs from unphased data
#'
#' @param y A sequence of lists, where each list contains the genetic data at
#'   each marker for one infection.
#'
#' @return Returns a vector of the (minimum) multiplicity of infection (MOI) for
#'   each infection.
#'
#' @examples
#' # two infections, three markers
#' y <- list(list(m1=c("A","B"), m2=c("A"), m3=c("C")),
#'           list(m1=c("B"), m2=c("B","C"), m3=c("A","B","C")))
#' determine_MOIs(y) # should be c(2, 3)
#'
#' @export
determine_MOIs <- function(y) {
  unname(sapply(y, function(x) max(sapply(x, length))))
}

#' Find all allele assignments for genotypes within the same infection
#'
#' Due to permutation symmetry of intra-infection genotypes, we only one
#' assignment for one of the markers whose number of alleles observed is equal
#' to the MOI.
#'
#' @param y.inf List of alleles observed across markers for genotypes within
#'   one infection.
#' @param gs.inf Vector of genotype names for genotypes within one infection.
#'
#' @return List of dataframes, one for each marker. The number of columns for
#'   each dataframe is the number of genotypes. The number of rows for each
#'   dataframe is the number of allele assignments for that marker.
#'
#' @examples
#' # three markers
#' y.inf <- list(m1=c("A","B"), m2=c("B","C","D"), m3=c("C"))
#' enumerate_alleles(y.inf, c("g1", "g2", "g3"))
#' # 6 assignments for m1 (BAA, ABA, BBA, AAB, BAB, ABB)
#' # 1 assignment for m2 (accounting for permutation symmetry)
#' # 1 assignment for m3 (CCC)
#'
#' @export
enumerate_alleles <- function(y.inf, gs.inf){
  # only one genotype
  if(length(gs.inf)==1) {
    return(lapply(y.inf, function(x) {
      x <- as.data.frame(x)
      colnames(x) <- gs.inf
      x
    }))
  }

  # find marker for fixed allele assignment
  n.uniq <- sapply(y.inf, length)
  MOI <- max(n.uniq)
  for(m.fix in names(n.uniq)) {
    if(n.uniq[m.fix] == MOI) break
  }

  comb_per_m <- list() # stores allele assignments for each marker
  for(m in names(n.uniq)) {
    if(m==m.fix) {
      comb_per_m[[m]] <- as.data.frame(t(y.inf[[m]]))
      colnames(comb_per_m[[m]]) <- gs.inf
    }
    else {
      # get all combinations
      combs <- expand.grid(rep(list(y.inf[[m]]), MOI), stringsAsFactors=F)
      # then filter out that are not consistent with observed data
      comb_per_m[[m]] <- combs[apply(combs, 1,
                                     function(row) length(unique(row)))==n.uniq[m],]
      colnames(comb_per_m[[m]]) <- gs.inf
    }
  }
  comb_per_m
}

#' Determine allele assignments that are consistent with each IBD pair
#'
#' Pre-process the allele assignments such that the enumeration required to
#' compte p(y|IBD) can be performed more efficiently. For each pair of
#' genotypes, find the rows of allele assignments that allow the pair of
#' genotypes to share the same allele. This function handles one marker only.
#'
#' @param al.df Dataframe of allele assignments (one assignment per
#'   row) across all genotypes (columns)
#'
#' @return List (hash table) that maps pairs of genotypes (represented as a
#' hyphen-separated string) to row numbers of \code{al.df} that have the
#' pair of genotypes sharing the same allele.
#'
#' @examples
#' al.df <- as.data.frame(matrix(c("A","A","A","B","B","A","B","B"), nrow=4))
#' colnames(al.df) <- c("g1", "g2")
#' allele_filter(al.df) # list where "g1-g2" maps to (2, 4)
#'
#' @export
allele_filter <- function(al.df) {
  stopifnot(class(al.df)[1] == "data.frame")
  allowed <- list()
  gs <- colnames(al.df)
  gs_count <- length(gs)
  for(g.i in tail(1:gs_count, -1)) {
    for(g.j in 1:(g.i-1)) {
      g1 <- gs[g.j]
      g2 <- gs[g.i]
      rows <- which(al.df[,g1] == al.df[,g2])
      allowed[[paste(g1, g2, sep="-")]] <- rows
      allowed[[paste(g2, g1, sep="-")]] <- rows
    }
  }
  allowed
}


#' Find all vectors of recurrence states compatible with relationship graph
#'
#' @param RG Relationship graph; see \code{\link{enumerate_RGs_alt}}.
#' @param gs_per_ts List of vectors of genotypes for each infection.
#'
#' @return Vector of strings (consisting of "C", "L", "I" for recrudescence,
#'   relapse, reinfection respectively) compatible with relationship graph.
#'
#' @examples
#' MOIs <- c(2, 2, 1)
#' RG <- enumerate_RGs_alt(MOIs, igraph=T)[[175]]
#' gs_per_ts <- split(paste0("g", 1:sum(MOIs)), rep(1:length(MOIs), MOIs))
#' # first recurrence cannot be recrudescence, second recurrence cannot be reinfection
#' plot_RG(RG, edge.curved=0.2)
#' compatible_rstrs(RG, gs_per_ts) # "LL" "IL" "LC" "IC"
#'
#' @export
compatible_rstrs <- function(RG, gs_per_ts) {
  infection_count <- length(gs_per_ts)
  n_recur <- infection_count-1
  r_by_recur <- lapply(1:n_recur, function(x) c("L"))
  clone.vec <- setNames(RG$clone.vec, paste0("g", 1:length(RG$clone.vec)))
  sib.vec <- RG$sib.vec

  for(i in 1:n_recur) {
    # can it be recrudescence
    recru <- TRUE
    for(g2 in gs_per_ts[[i+1]]) {
      if(!any(clone.vec[gs_per_ts[[i]]] == clone.vec[g2])) {
        recru <- FALSE
        break
      }
    }
    if(recru) r_by_recur[[i]] = c(r_by_recur[[i]], "C")

    # can it be reinfection
    reinf <- TRUE
    for(g2 in gs_per_ts[[i+1]]) {
      if(any(sib.vec[clone.vec[unlist(gs_per_ts[1:i])]] == sib.vec[clone.vec[g2]])) {
        reinf <- FALSE
        break
      }
    }
    if(reinf) r_by_recur[[i]] = c(r_by_recur[[i]], "I")
  }

  do.call(paste0, expand.grid(r_by_recur))
}

#' Convert IBD partition to a unique string for hashing
#'
#' This is used for building a hash table for p(y at marker m|IBD). Works up to
#' 10 genotypes. This seems to be marginally faster than
#' \code{\link{convert_IP_to_string}}.
#'
#' @param IP List containing vectors of genotype names corresponding to IBD
#'   clusters.
#' @param gs Vector containing all genotype names.
#'
#' @return Membership vector in string format.
#'
#' @export
hash.IP <- function(IP, gs) {
  ibd_vec <- setNames(gs, gs)
  ibd_i <- 0
  for(ibd_idx in order(sapply(IP, min))) {
    for(g in IP[[ibd_idx]]) ibd_vec[g] <- ibd_i
    ibd_i <- ibd_i + 1
  }
  paste(ibd_vec, collapse="")
}

#' Convert equivalence object to canonical form, i.e. each subset is sorted,
#' and the subsets are presented in lexicographical order.
#'
#' @param eq A list or equivalence object representing a partition.
#'
#' @return Equivalence object in canonical form.
#'
#' @examples
#' parts <- partitions::listParts(4)
#' lapply(parts, eq_canonical)
#'
#' @export
eq_canonical <- function(eq) {
  inner.sorted <- lapply(eq, sort)
  out <- unname(inner.sorted[order(sapply(inner.sorted,'[',1))])
  class(out) <- c(class(out), "equivalence")
  out
}

#' Inverse of partitions::vec_to_eq.
#'
#' @noRd
part_to_vec <- function(part, count) {
  i <- 1
  vec <- setNames(rep(NA, count), unlist(part))
  for(s in part) {
    for(num in s) {
      vec[num] <- i
    }
    i <- i+1
  }
  return(vec)
}

#' Given a vector, return a list consisting of the vector itself and all possible
#' pairs of disjoint subsets whose union covers the vector.
#'
#' @noRd
split_two <- function(s) {
  n <- length(s)
  if(n == 1) return(list(list(s)))
  masks <- 2^(1:n-1)
  c(list(list(s)),
    lapply(seq(1,2^n-2,2), function(u) list(s[bitwAnd(u, masks)!=0],
                                            s[bitwAnd(u, masks)==0])))
}

#' Extract information from weights matrix of an \code{igraph} object.
#'
#' Return a list containing (a) list of clonal units and (b) a list of sibling
#' clusters. Each clonal unit stores a vector of genotype names. Each sibling
#' cluster stores a vector of clonal unit names.
#'
#' @noRd
igraph_to_parts <- function(RG) {
  gs_count <- igraph::vcount(RG)
  W <- RG[]

  clone.part <- 1:gs_count
  for(i in tail(1:gs_count, -1)) {
    for(j in 1:(i-1)) {
      if(W[i,j]==1) {
        clone.part[i] <- clone.part[j]
        break
      }
    }
  }
  clone.list <- split(igraph::V(RG)$name, clone.part)
  n.clone <- length(clone.list)
  clone.names <- paste0("c", 1:n.clone)

  sib.part <- 1:n.clone
  for(i in tail(1:n.clone, -1)) {
    row <- W[clone.list[[i]][1]]
    for(j in 1:(i-1)) {
      if(row[clone.list[[j]][1]]==0.5) {
        sib.part[i] <- sib.part[j]
        break
      }
    }
  }
  sib.list <- split(clone.names, sib.part)
  sib.names <- paste0("s", 1:length(sib.list))

  return(list(clone=setNames(clone.list, clone.names),
              sib=setNames(sib.list, sib.names)))
}
