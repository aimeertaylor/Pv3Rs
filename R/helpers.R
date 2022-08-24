# inverse of partitions::vec_to_eq
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

# given a vector, return a list consisting of the vector itself and all possible
# pairs of disjoint subsets whose union covers the vector
split_two <- function(s) {
  n <- length(s)
  if(n == 1) return(list(list(s)))
  masks <- 2^(1:n-1)
  c(list(list(s)),
    lapply(seq(1,2^n-2,2), function(u) list(s[bitwAnd(u, masks)!=0],
                                            s[bitwAnd(u, masks)==0])))
}

#' Convert equivalence object to canonical form (ie sort lists)
#'
#' @export
eq_canonical <- function(eq) {
  inner.sorted <- lapply(eq, sort)
  out <- unname(inner.sorted[order(sapply(inner.sorted,'[',1))])
  class(out) <- c(class(out), "equivalence")
  out
}

# weights matrix to list containing (a) list of clonal units and (b) a list of
# sibling clusters. Each clonal unit stores a vector of genotype names. Each
# sibling cluster stores a vector of clonal unit names.
RG_to_parts <- function(RG) {
  gs_count <- igraph::vcount(RG)

  clone.part <- 1:gs_count
  for(i in 2:gs_count) {
    for(j in 1:(i-1)) {
      if(RG[i,j]==1) {
        clone.part[i] <- clone.part[j]
        break
      }
    }
  }
  clone.list <- split(igraph::V(RG)$name, clone.part)
  n.clone <- length(clone.list)
  clone.names <- paste0("c", 1:n.clone)

  sib.part <- 1:n.clone
  for(i in 2:n.clone) {
    row <- RG[clone.list[[i]][1]]
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
