part.list <- list()
PART_MAX <- 10
for(i in 1:PART_MAX) {
  part.list[[i]] <- partitions::setparts(i)
}

#' Enumerate transitive relationship graphs, alternate version
#'
#' This alternate version of \code{enumerate_RGs} is based on generating all
#' set partitions. This is done twice in a nested fashion: once for determining
#' clonal relationships, once for determining sibling relationships. A list of
#' all partitions for sets of difference sizes is pre-computed.
#' @export
enumerate_RGs_alt <- function(MOIs) {

  # Check MOIs are positive whole numbers
  if (!all(is.wholenumber(MOIs)) | any(MOIs < 1)) stop("MOIs should be positive integers")

  # Hard code relationship types to satisfy the is.transitive function
  relationship_types <- c(stranger = 0, sibling = 0.5, clone = 1)
  intra_relationship_types <- relationship_types[setdiff(names(relationship_types), "clone")]
  infection_count <- length(MOIs) # Number of time points
  gs_count <- sum(MOIs) # Number of genotypes

  # remove restriction for now as this method should scale better
  # if (gs_count > 6 | infection_count > 3) stop("Sorry, too many RGs")
  # if (gs_count <= 1) stop("Sorry, no RGs")

  gs <- paste0("g", 1:gs_count) # Genotype names
  ts_per_gs <- rep(1:infection_count, MOIs)
  gnum_per_ts <- split(1:gs_count, ts_per_gs)
  gstarts <- head(c(1, cumsum(MOIs) + 1), -1) # first genotype of each infection

  RG_i <- 0
  RGs <- list()

  # find all possible clonal relationships
  # initialise with lexicographically 'smallest' possible clone partition
  clone_init <-  unlist(lapply(MOIs, function(x) 1:x))
  clone_part <- clone_init
  # prefix_max stores largest element of each prefix of clone_part (excluding current index)
  prefix_max <- rep(0, gs_count)
  for(i in 1:gs_count) {
    if(i == 1) next
    prefix_max[i] <-  max(prefix_max[i-1], clone_part[i-1])
  }
  clone_list <- list(clone_part)
  clone_i <- 1
  while(clone_part[gs_count] != gs_count) {
    # find genotype to increment
    for(i in gs_count:1) {
      if(clone_part[i] <= prefix_max[i]) break
    }
    start_i <- i
    infection <- ts_per_gs[i]
    # clonal units to avoid such that there are no clonal relationships within an infection
    avoid <- head(clone_part[gstarts[infection]:i], -1)
    candidate <- clone_part[i] + 1
    first <- TRUE
    while(i <= gs_count & ts_per_gs[i] == infection) { # while in current infection
      while(candidate %in% avoid) {
        candidate <- candidate + 1
      }
      clone_part[i] <- candidate
      if(first) {
        avoid <- c(avoid, candidate)
        candidate <- 1
        first <- FALSE
      } else candidate <- candidate + 1
      i <- i+1
    }
    if(i <= gs_count) {
      # complete with lexicographically 'smallest' possible clone partition
      clone_part[i:gs_count] <- clone_init[i:gs_count]
    }
    for(j in (start_i+1):gs_count) {
      prefix_max[j] <-  max(prefix_max[j-1], clone_part[j-1])
    }
    clone_i <- clone_i + 1
    clone_list[[clone_i]] <- clone_part
  }

  # Count number of valid graphs and create progress bar
  n.RG <- sum(sapply(clone_list, function(x) ncol(part.list[[max(x)]])))
  pbar <- txtProgressBar(min=0, max=n.RG) # min=0 in case n.RG is 1
  writeLines(paste("\nnumber of valid graphs is", n.RG))

  for(i in 1:clone_i) {
    n.units <- max(clone_list[[i]])
    units <- split(1:gs_count, clone_list[[i]])

    # given clonal relationships, generate all compatible sibling relationships
    if (n.units > 10) stop("Currently supports up to 10 genotypes")
    sib.parts <- part.list[[n.units]]
    n.parts <- ncol(sib.parts)
    for(j in 1:n.parts) {
      adj_all <- array(0, dim = rep(gs_count, 2), dimnames = list(gs,gs))

      n.sib.units <- max(sib.parts[,j])
      sib.units <- split(1:n.units, sib.parts[,j])
      for(k in 1:n.sib.units) {
        selection <- unlist(units[sib.units[[k]]])
        adj_all[selection, selection] <- 0.5;
      }
      for(u in 1:n.units) {
        adj_all[units[[u]], units[[u]]] <- 1;
      }
      diag(adj_all) <- 0

      # Convert adjacency matrix into an igraph item
      RG <- igraph::graph_from_adjacency_matrix(adjmatrix = adj_all,
                                                mode = "lower",
                                                diag = F,
                                                weighted = T)
      # Add a time-point attribute
      RG <- igraph::set_vertex_attr(RG, "group", value = ts_per_gs)

      RG_i <- RG_i + 1
      RGs[[RG_i]] <- RG
      setTxtProgressBar(pbar, RG_i)
    }
  }

  return(RGs)
}
