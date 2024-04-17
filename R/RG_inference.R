#' Computes the likelihood of all relationship graphs
#'
#' The likelihood of a relationship graph (RG) can be decomposed over markers,
#' as we assume that marker data is conditionally independent across markers
#' given the RG. The likelihood of a RG for one marker is then found by
#' integrating over all possible IBD (identity-by-descent) partitions that are
#' consistent with the RG. The probability distribution of these IBD partitions
#' are always uniform. Since different RGs may use the same IBD, the likelihood
#' of IPs for each marker are stored in hash tables for improved computational
#' efficiency.
#'
#' @param MOIs A numeric vector specifying, for each infection, the number of
#'   distinct parasite genotypes, a.k.a. the multiplicity of infection (MOI).
#' @param fs List of allele frequencies as vectors. Names of the list must
#'   match with the marker names in `alleles_per_m`.
#' @param alleles_per_m List of allele assignments as dataframes for each
#'   marker. Each column corresponds to a genotype and each row corresponds to
#'   an allele assignment.
#'
#' @return List of all relationship graph objects, including their
#'   log-likelihoods as a variable under each object.
#'
#' @examples
#' # 1 marker, 2 infections, MOIs = 2, 1
#' MOIs <- c(2, 1)
#' fs <- list(
#'   m1 = setNames(c(0.4, 0.6), c("A", "B")),
#'   m2 = setNames(c(0.2, 0.8), c("C", "D"))
#' )
#' al_df1 <- as.data.frame(matrix(c("A", "B", "B"), nrow = 1))
#' al_df2 <- as.data.frame(matrix(c(
#'   "C", "D", "C",
#'   "D", "C", "C"
#' ), nrow = 2, byrow = T))
#' colnames(al_df1) <- colnames(al_df2) <- paste0("g", 1:3)
#' alleles_per_m <- list(m1 = al_df1, m2 = al_df2)
#' RGs <- RG_inference(MOIs, fs, alleles_per_m)
#'
#' @export
RG_inference <- function(MOIs, fs, alleles_per_m) {
  # check that allele assignments are given as data frames
  for (al.df in alleles_per_m) stopifnot(class(al.df)[1] == "data.frame")

  if(sum(MOIs) > 10) warning(
    "Total MOI > 10 may lead to high memory use", immediate=T
  )

  ms <- names(alleles_per_m) # marker names
  log_fs <- lapply(fs, log) # allele log-frequencies
  # number of allele assignments for each marker
  a_sizes <- lapply(alleles_per_m, nrow)

  # hash tables to store p(marker m observed data | IBD partition)
  IP_lookups <- setNames(lapply(ms, function(m) {
    new.env(
      hash = T,
      size = multicool::Bell(sum(MOIs)),
      parent = emptyenv()
    )
  }), ms)
  # enumerate all relationship graphs
  RGs <- enumerate_RGs(MOIs, igraph = FALSE)
  gs <- paste0("g", 1:sum(MOIs)) # genotype names

  RG_i <- 0
  n.RG <- length(RGs)
  pbar <- txtProgressBar(min = 0, max = n.RG) # min=0 in case n.RG is 1
  writeLines(paste("Computing log p(Y|RG) for", n.RG, "RGs"))

  for (RG in RGs) { # for each relationship graph
    # check if there is some cell of clones cannot have the same genotype
    incompatible <- FALSE
    for (clones in RG$clone) { # for each cell of clones
      if (length(clones) == 1) next
      for (assignment_df in alleles_per_m) { # for each marker
        # check if there is no assignment where all of their alleles are equal
        if (!any(apply(
          assignment_df[clones], 1,
          function(row) length(unique(na.omit(row))) <= 1
        ))) {
          incompatible <- TRUE
          break
        }
      }
      if (incompatible) break
    }
    # set logp to -Inf in case of any clonal edges that are impossible
    if (incompatible) {
      RG_i <- RG_i + 1
      setTxtProgressBar(pbar, RG_i)
      RGs[[RG_i]]$logp <- -Inf
      next
    }

    # enumerate all IBD partitions consistent with relationship graph
    IPs_RG <- enumerate_IPs_RG(RG)
    n.IPs <- length(IPs_RG) # number of IPs (IBD partitions)
    # stores log p(y for each marker m | IP)
    IP_logps <- matrix(0, n.IPs, length(ms))
    colnames(IP_logps) <- ms

    IP_i <- 1
    for (IP in IPs_RG) {
      IP_str <- hash.IP(IP, gs)
      # p(y | IP) is sum of p(a | IP), where p(a | IP) is decomposed by marker
      for (m in ms) { # for each marker m
        IP.prob <- IP_lookups[[m]][[IP_str]] # look up p(a_m | IP) in hash table
        if (is.null(IP.prob)) {
          logf <- log_fs[[m]] # allele frequencies for marker m
          logf_sums <- apply(
            alleles_per_m[[m]], 1,
            function(row) { # for each allele assignment
              logf_sum <- 0 # sum of log allele frequencies across IBD cells
              for (g_seq in IP) { # for genotypes in each IBD cell
                repr <- NA # allele representative for current cell
                for (allele in row[g_seq]) { # for allele of each genotype
                  if (!is.na(allele)) { # allele is not missing
                    if (is.na(repr)) {
                      # record as representative if it does not exist yet
                      repr <- allele
                    } else if (repr != allele) {
                      # return NA if IBD cell has different alleles
                      return(NA)
                    }
                  }
                }
                if (!is.na(repr)) logf_sum <- logf_sum + logf[repr]
              }
              return(logf_sum)
            }
          )
          # p(y | IP) is sum of p(a | IP), but compute in log space
          IP.prob <- matrixStats::logSumExp(logf_sums, na.rm = T)
          IP_lookups[[m]][[IP_str]] <- IP.prob # store in hash table
        }
        IP_logps[IP_i, m] <- IP.prob
      }
      IP_i <- IP_i + 1
    }

    RG_i <- RG_i + 1
    setTxtProgressBar(pbar, RG_i)

    # p(y | RG) = sum of p(y | IP) p(IP | RG), where p(IP | RG) is uniform
    RGs[[RG_i]]$logp <- (sum(matrixStats::colLogSumExps(IP_logps))
    - length(ms) * log(n.IPs))
  }
  writeLines("")
  RGs
}
