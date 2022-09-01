#' Computes the likelihood of all relationship graphs
#'
#' @param MOIs A numeric vector specifying, for each infection, the number of
#'   distinct parasite genotypes, a.k.a. the multiplicity of infection (MOI).
#' @param fs List of frequencies of each allele at each marker.
#' @param alleles_per_m List of allele assignments for each marker (see the
#'   \code{al.mat} parameter to \code{\link{allele.filter}}).
#' @param als_by_edge List of outputs from \code{\link{allele_filter}} for
#'   each marker.
#'
#' @return List of all relationship graphs, including their log-likelihoods as
#'   a variable under each relationship graph.
#'
#' @examples
#' # 1 marker, 2 infections, MOIs = 2, 1
#' MOIs <- c(2, 1)
#' fs <- list(m1=setNames(c(0.4, 0.6), c("A", "B")),
#'            m2=setNames(c(0.2, 0.8), c("C", "D")))
#' al_df1 <- as.data.frame(matrix(c("A","B","B"), nrow=1))
#' al_df2 <- as.data.frame(matrix(c("C","D","C",
#'                                   "D","C","C"), nrow=2, byrow=T))
#' colnames(al_df1) <- colnames(al_df2) <- paste0("g", 1:3)
#' alleles_per_m <- list(m1=al_df1, m2=al_df2)
#' als_by_edge <- lapply(alleles_per_m, allele_filter)
#' RGs <- RG_inference(MOIs, fs, alleles_per_m, als_by_edge)
#'
#' @export
RG_inference <- function(MOIs, fs, alleles_per_m, als_by_edge) {
  for(al.df in alleles_per_m) stopifnot(class(al.df)[1] == "data.frame")
  ms <- names(fs)
  log_fs <- lapply(fs, log) # allele log-frequencies
  # number of allele assignments for each marker
  a_sizes <- lapply(alleles_per_m, nrow)
  # hash tables to store p(marker m observed data|IBD)
  IP_lookups <- setNames(lapply(ms, function(m) new.env(hash=T,
                                                        parent=emptyenv())), ms)
  RGs <- enumerate_RGs_alt(MOIs, igraph=FALSE)
  gs <- paste0("g", 1:sum(MOIs))

  RG_i <- 0
  n.RG <- length(RGs)
  pbar <- txtProgressBar(min=0, max=n.RG) # min=0 in case n.RG is 1
  writeLines(paste("\nComputing log p(Y|RG) for", n.RG, "RGs"))

  for(RG in RGs) {
    # TODO: for every pair of clones, if they cannot have same genotype, set logp to -Inf

    IPs_RG <- enumerate_IPs_RG(RG)
    n.IPs <- length(IPs_RG)
    IP_logps <- matrix(0, n.IPs, length(log_fs)) # stores p(y for each marker m | IP)
    colnames(IP_logps) <- ms
    IP_i <- 1
    for(IP in IPs_RG) {
      IP_str <- hash.IP(IP, gs)
      # p(y | IP) is sum of p(a | IP)
      # p(a | IP) can be decomposed by marker
      for(m in ms) {
        res <- IP_lookups[[m]][[IP_str]]
        if(is.null(res)) {
          rows <- 1:a_sizes[[m]]
          # find rows of alleles_per_m[[m]]
          for(g_seq in IP) {
            if(length(g_seq) > 1) {
              for(edge in paste(head(g_seq, -1), tail(g_seq, -1), sep='-')) {
                rows <- intersect(rows, als_by_edge[[m]][[edge]])
              }
            }
          }
          if(length(rows) == 0) res <- -Inf
          else {
            log_fs_mat <- apply(alleles_per_m[[m]][rows, sapply(IP, '[', 1)], 2,
                                function(x) log_fs[[m]][x])
            if(is.null(dim(log_fs_mat))) res <- sum(log_fs_mat) # only one allele assignment
            else res <- matrixStats::logSumExp(rowSums(log_fs_mat))
          }
          IP_lookups[[m]][[IP_str]] <- res
        }
        IP_logps[IP_i,m] <- res
      }
      IP_i <- IP_i + 1
    }

    RG_i <- RG_i + 1
    setTxtProgressBar(pbar, RG_i)

    RGs[[RG_i]]$logp <- sum(matrixStats::colLogSumExps(IP_logps)) - length(ms)*log(n.IPs)

  }
  writeLines("")
  RGs
}