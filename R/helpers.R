#' Test if x is an integer number.
#'
#' Internal function copied from [base::integer()]. It is differs from
#' [base::is.integer()] because [base::is.integer()] tests if x is an integer
#' type (for compatibility with C and Fortran code).
#'
#' @noRd
is_wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}

#' Find all vectors of recurrence states compatible with relationship graph
#'
#' Finds all possible recurrence states for each recurrence compatible with the
#' relationship graph, then takes the Cartesian product to get all vectors of
#' recurrence states. For a recurrence to be a recrudescence, all edges between
#' the recurrent infection and the immediately preceding infection must be
#' clonal edges. For a recurrence to be a reinfection, all edges between the
#' recurrent infection and any preceding infection must be stranger edges. All
#' recurrences may possibly be relapses.
#'
#' @param RG Relationship graph; see \code{\link{enumerate_RGs}}.
#' @param gs_per_ts List of vectors of genotypes for each infection.
#'
#' @return Vector of strings (consisting of "C", "L", "I" for recrudescence,
#'   relapse, reinfection respectively) compatible with relationship graph.
#'
#' @examples
#' MOIs <- c(2, 2, 1)
#' RG <- enumerate_RGs(MOIs, igraph = TRUE)[[175]]
#' gs_per_ts <- split(paste0("g", 1:sum(MOIs)), rep(1:length(MOIs), MOIs))
#' # 1st recurrence can't be recrudescence, 2nd recurrence can't be reinfection
#' plot_RG(RG, edge.curved = 0.2)
#' compatible_rstrs(RG, gs_per_ts) # "LL" "IL" "LC" "IC"
#'
#' @noRd
compatible_rstrs <- function(RG, gs_per_ts) {
  infection_count <- length(gs_per_ts)
  n_recur <- infection_count - 1

  # `r_by_recur` is a vector storing possible recur. states for each recurrence
  r_by_recur <- lapply(1:n_recur, function(x) c("L")) # relapse always possible
  # prepend clonal partition vector with 'g'
  clone.vec <- stats::setNames(RG$clone.vec, paste0("g", 1:length(RG$clone.vec)))
  sib.vec <- RG$sib.vec

  for (i in 1:n_recur) { # for each recurrence
    # can it be recrudescence
    recru <- TRUE
    for (g2 in gs_per_ts[[i + 1]]) { # for each genotype in recurrent infection
      if (!any(clone.vec[gs_per_ts[[i]]] == clone.vec[g2])) {
        # cannot be recrudescence if it is not in the same clonal unit as any
        # genotype in the immediately preceding infection
        recru <- FALSE
        break
      }
    }
    if (recru) r_by_recur[[i]] <- c(r_by_recur[[i]], "C")

    # can it be reinfection
    reinf <- TRUE
    for (g2 in gs_per_ts[[i + 1]]) { # for each genotype in recurrent infection
      if (any(sib.vec[clone.vec[unlist(gs_per_ts[1:i])]]
      == sib.vec[clone.vec[g2]])) {
        # cannot be reinfection if any genotype from the any preceding
        # infection is in the same sibling unit
        reinf <- FALSE
        break
      }
    }
    if (reinf) r_by_recur[[i]] <- c(r_by_recur[[i]], "I")
  }

  # take cartesian product
  do.call(paste0, expand.grid(r_by_recur))
}

#' Convert IBD partition to a unique string for hashing
#'
#' This is used for building a hash table for p(y at marker m|IBD).
#'
#' @param IP List containing vectors of genotype names with each vector
#' corresponding to an IBD cell.
#' @param gs Vector containing all genotype names.
#'
#' @return String where the integers in the IBD membership vector have been
#'   converted to ASCII characters.
#'
#' @examples
#' gs <- paste0("g", 1:3)
#' IP1 <- list(c("g1", "g3"), c("g2"))
#' IP2 <- list(c("g2"), c("g3", "g1"))
#' hash1 <- hash_IP(IP1, gs)
#' hash2 <- hash_IP(IP2, gs)
#' hash1 == hash2 # TRUE, even though the order is different
#'
#' @noRd
hash_IP <- function(IP, gs) {
  ibd_vec <- stats::setNames(gs, gs)
  ibd_i <- 1
  # use of `order` on the 'min' genotype name of each IBD cell ensures that the
  # for loop runs in the same order even if genotype names within an IBD cell
  # are arranged differently
  for (ibd_idx in order(sapply(IP, min))) {
    for (g in IP[[ibd_idx]]) ibd_vec[g] <- ibd_i
    ibd_i <- ibd_i + 1
  }
  intToUtf8(ibd_vec)
}

#' Partition a vector into at most two subvectors
#'
#' @param s Vector to be split.
#'
#' @return
#' Given a vector with no repeats, returns a list consisting of
#' \itemize{
#'   \item a list that contains the original vector as its only element
#'   \item other lists where each list contains two disjoint vectors whose
#'     union covers the vector. All possible unordered pairs are included.
#' }
#'
#' @examples
#' gs <- paste0("g", 1:3)
#' # 4 possibilities in total
#' # either all genotypes in one vector (1 possibility)
#' # or 2 genotypes in one vector and the last in one vector (3 possibilties)
#' split_two(gs)
#'
#' @noRd
split_two <- function(s) {
  n <- length(s)
  # special case of one genotype only
  if (n == 1) {
    return(list(list(s)))
  }

  masks <- 2^(1:n - 1) # bitmasks 00..001, 00..010, etc.
  c(
    list(list(s)), # list that contains the original vector as its only element
    lapply(seq(1, 2^n - 2, 2), function(u) { # u = odd numbers from 1 to 2^n-3
      # 2^n - 1 is excluded as its binary representation is all 1s
      # bits of u that are 1 go into first vector, rest go to second vector
      list(
        s[bitwAnd(u, masks) != 0],
        s[bitwAnd(u, masks) == 0]
      )
    })
  )
}

#' Pre-process data to remove repeats and \code{NA}s
#'
#' Removes repeat alleles and all \code{NA}s from allelic vectors with non-\code{NA} values.
#' Removes repeat \code{NA}s from allelic vectors with only \code{NA} values.
#'
#' @param y Observed data in the form of a list of lists. The outer list is a
#'   list of episodes in chronological order. The inner list is a list of named
#'   markers per episode. Episode names can be specified, but they are not used.
#'   Markers must be named. Each episode must list the same markers. If not all
#'   markers are typed per episode, data on untyped markers can be encoded as
#'   missing (see below). For each marker, one must specify an allelic vector: a
#'   set of distinct alleles detected at that marker. \code{NA}s encode missing
#'   per-marker data, i.e., when no alleles are observed for a given marker.
#'
#' @examples
#'
#' y <- list(list(m1 = c("A", "A", NA, "B"), m2 = c("A"), m3 = c("C")),
#'           list(m1 = c(NA, NA), m2 = c("B", "C"), m3 = c("A", "B", "C")))
#'
#' prep_data(y)
#'
#' @noRd
prep_data <- function(y) {
  warned_rep <- F
  warned_na <- F
  for(episode in 1:length(y)) {
    for(m in names(y[[episode]])) {
      alleles <- y[[episode]][[m]]
      # Collapse repeated alleles to single occurrence
      if(anyDuplicated(alleles)) {
        if(!warned_rep) {
          warning("Ignoring allele repeats at markers with observed data (or NA repeats at markers with missing data).")
          warned_rep <- T
        }
        alleles <- unique(alleles)
      }
      # Check NA is not mixed with named alleles for one marker + episode
      if(any(!is.na(alleles)) & any(is.na(alleles))) {
        if(!warned_na) {
          warning("Ignoring NAs among alleles detected at markers with observed data.")
          warned_na <- T
        }
        alleles <- alleles[!is.na(alleles)]
      }
      y[[episode]][[m]] <- alleles
    }
  }
  return(y)
}

#' Message progress bar
#'
#' @description
#' A simple progress bar to use in R packages where messages are preferred to
#' console output. See https://gist.github.com/MansMeg/1ec56b54e1d9d238b4fd.
#'
#' @field iter Total number of iterations
#' @field i Current iteration
#' @field width Width of the R console
#' @field width_bar Width of the progress bar
#' @field progress The number of character printed (continous)
#' @field progress_step Addition to progress per iteration
#'
#' @author Mans Magnusson (MansMeg @ github)
#' @importFrom methods new
#' @noRd
msg_progress_bar <-
  methods::setRefClass(
    Class = "msg_progress_bar",
    fields = list(iter = "numeric",
                  i = "numeric",
                  progress = "numeric",
                  progress_step = "numeric",
                  width = "numeric",
                  width_bar = "numeric"),

    methods = list(
      initialize = function(iter){
        'Initialize a messagebar object'
        .self$width <- getOption("width")
        .self$iter <- iter
        .self$i <- 0
        .self$progress <- 0
        white_part <- paste(rep(" ", (.self$width - 11) %/% 4), collapse="")
        init_length <- .self$width - ((.self$width - 11) %/% 4) * 4 - 11
        white_init <- paste(rep(" ", init_length), collapse="")
        .self$width_bar <- .self$width - init_length - 2 + 0.1
        .self$progress_step <- .self$width_bar / .self$iter
        #message(paste(white_init, "|", white_part, "25%", white_part, "50%", white_part, "75%", white_part, "|","\n", white_init, "|", sep=""), appendLF = FALSE)
      },

      increment = function(){
        'A messagebar object.'
        if(.self$i > .self$iter) return(invisible(NULL))
        new_progress <- .self$progress + .self$progress_step
        diff_in_char <- floor(new_progress) - floor(.self$progress)
        if(diff_in_char > 0) {
          message(paste(rep("=", diff_in_char),collapse=""), appendLF = FALSE)
        }

        .self$progress <- new_progress
        .self$i <- .self$i + 1
        if(.self$i == .self$iter) message("|\n", appendLF = FALSE)

      }
    )
  )
