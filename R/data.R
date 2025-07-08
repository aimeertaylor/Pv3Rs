#' Example *Plasmodium vivax* data
#'
#' Previously-published microsatellite data on *P. vivax* parasites
#' extracted from study participants enrolled in the Best Primaquine Dose (BPD) and Vivax
#' History (VHX) trials; see
#' [Taylor & Watson et al. 2019](https://www.nature.com/articles/s41467-019-13412-x)
#' for more details of the genetic data; for more details of the trials, see
#' [Chu et al. 2018a](https://doi.org/10.1093/cid/ciy319) and
#' [Chu et al. 2018b](https://doi.org/10.1093/cid/ciy735).
#'
#' @format A list of 217 study participants; for each study participant, a list of one or more
#' episodes; for each episode, a list of three or more microsatellite markers;
#' for each marker, a vector of observed alleles (repeat lengths). For example:
#' \describe{
#'   \item{BPD_103}{Study participant identifier: study participant 103 in the BPD trial}
#'   \item{BPD_103_1}{Episode identifier: episode one of study participant 103 in the BPD trial}
#'   \item{PV.3.27}{Marker identifier: *P. vivax* 3.27}
#'   \item{18}{Allele identifier: 18 repeat lengths}
#' }
#'
#' @source
#' * <https://github.com/jwatowatson/RecurrentVivax/blob/master/RData/GeneticModel/MS_data_PooledAnalysis.RData>
#' * <https://github.com/aimeertaylor/Pv3Rs/blob/main/data-raw/ys_VHX_BPD.R>
"ys_VHX_BPD"



#' Allele frequencies computed using example *Plasmodium vivax* data
#'
#' The posterior mean of a multinomial-Dirichlet model with uniform prior fit to
#' data on allele prevalence in initial episodes of [ys_VHX_BPD]. Because the
#' model is fit to allele prevalence (observed) not allele frequency ( requires
#' integrating-out unknown multiplicities of infection) it is liable to
#' underestimate the frequencies of common alleles and overestimate those of
#' rare but detected alleles.
#'
#' @format
#' A list of nine markers; for each marker a named vector of allele frequencies
#' that sum to one.
#'
#' @source
#' * <https://github.com/jwatowatson/RecurrentVivax/blob/master/RData/GeneticModel/MS_data_PooledAnalysis.RData>
#' * <https://github.com/aimeertaylor/Pv3Rs/blob/main/data-raw/fs_VHX_BPD.R>
"fs_VHX_BPD"


#' Upper bounds on posterior probabilities
#'
#' Upper bounds on posterior probabilities of 1st-recurrence recrudescence and
#' reinfection (rows) for various MOI vectors (columns) that sum to at most
#' eight.
#' #' These upper bounds feature in
#' [vignette("understand-graph-prior")](../doc/understand-graph-prior.pdf) and
#' [vignette("intra-episode-siblings")](../doc/intra-episode-siblings.html).
#' For a single recurrence (MOI vectors of length two), upper bounds are
#' induced by the prior. For more than one recurrence (MOI vectors of length
#' three or more), upper bounds assume all but the first recurrence has data;
#' these bounds are only valid for monoclonal episodes; see XXX. Two upper
#' bounds for each recurrence state are reported: one where summation over
#' graphs includes graphs with cliques of three or more intra-episode siblings,
#' another where summation over graphs excludes graphs with cliques of three or
#' more intra-episode siblings. Since graphs with cliques of three or more
#' intra-episode siblings only apply when MOI vectors include MOIs of three or
#' more, these two bounds only differ for MOI vectors that include MOIs of three
#' or more; see **Examples**.
#'
#' @format
#' A matrix with 4 rows and 247 columns:
#' \describe{
#'   \item{C_with}{Row name; upper bound on the posterior probability of recrudescence
#'   when summation includes graphs with cliques of three or more intra-episode siblings.}
#'   \item{C_wout}{Row name; upper bound on the posterior probability of recrudescence
#'   when summation excludes graphs with cliques of three or more intra-episode siblings.}
#'   \item{I_with}{Row name; upper bound on the posterior probability of reinfection when
#'   summation includes graphs with cliques of three or more intra-episode siblings.}
#'   \item{I_wout}{Row name; upper bound on the posterior probability of reinfection when
#'   summation excludes graphs with cliques of three or more intra-episode siblings.}
#'   \item{11}{Column name; MOI vector (1, 1).}
#'   \item{12}{Column name; MOI vector (1, 2).}
#'   ...
#'   \item{11111111}{Column name; MOI vector (1, 1, 1, 1, 1, 1, 1, 1).}
#' }
#' @source <https://github.com/aimeertaylor/Pv3Rs/blob/main/data-raw/maxima.R>
#'
#' @examples
#' # Convert column names to MOI character vectors
#' MOIvec <- strsplit(colnames(maxima), split = "")
#'
#' # Get MOI vectors with MOIs that exceed two
#' Exceed2 <- sapply(MOIvec, function(x) any(as.numeric(x) > 2))
#'
#' # Compare reinfection upper bounds with and with summation for MOI vectors that exclude
#' # and include MOIs greater than two
#' all(maxima["I_with", !Exceed2] == maxima["I_wout", !Exceed2])
#' any(maxima["I_with", Exceed2] == maxima["I_wout", Exceed2])
"maxima"
