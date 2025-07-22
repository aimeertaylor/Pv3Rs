#' Example *Plasmodium vivax* data
#'
#' Previously-published microsatellite data on *P. vivax* parasites
#' extracted from study participants enrolled in the Best Primaquine Dose (BPD) and Vivax
#' History (VHX) trials; see
#' [Taylor & Watson et al. 2019](https://www.nature.com/articles/s41467-019-13412-x)
#' for more details of the genetic data; for more details of the VHX and BPD trials, see
#' [Chu et al. 2018a](https://academic.oup.com/cid/article/67/10/1543/5034800) and
#' [Chu et al. 2018b](https://academic.oup.com/cid/article/68/8/1311/5079011),
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
#' * <https://zenodo.org/records/3368828>
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
#' * <https://zenodo.org/records/3368828>
#' * <https://github.com/aimeertaylor/Pv3Rs/blob/main/data-raw/fs_VHX_BPD.R>
"fs_VHX_BPD"
