#' Example Plasmodium vivax data
#'
#' Previously-published microsatellite data on Plasmodium vivax parasites
#' extracted from patients enrolled in the Best Primaquine Dose (BPD) and Vivax
#' History (VHX) trials; see
#' <https://www.nature.com/articles/s41467-019-13412-x>.
#'
#' @format A list of 217 patients; for each patient, a list of one or more
#' episodes; for each episode, a list of three or more microsatellite markers;
#' for each marker, a vector of observed alleles (repeat lengths). For example:
#' \describe{
#'   \item{BPD_103}{Patient identifier: patient 103 in the BPD trial}
#'   \item{BPD_103_1}{Episode identifier: episode one of patient 103 in the BPD trial}
#'   \item{PV.3.27}{Marker identifier: Plasmodium vivax 3.27}
#'   \item{18}{Repeat length: 18}
#' }
#'
#' @source
#' * <https://github.com/jwatowatson/RecurrentVivax/blob/master/RData/GeneticModel/MS_data_PooledAnalysis.RData>
#' * <https://github.com/aimeertaylor/Pv3Rs/blob/main/data-raw/ys_VHX_BPD.R>
"ys_VHX_BPD"



#' Allele frequencies computed from example Plasmodium vivax data
#'
#' The posterior mean of a multinomial-Dirichlet model with uniform prior fit to
#' data on allele prevalence in initial episodes. Because the model is fit to
#' allele prevalence (observed) not allele frequency (would requires
#' integrating-out unknown multiplicities of infection) it is liable to
#' underestimate the frequencies of common alleles and overestimate those of
#' rare but detected alleles.
#'
#' @format
#' A list of nine markers; for each marker a named vector of allele frequencies that sum to one.
#'
#' @source
#' * <https://github.com/jwatowatson/RecurrentVivax/blob/master/RData/GeneticModel/MS_data_PooledAnalysis.RData>
#' * <https://github.com/aimeertaylor/Pv3Rs/blob/main/data-raw/fs_VHX_BPD.R>
"fs_VHX_BPD"
