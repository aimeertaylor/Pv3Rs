################################################################################
# Code to prepare fs_VHX_BPD dataset following best practice outlined in 7.1.1
# of https://r-pkgs.org/data.html. This script only features in the source
# version of the package (it is listed under .Rbuildignore). Path to real data:
# https://github.com/jwatowatson/RecurrentVivax/blob/master/RData/GeneticModel/
#################################################################################

rm(list = ls())
load('MS_data_PooledAnalysis.RData')
MSs_all <- tail(names(MS_pooled), 9) # marker names
D_weight_Prior <- 1 # uniform Dirichlet prior for allele frequencies
ind_enrol = which(MS_pooled$Episode == 1) # Indices of enrolment episodes

# Compute the posterior Dirichlet parameter vector
Alpha_Posteriors <- apply(MS_pooled[,MSs_all], 2, function(x, ind_enrol){

  # number of possible alleles for this marker
  xmax <-  max(x, na.rm=T)

  # prior concentration parameters
  param_vector <- array(D_weight_Prior, dim = xmax, dimnames = list(1:xmax))

  # observed data summarised as counts
  obs_counts <- table(x[ind_enrol])

  # posterior concentration parameters
  param_vector[names(obs_counts)] <-  param_vector[names(obs_counts)] + obs_counts

  return(param_vector)

}, ind_enrol)

# Normalise
fs_VHX_BPD <- sapply(Alpha_Posteriors, function(x){x/sum(x)})

# Save fs_VHX_BPD as exported data
usethis::use_data(fs_VHX_BPD, overwrite = TRUE)
