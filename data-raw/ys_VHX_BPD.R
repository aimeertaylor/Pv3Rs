################################################################################
# Code to prepare ys_VHX_BPD dataset following best practice outlined in 7.1.1
# of https://r-pkgs.org/data.html. This script only features in the source
# version of the package (it is listed under .Rbuildignore). Path to real data:
# https://github.com/jwatowatson/RecurrentVivax/blob/master/RData/GeneticModel/
#################################################################################

rm(list = ls())
load('MS_data_PooledAnalysis.RData')
MSs_all <- tail(names(MS_pooled), 9) # marker names
y.dfs <- plyr::dlply(MS_pooled, 'ID') # Group genetic data by patient ID

# ==============================================================================
# Reformat data
# ==============================================================================
ys_VHX_BPD <- lapply(names(y.dfs), function(pid) {

  # Extract data for a given patient IDs
  y.df <- y.dfs[[pid]]

  # Extract markers for which there is at least one non-NA
  ms <- MSs_all[apply(!is.na(y.df[MSs_all]), 2, any)]

  # Group data by episode
  y.by.episode <- plyr::dlply(y.df, 'Episode_Identifier')

  # Transform data frame format to format taken by 'compute_posterior'
  y <- lapply(y.by.episode, function(episode) { # For each episode
    setNames(lapply(ms, function(m) { # For each marker
      alleles <- episode[m][!is.na(episode[m])] # Extract non-NAs
      if (length(alleles) > 0) {
        return(alleles)
      } else {
        return(NA)
      }
    }), ms)
  })

  # Covert episode IDs into episode indices
  names(y) <- do.call(rbind, strsplit(names(y), split = "_"))[,3]

  # Check episodes ordered correctly and if not re-order and re-check
  if(!all(as.numeric(names(y)) == cummax(as.numeric(names(y))))){
    y <- y[order(as.numeric(names(y)))]
    if (!all(as.numeric(names(y)) == cummax(as.numeric(names(y))))) stop()
  }

  # Return patient data
  return(y)
})

names(ys_VHX_BPD) <- names(y.dfs)

# Save ys_VHX_BPD as exported data
usethis::use_data(ys_VHX_BPD, overwrite = TRUE)
