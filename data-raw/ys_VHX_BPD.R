################################################################################
# Code to prepare ys_VHX_BPD dataset following best practice outlined in 7.1.1
# of https://r-pkgs.org/data.html. This script only features in the source
# version of the package (it is listed under .Rbuildignore)
#################################################################################

rm(list = ls())
load('../data/MS_data_PooledAnalysis.RData')
MSs_all <- tail(names(MS_pooled), 9) # marker names
y.dfs <- plyr::dlply(MS_pooled, 'ID') # Group genetic data by patient ID

# ==============================================================================
# Reformat data
# ==============================================================================
ys_VHX_BPD <- lapply(names(y.dfs), function(i) {

  # Extract data for a given patient
  y.df <- y.dfs[[i]]

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

  # Return patient data
  return(y)
})

names(ys_VHX_BPD) <- names(y.dfs)

# Save ys_VHX_BPD as exported data
usethis::use_data(ys_VHX_BPD, overwrite = TRUE)
