# To integrate into tests

# Test that fails when markers different
fs = list(m1 = c("1" = 1), m2 = c("1" = 1))
y <- list(epi1 = list(m1 = "1", m2 = NA),
          epi2 = list(m1 = NA))
compute_posterior(y, fs)

# Test that no diversity returns the prior
fs = list(m1 = c("1" = 1))
y <- list(list(m1 = "1"), list(m1 = "1"))
compute_posterior(y, fs)


# Test that no common allele returns the prior (only if MOIs not supplied!)
fs = list(m1 = c("1" = 0.6, "other" = 0.4),
          m2 = c("1" = 0.6, "other" = 0.4),
          m3 = c("1" = 0.6, "other" = 0.4))
y <- list(list(m1 = "1", m2 = NA, m3 = "1"),
          list(m1 = "1", m2 = NA, m3 = "1"),
          list(m1 = NA, m2 = "1", m3 = NA))
compute_posterior(y, fs, MOIs = c(2,3,1))


# Ensure for all episodes there are inter-episode data on at least one common marker
y_by_episode <- plyr::dlply(y.df, 'Episode') # Group data by episode

# Check there's at least one common marker with inter-episode data per episode
some_common_per_epi <- sapply(y, function(y_epi) {
  is.na(y_epi)

  other_epis <- setdiff(names(y), epi) # Other episode names

  some_common <- any(sapply(MSs_all, function(ms){

    # Are the alleles at marker ms in epi all NAs?
    per_marker_epi_data_all_na <- all(is.na(y_by_episode[[epi]][ms]))

    # Are the alleles at marker ms in the other episodes all NAs?
    per_marker_other_epis_data_all_na <- all(sapply(other_epis, function(other_epi) {
      all(is.na(y_by_episode[[other_epi]][ms]))}))

    return(!per_marker_epi_data_all_na & !per_marker_other_epis_data_all_na)
  }))

  return(some_common)
})
