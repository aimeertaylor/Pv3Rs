################################################################################
# When alleles are equifrequent, why does the Pr(L | D) more frequently equal
# one as the number of markers increases?
#
# If there are three distinct parental alleles for children to draw from, there
# are 27 ways for children alleles to be drawn (see below), among which 6/27
# result in intra episode matches (suggestive of reinfection), and 12/27 result
# in inter episode matches (suggestive of relapse).
#
# If there are two distinct parental alleles for children to draw from, there
# are 8 ways for children alleles to be drawn (see below), among which 2/8
# result in intra episode matches (suggestive of reinfection), and 4/8 result in
# inter episode matches (suggestive of relapse).
#
# Despite more ways to match inter versus intra alleles, the posterior
# probability of relapse versus reinfection when all ways are represented
# depends on the allele frequencies, which in turn depend heavily on the
# concentration parameter when marker cardinality is low.
################################################################################
library(Pv3Rs)
library(MCMCpack)
set.seed(2)

# Function to draw allele frequencies
fs_function <- function(c = NULL, allele_count) {
  if (is.na(c)) { # If no concentration parameter
    p <- rep(1/allele_count, allele_count)
  } else {
    p <- rdirichlet(1, alpha = rep(c, allele_count))
  }
  return(p)
}

# Example data assuming there are two parental alleles to choose from:
y2 <- list(
  init=list(
    m1="A",  # All match
    m2="B",  # All match
    #
    m3="A",  # Intra match
    m4="B",  # Intra match
    #
    m5="A",  # Inter match
    m6="A",  # Inter match
    m7="B",  # Inter match
    m8="B"), # Inter match
  #
  recur=list(
    m1="A",
    m2="B",
    #
    m3="B",
    m4="A",
    #
    m5=c("A","B"),
    m6=c("B","A"),
    m7=c("A","B"),
    m8=c("B","A")
  )
)

# Example assuming there are three parental alleles to choose from:
y3 <- list(
  init=list(
    m1="A", # All different
    m2="B", # All different
    m3="C", # All different
    m4="A", # All different
    m5="B", # All different
    m6="C", # All different
    #
    m7="A", # All match
    m8="B", # All match
    m9="C", # All match
    #
    m10="A", # Intra match
    m11="A", # Intra match
    m12="B", # Intra match
    m13="B", # Intra match
    m14="C", # Intra match
    m15="C", # Intra match
    #
    m16="A", # Inter match
    m17="A", # Inter match
    m18="B", # Inter match
    m19="B", # Inter match
    m20="C", # Inter match
    m21="C", # Inter match
    m22="A", # Inter match
    m23="A", # Inter match
    m24="B", # Inter match
    m25="B", # Inter match
    m26="C", # Inter match
    m27="C"), # Inter match
  #
  recur=list(
    m1=c("B","C"),
    m2=c("A","C"),
    m3=c("A","B"),
    m4=c("B","C"),
    m5=c("A","C"),
    m6=c("A","B"),
    #
    m7="A",
    m8="B",
    m9="C",
    #
    m10="B",
    m11="C",
    m12="A",
    m13="C",
    m14="A",
    m15="B",
    #
    m16=c("A","B"),
    m17=c("A","C"),
    m18=c("B","A"),
    m19=c("B","C"),
    m20=c("C","A"),
    m21=c("C","B"),
    m22=c("A","B"),
    m23=c("A","C"),
    m24=c("B","A"),
    m25=c("B","C"),
    m26=c("C","A"),
    m27=c("C","B")
  )
)

plot_data(ys = list(pid1 = y3), fs = NULL, marker_annotate = F)
allele_counts <- c(2,3,4,5) # Cardinality of marker
conc_params <- c(0.1, 1, 1000) # Concentration parameter of allele frequency distribution
fs_store <- list()

# Generate frequencies
for(allele_count in allele_counts) {
  for(conc_param in conc_params) {

    # Draw frequencies
    fs <- sapply(paste0("m", 1:27),
                 function(x) setNames(fs_function(conc_param, allele_count), LETTERS[1:allele_count]),
                 simplify = F)

    # Store frequencies
    fs_store[[as.character(allele_count)]][[as.character(conc_param)]] <- fs
  }
}

# Compute posteriors
Two_alleles <- sapply(as.character(allele_counts), function(allele_count) {
  sapply(as.character(conc_params), function(conc_param) {
    fs <- fs_store[[allele_count]][[conc_param]] # Unpack frequencies
    post2 <- compute_posterior(y2, fs[1:8]) # Compute posterior state probabilities
    return(post2$marg[,"L"])  # Unpack posterior relapse probability
  })
})

# Compute posteriors
Three_alleles <- sapply(as.character(allele_counts)[-1], function(allele_count) {
  sapply(as.character(conc_params), function(conc_param) {
    fs <- fs_store[[allele_count]][[conc_param]] # Unpack frequencies
    post3 <- compute_posterior(y3, fs) # Compute posterior state probabilities
    return(post3$marg[,"L"]) # Unpack posterior relapse probability
  })
})


Two_alleles
Three_alleles

