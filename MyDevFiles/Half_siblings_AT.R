################################################################################
# This script contains
#
# 1) The trajectory of the posterior relapse probability upon the addition of an
# inter- versus intra-match upon a base of all different plus all match
#
# 2) Posterior relapse probabilities as a function of the intra- versus
# inter-match allele frequencies for an example with one intra and two
# inter-matches
#
# 3) Posterior relapse probabilities for examples with all 27 and 8 ways for
# children alleles to be drawn give 3 and 2 parental alleles to draw from
################################################################################

library(Pv3Rs)
library(MCMCpack) # for rdirichlet
set.seed(2) # for rdirichlet

# ==============================================================================
# Example showing the trajectory of the posterior relapse probability upon the
# addition of an inter- versus intra-match (adapted from Yong See's example),
# where m1 is tri-allelic, and m2 and m3 are bi-allelic
# ==============================================================================

# Example with all different
y_base1 <- list(
  init=list(m1=c("A")),
  recur=list(m1=c("B","C"))
)

# Example with all different, and all same
y_base2 <- list(
  init=list(m1=c("A"), m2=c("D")),
  recur=list(m1=c("B","C"), m2=c("D"))
)

# Example with all different, all same, and intra-match
y_intra <- list(
  init=list(m1=c("A"), m2=c("D"), m3=c("notP")),
  recur=list(m1=c("B","C"), m2=c("D"), m3=c("P"))
)

# Example with all different, all same, and inter-match
y_inter <- list(
  init=list(m1=c("A"), m2=c("D"), m3=c("P")),
  recur=list(m1=c("B","C"), m2=c("D"), m3=c("P","notP"))
)

d <- runif(1, min=0, max=0.1) # sample rare allele that suggests all genotypes are siblings
p <- runif(1, min=0, max=0.1) # sample rare allele that suggests two genotypes are siblings

fs <- list(
  # note: m1 allele freqs do not affect the posterior
  m1=setNames(rdirichlet(1, rep(1,3))[1,], c("A","B","C")),
  m2=c(D=d, notD=1-d),
  m3=c(P=p, notP=1-p)
)

# Compute posterior relapse probability and odds
post_base2 <- function(fs) {

  d <- fs$m2["D"]

  rg1 <- d
  rg2 <- (d+1)/4
  rg4 <- rg3 <- rg2

  Like_L <- (1/9) * (rg1 + rg2 + rg3 + rg4)
  Like_I <- (1/2) * (rg1 + rg2)

  post_L = Like_L / (Like_L + Like_I)
  post_odds = (2/9) * ((7*d + 3)/(5*d + 1))
  return(c(post_L = post_L, post_odds = post_odds))
}

# Compute posterior relapse probability and odds
post_inter <- function(fs) {

  d <- fs$m2["D"]
  p <- fs$m3["P"]

  rg1 <- d*p
  rg2 <- (d+1)*p/8
  rg3 <- (d+1)*(p+1)/8
  rg4 <- rg2

  Like_L <- (1/9) * (rg1 + rg2 + rg3 + rg4)
  Like_I <- (1/2) * (rg1 + rg2)

  post_L = Like_L / (Like_L + Like_I)
  post_odds = 2/9 * (1 + ((2*p + 1)*(d+1)/8) / (p*d + p*(d+1)/8))
  return(c(post_L = post_L, post_odds = post_odds))
}


# Compute posterior relapse probability and odds
post_intra <- function(fs) {

  d <- fs$m2["D"]
  p <- fs$m3["P"]

  rg1 <- p*d
  rg2 <- (p+1)*(d+1)/8
  rg3 <- p*(d+1)/8
  rg4 <- rg3

  Like_L <- (1/9) * (rg1 + rg2 + rg3 + rg4)
  Like_I <- (1/2) * (rg1 + rg2)

  post_L = Like_L / (Like_L + Like_I)
  post_odds = 2/9 * (1 + (p*(d+1)/4) / (p*d + (p+1)*(d+1)/8))
  return(c(post_L = post_L, post_odds = post_odds))
}

# Compute posterior state probabilities under Pv3Rs model
p_base1 <- compute_posterior(y_base1, fs[1])
p_base2 <- compute_posterior(y_base2, fs[1:2])
p_inter <- compute_posterior(y_inter, fs)
p_intra <- compute_posterior(y_intra, fs)

# Check results match
c(p_base1$marg[,"L"], p_base1$marg[,"L"] / p_base1$marg[,"I"])
(2/9)*(5/3) # Odds for one all different marker

c(p_base2$marg[,"L"], p_base2$marg[,"L"] / p_base2$marg[,"I"])
post_base2(fs)

c(p_inter$marg[,"L"], p_inter$marg[,"L"] / p_inter$marg[,"I"])
post_inter(fs)

c(p_intra$marg[,"L"], p_intra$marg[,"L"] / p_intra$marg[,"I"])
post_intra(fs)

# Plot trajectory
plot(NULL, xlim = c(1,3), ylim = c(0,1), bty = "n", xaxt = "n",
     las = 1, ylab = "Posterior relapse probability", xlab = "Marker count")
axis(side = 1, at = c(1,2,3))
abline(h = 0.5, lty = "dashed")
arrows(x0 = c(1,2,2), x1 = c(2,3,3), length = 0.075,
       y0 = c(p_base1$marg[,"L"], p_base2$marg[,"L"], p_base2$marg[,"L"]),
       y1 = c(p_base2$marg[,"L"], p_intra$marg[,"L"], p_inter$marg[,"L"]))
text(x = c(1,2,3,3),
     y = c(p_base1$marg[,"L"], p_base2$marg[,"L"], p_intra$marg[,"L"], p_inter$marg[,"L"]),
     labels = c("\nAll different", "\nAll match", "\nIntra-match", "Inter-match"),
     pos = c(4, 1, 2, 2))


# ==============================================================================
# Posterior relapse probabilities as a function of the intra- versus inter-match
# allele frequencies for an example with one intra and two inter-matches
# (requires phasing!), where m1 is tri-allelic, and m2 to m5 are bi-allelic.
# ==============================================================================
q <- 0.4 # inter-match allele frequency

# Data to pass to compute_posterior
y <- list(
  init=list(m1="A", m2 = "D", m3="notP", m4="Q", m5="Q"),
  recur=list(m1=c("B","C"), m2 = "D", m3="P", m4=c("Q", "notQ"), m5=c("Q", "notQ"))
)

fs <- list(
  # note: m1 allele freqs do not affect the posterior
  m1=setNames(rdirichlet(1, rep(1,3))[1,], c("A","B","C")),
  m2=c(D=d, notD=1-d),
  m3=c(P=p, notP=1-p),
  m4=c(Q=q, notQ=1-q),
  m5=c(Q=q, notQ=1-q)
)

# Plot the data
plot_data(ys = list(pid1 = y), fs = fs)

# Compute posterior relapse probability, summing over phase
postL_function <- function(fs, p, q, prob_only = FALSE, odds_only = FALSE) {

  rg1 <-d*p*q*(1-q)*q^2
  rg2 <- (1/32)*(d+1)*(p+1)*q*(1-q)*q^2

  rg3.1 <- (1/32)*(d+1)*p*(q+1)*((1-q)*q^2 + (1-q)*q)
  rg4.1 <- (1/32)*(d+1)*p*q*(1-q)*q^2

  rg3.2 <- (1/32)*(d+1)*p*(q+1)*(1-q)*q^2
  rg4.2 <- (1/32)*(d+1)*p*(q+1)*(1-q)*q^2

  Like_L <- (1/9) * (2*(rg1 +rg2) +rg3.1 +rg4.1 +rg3.2 +rg4.2)
  Like_I <- (1/2) * (2*(rg1 +rg2))

  post_L = Like_L / (Like_L + Like_I)
  post_odds = Like_L / Like_I

  if(odds_only & prob_only) stop("Choose probs, odds, or neither")
  if(odds_only) return(post_odds)
  if(prob_only) return(post_L)
  if(!prob_only & !odds_only) {
    return(c(post_L = post_L, post_odds = post_odds))
  }
}

# Check results agree with those computed under the Pv3Rs model
post <- compute_posterior(y, fs) #Compute posterior state probabilities
c(post$marg[,"L"], post$marg[,"L"] / post$marg[,"I"])
postL_function(fs, p, q)

# Plot posterior relapse probability as a function of the intra- versus
# inter-match allele frequency
ps <- seq(0,1,0.005)
qs <- seq(0.001,1,0.005)

fields::image.plot(outer(ps, qs, postL_function, fs = fs, prob_only = TRUE),
                   ylab = "Inter-match allele frequency",
                   xlab = "Intra-match allele frequency",
                   main = "Posterior relapse probability",
                   col = RColorBrewer::brewer.pal(n = 10, "RdBu"),
                   breaks = seq(0,1,length.out = 11))

abline(v = p, h = q); postL_function(fs, p, q) # Something odd about colours



#===============================================================================
# Posterior relapse probabilities for examples with all 27 and 8 ways for
# children alleles to be drawn give 3 and 2 parental alleles to draw from
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
#===============================================================================
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



