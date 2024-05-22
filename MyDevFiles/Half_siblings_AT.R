################################################################################
# This script contains
#
# 1) The trajectory of the posterior relapse probability upon the addition of an
# inter- versus intra-match upon a base of all different plus all match (3
# markers max)
#
# 2) Posterior relapse probabilities as a function of the intra- versus
# inter-match allele frequencies for an example with one intra and two
# inter-matches, as well as one all different plus one all match (5 markers)
#
# 3) Posterior relapse probabilities for examples with all ways for children
# alleles to be drawn from parental alleles given different marker cardinalities
# assuming equifrequent alleles
################################################################################
rm(list = ls())
library(Pv3Rs)
library(MCMCpack) # for rdirichlet
set.seed(2) # for reproducibility

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

  Like_L <- (1/9) * (2*(rg1 +rg2) + rg3.1 + rg4.1 + rg3.2 + rg4.2)
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

abline(v = p, h = q); postL_function(fs, p, q)



#===============================================================================
# Posterior relapse probabilities for data in which of all ways to draw alleles
# for a trio of half siblings are represented assuming equifrequent alleles.
# Marker cardinality impacts the All different to All match ratio but not the
# intra-to-inter match ratio. Half siblings cannot be all different if markers
# are biallelic.
#===============================================================================
allele_counts <- c(2,3,4,5) # Cardinality of marker
exp_locus_type_props <- c("All diff." = NA,
                          "All match" = NA,
                          "Intra-match" = NA,
                          "Inter-match" = NA)

# Does cardinality matter to proportions
sapply(allele_counts, function(allele_count){

  # Generate alleles
  halfsib_alleles <- enumerate_halfsib_alleles(allele_count)
  row.names(halfsib_alleles) <- paste0("m", 1:nrow(halfsib_alleles))

  # Format into a list to pass to locus_type_summary
  halfsib_y <- list(init = as.list(halfsib_alleles[,1]),
                    recur = apply(halfsib_alleles[,2:3], 1, unique))

  # Generate locus types
  halfsib_locus_types <- sapply(1:nrow(halfsib_alleles), locus_type_summary, y = halfsib_y)

  # Compute locus type proportions
  x <- table(halfsib_locus_types)/nrow(halfsib_alleles)
  exp_locus_type_props[names(x)] <- x

  # Draw frequencies
  fs <- sapply(row.names(halfsib_alleles), function(m){
    setNames(rep(1/allele_count, allele_count), LETTERS[1:allele_count])
  }, simplify = F)

  # Compute posterior state probabilities
  post <- compute_posterior(halfsib_y, fs)

  # Unpack posterior relapse probability
  return(c(exp_locus_type_props,
           marker_count = nrow(halfsib_alleles),
           posterior_relapse_pr = post$marg[,"L"]))
})




