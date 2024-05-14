library(Pv3Rs)
library(MCMCpack) # for rdirichlet
set.seed(2)

r <- 0.2
p <- r #runif(1, min=0, max=0.1) # sample rare allele that suggests intra-inf. genotypes are siblings
q <- r #runif(1, min=0, max=0.1) # sample very rare allele that suggests inter-inf. genotypes are siblings


# Data to pass to compute_posterior
y <- list(
  init=list(m1="A", m2="notP", m3="Q", m4="Q", m5="R"),
  recur=list(m1=c("B","C"), m2="P", m3=c("Q", "notQ"), m4=c("Q", "notQ"),  m5="R")
)

# Plot the data
plot_data(ys = list(pid1 = y), fs = fs)


fs <- list(
  # note: m1 allele freqs do not affect the posterior
  m1=setNames(rdirichlet(1, rep(1,3))[1,], c("A","B","C")),
  m2=c(P=p, notP=1-p),
  m3=c(Q=q, notQ=1-q),
  m4=c(Q=q, notQ=1-q),
  m5=c(R=r, notR=1-r)
)


postL_function <- function(fs, p, q, prob_only = FALSE, odds_only = FALSE) {

  rg1 <- r*p*q*(1-q)*q^2
  rg2 <- (1/32)*(r+1)*(p+1)*q*(1-q)*q^2

  rg3.1 <- (1/32)*(r+1)*p*(q+1)*((1-q)*q^2 + (1-q)*q)
  rg4.1 <- (1/32)*(r+1)*p*q*(1-q)*q^2

  rg3.2 <- (1/32)*(r+1)*p*(q+1)*(1-q)*q^2
  rg4.2 <- (1/32)*(r+1)*p*(q+1)*(1-q)*q^2

  Like_L <- (1/9) * (2*(rg1 + rg2) + rg3.1 + rg4.1 + rg3.2 + rg4.2)
  Like_I <- (1/2) * (2*(rg1 + rg2))

  post_L = Like_L / (Like_L + Like_I)
  post_odds = Like_L / Like_I

  if(odds_only & prob_only) stop("Choose probs, odds, or neither")
  if(odds_only) return(post_odds)
  if(prob_only) return(post_L)
  if(!prob_only & !odds_only) {
    return(c(post_L = post_L, post_odds = post_odds))
  }
}

# Compare results
post <- compute_posterior(y, fs) #Compute posterior state probabilities
c(post$marg[,"L"], post$marg[,"L"] / post$marg[,"I"])
postL_function(fs, p, q)


# Plots
ps <- seq(0,0.999,0.005)
qs <- seq(0.001,1,0.005)

fields::image.plot(outer(ps, qs, postL_function, fs = fs, prob_only = TRUE),
                   ylab = "q",
                   xlab = "p",
                   main = "Posterior relapse probability",
                   col = RColorBrewer::brewer.pal(n = 10, "RdBu"),
                   breaks = seq(0,1,length.out = 11))


fields::image.plot(outer(ps, qs, postL_function, fs = fs, odds_only = TRUE),
                   ylab = "q",
                   xlab = "p",
                   main = "Posterior odds of relapse:reinfection",
                   col = RColorBrewer::brewer.pal(n = 10, "RdBu"),
                   breaks = seq(0,1,length.out = 11))



