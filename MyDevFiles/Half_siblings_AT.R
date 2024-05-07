library(Pv3Rs)
library(MCMCpack) # for rdirichlet
set.seed(2)

p <- runif(1, min=0, max=0.1) # sample rare allele that suggests intra-inf. genotypes are siblings
q <- p #runif(1, min=0, max=0.1) # sample very rare allele that suggests inter-inf. genotypes are siblings

# Data to pass to compute_posterior
y <- list(
  init=list(m1="A", m2="notP", m3="Q", m4="Q", m5 = "Q"),
  recur=list(m1=c("B","C"), m2="P", m3=c("Q", "notQ"), m4=c("Q", "notQ"), m5 = "Q")
)

# y <- list(
#   init=list(m1="A", m2="notP", m3="Q"),
#   recur=list(m1=c("B","C"), m2="P", m3=c("Q", "notQ"))
# )

fs <- list(
  # note: m1 allele freqs do not affect the posterior
  m1=setNames(rdirichlet(1, rep(1,3))[1,], c("A","B","C")),
  m2=c(P=p, notP=1-p),
  m3=c(Q=q, notQ=1-q),
  m4=c(Q=q, notQ=1-q),
  m5=c(Q=q, notQ=1-q)
)

# Posterior relapse probability function assuming uniform priors on recurrent states
# postL_function <- function(fs, p, q) {
#   Like_L <- (1/9) * (p*q + q*(p+1)/8 + p*(q+1)/8 + p*q/8)
#   Like_I <- (1/2) * (p*q + q*(p+1)/8)
#   post_L = Like_L / (Like_L + Like_I)
#   return(post_L)
# }

postL_function <- function(fs, p, q) {
  Like_L <- (1/9) * (p*q + q*(p+1)/16 + p*(q+1)/8)
  Like_I <- (1/2) * (p*q + q*(p+1)/16)
  post_L = Like_L / (Like_L + Like_I)
  return(post_L)
}

# Posterior relapse:reinfection odds function assuming uniform priors on recurrent states
# odds_func <- function(p, q) 2/9 * (1 + ((p*(q+1) + p*q)/8) / (p*q + q*(p+1)/8))
odds_func <- function(p, q) 2/9 * (1 + (p*(q+1)/8) / (p*q + q*(p+1)/16))

# Plot the data
plot_data(ys = list(pid1 = y), fs = fs)


# Compute posterior state probabilities
post <- compute_posterior(y, fs)
post$marg

postL_function(fs, p, q)*1/2

# Compute posterior odds of relapse to reinfection
writeLines(sprintf("Odds of relapse to reinfection: %s",
                   round(post$marg[,"L"] / post$marg[,"I"], 4)))
print(odds_func(p, q))

# Plots
ps <- seq(0,0.999,0.001)
qs <- seq(0.001,1,0.001)

fields::image.plot(outer(ps, qs, postL_function, fs = fs),
                   ylab = "q",
                   xlab = "p",
                   main = "Posterior relapse probability",
                   col = RColorBrewer::brewer.pal(n = 11, "RdBu"),
                   breaks = seq(0,1,length.out = 12))

fields::image.plot(outer(ps, qs, odds_func),
                   ylab = "q",
                   xlab = "p",
                   main = "Posterior odds of relapse:reinfection",
                   col = RColorBrewer::brewer.pal(n = 11, "RdBu"),
                   breaks = seq(0,1,length.out = 12))



