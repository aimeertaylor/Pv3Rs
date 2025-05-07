################################################################################
# Generate the figure `cond_allele_freq.png` in `half_siblings.ltx`
# This aims to theoretically corroborate results pertaining to a data-rich
# limit for 3 half siblings (with parents AB, BC, CA).
################################################################################

library(MCMCpack)

# eq:cond_dens in half_siblings.ltx
myden <- function(x) (3*lambda+2)/(4*lambda+3) * (1+x) * dbeta(x, lambda+1,2*lambda+1)

png("../vignettes/figures/cond_allele_freq.png", width = 6, height = 4, units = "in", res = 75)

lambda <- 10
curve(myden, from=1e-4, to=1-1e-4, add=F, col=1,
      xlab="Allele frequency", ylab="", main="Conditional density given obs. case 3 or 4")
lambda <- 3
curve(myden, from=1e-4, to=1-1e-4, add=T, col=2)
lambda <- 1
curve(myden, from=1e-4, to=1-1e-4, add=T, col=3)
lambda <- 0.3
curve(myden, from=1e-4, to=1-1e-4, add=T, col=4)
lambda <- 0.1
curve(myden, from=1e-4, to=1-1e-4, add=T, col=5)
legend(0.76, 4.5, legend=c(
  expression(paste(lambda, " = 10   ")),
  expression(paste(lambda, " = 3   ")),
  expression(paste(lambda, " = 1   ")),
  expression(paste(lambda, " = 0.3   ")),
  expression(paste(lambda, " = 0.1   "))
), col=1:5, lty=rep(1,5), cex=0.8)

dev.off()
