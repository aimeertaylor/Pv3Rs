##########
# Simple illustration of what happens when reinfection is excluded a priori
#
# There are one recurrent episode, and the same allele is observed at enrollment
# and recurrence. By reinfection is excluded from the prior, the posterior
# probability of recrudescence is less sensitive to the frequency of the allele
# observed.
##########

library(Pv3Rs)

f <- runif(1)
fs <- list(m1=c("A" = f, "B" = 1-f))
y <- list(enroll = list(m1="A"), recur = list(m1="A"))

# prior includes reinfection: check that numerical and analytical results agree for random f
suppressMessages(compute_posterior(y, fs))$marg # # numerical marginal posterior
c(2, 1+f, 2*f) / (3+3*f) # analytical formula for marginal posterior, see preprint 4.2


# prior excludes reinfection: check that numerical and analytical results agree for random f
no_reinf_prior <- matrix(c(0.5,0.5,0), nrow=1, dimnames=list(NULL, c("C","L","I")))
suppressMessages(compute_posterior(y, fs, prior=no_reinf_prior))$marg # marginal posterior
c(2, 1+f) / (3+f) # rescale previous analytical formula to exclude


# plot comparison
par(mfrow=c(1,2))
fgrid <- seq(0, 1, 0.01)

plot(fgrid, 2/(3+fgrid), type="l", col="red",
     ylim=c(0, 1),
     main="Posterior probability of recrudescence", xlab="Allele frequency", ylab="")
lines(fgrid, 2/(3+3*fgrid))
legend("topright", legend = c("include reinfection", "exclude reinfection"),
       col=c("black","red"), lty=c(1,1))

plot(fgrid, (1+fgrid)/(3+fgrid), type="l", col="red",
     ylim=c(0, 1),
     main="Posterior probability of relapse", xlab="Allele frequency", ylab="")
lines(fgrid, (1+fgrid)/(3+3*fgrid))
legend("topright", legend = c("include reinfection", "exclude reinfection"),
       col=c("black","red"), lty=c(1,1))
