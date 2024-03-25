library(Pv3Rs)
library(MCMCpack) # for rdirichlet

p <- runif(1, min=0, max=0.1) # rare allele that suggests initial genotypes are siblings
q <- runif(1, min=0, max=0.001) # very rare allele that suggests all genotypes are siblings

y <- list(
  init=list(m1=c("A","B"), m2=c("P"), m3=c("Q")),
  recur=list(m1=c("C"), m2=c("notP"), m3=c("Q"))
)

fs <- list(
  #note: m1 allele freqs do not affect the posterior
  m1=setNames(rdirichlet(1, rep(1,3))[1,], c("A","B","C")),
  m2=c(P=p, notP=1-p),
  m3=c(Q=q, notQ=1-q)
)

post <- compute_posterior(y, fs)

# posterior odds of relapse to reinfection
print(post$marg[,"L"] / post$marg[,"I"])

# for all values of p and q, these odds are bounded between 2/9 (as p->0)
# and 4/9 (as p->1, q->0)

# desired behaviour: a rare Q allele ought to suggest that all genotypes are
#                    siblings, leading to a posterior that favours relapse

# observed outcome: even in the upper bound, reinfection is still more likely

# extension: adding more markers like m2 can decrease the lower bound to 0,
#            but adding more markers like m3 does not increase the upper bound

# odds can be analytically verified
odds_func <- function(p, q) 2/9 * (1 + p*(1+q)/4 / ((1+p)*(1+q)/8 + p*q))
print(odds_func(p, q))

# note that a small p has a stronger effect in lowering the odds than a small q
ps <- 0:0.001:1
qs <- 0:0.001:1
odds <- outer(ps, qs, odds_func)
filled.contour(odds, color = function(x) heat.colors(x),
               plot.title=title(main="Post. odds of relapse:reinfection",
                                xlab="p", ylab="q")
               )
