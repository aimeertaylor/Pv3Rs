library(Pv3Rs)
library(gtools) # for rdirichlet

set.seed(1)
causes <- c("C", "L", "I")
alleles <- c("A", "C", "G", "T")
n.as <- length(alleles)
alphas <- rep(1, n.as) # Dirichlet param. vector

# Example 1
y <- list(list(m1="A"), list(m1="T"))
n.m <- 1
fmat <- rdirichlet(n.m, alphas)
colnames(fmat) <- alleles
fs <- setNames(lapply(1:n.m, function(x) fmat[x,]), paste0("m", 1:n.m))
post <- compute_posterior(y, fs)
expect <- matrix(c(0, 1/3, 2/3), ncol=3, byrow=T)
colnames(expect) <- causes
testthat::expect_equal(post$marg, expect)


# Example 2
y <- list(list(m1="A"), list(m1="A"))
n.m <- 1
fmat <- rdirichlet(n.m, alphas)
colnames(fmat) <- alleles
fs <- setNames(lapply(1:n.m, function(x) fmat[x,]), paste0("m", 1:n.m))
post <- compute_posterior(y, fs)
l.ntor <- fs$m1['A']+1
i.ntor <- 2*fs$m1['A']
c.ntor <- 2
dtor <- sum(c(l.ntor, i.ntor, c.ntor))
expect <- matrix(c(c.ntor, l.ntor, i.ntor)/dtor,
                 ncol=3, byrow=T)
colnames(expect) <- causes
testthat::expect_equal(post$marg, expect)


# Example 3
y <- list(list(m1="A", m2="T", m3="G"),
          list(m1="T", m2="T", m3="G"))
n.m <- 3
fmat <- rdirichlet(n.m, alphas)
colnames(fmat) <- alleles
fs <- setNames(lapply(1:n.m, function(x) fmat[x,]), paste0("m", 1:n.m))
post <- compute_posterior(y, fs)
l.ntor <- (fs$m2['T']*fs$m3['G']+(fs$m2['T']+1)*(fs$m3['G']+1)/8)/3
i.ntor <- fs$m2['T']*fs$m3['G']
c.ntor <- 0
dtor <- sum(c(l.ntor, i.ntor, c.ntor))
expect <- matrix(c(c.ntor, l.ntor, i.ntor)/dtor,
                 ncol=3, byrow=T)
colnames(expect) <- causes
testthat::expect_equal(post$marg, expect)


# Example 4
y <- list(list(m1="A", m2="T", m3="G"),
          list(m1="A", m2="T", m3="G"))
n.m <- 3
fmat <- rdirichlet(n.m, alphas)
colnames(fmat) <- alleles
fs <- setNames(lapply(1:n.m, function(x) fmat[x,]), paste0("m", 1:n.m))

post <- compute_posterior(y, fs)
l.ntor <- (fs$m1['A']+1)*(fs$m2['T']+1)*(fs$m3['G']+1)/8 + fs$m1['A']*fs$m2['T']*fs$m3['G'] + 1
i.ntor <- 3*fs$m1['A']*fs$m2['T']*fs$m3['G']
c.ntor <- 3
dtor <- sum(c(l.ntor, i.ntor, c.ntor))
expect <- matrix(c(c.ntor, l.ntor, i.ntor)/dtor,
                 ncol=3, byrow=T)
colnames(expect) <- causes
testthat::expect_equal(post$marg, expect)


# Example 5
y <- list(list(m1="A"), list(m1="T"), list(m1="T"))
n.m <- 1
fmat <- rdirichlet(n.m, alphas)
colnames(fmat) <- alleles
fs <- setNames(lapply(1:n.m, function(x) fmat[x,]), paste0("m", 1:n.m))

post <- compute_posterior(y, fs)
l1.ntor <- 17/24*fs$m1['A']*fs$m1['T']^2 + 11/16*fs$m1['A']*fs$m1['T']
i1.ntor <- 7/5*fs$m1['A']*fs$m1['T']^2 + 13/10*fs$m1['A']*fs$m1['T']
c1.ntor <- 0
dtor1 <- sum(c(l1.ntor, i1.ntor, c1.ntor))
l2.ntor <- 73/120*fs$m1['A']*fs$m1['T']^2 + 39/80*fs$m1['A']*fs$m1['T']
i2.ntor <- 3/2*fs$m1['A']*fs$m1['T']^2
c2.ntor <- 3/2*fs$m1['A']*fs$m1['T']
dtor2 <- sum(c(l2.ntor, i2.ntor, c2.ntor))
expect <- matrix(c(c(c1.ntor, l1.ntor, i1.ntor)/dtor1,
                   c(c2.ntor, l2.ntor, i2.ntor)/dtor2),
                 ncol=3, byrow=T)
colnames(expect) <- causes
testthat::expect_equal(post$marg, expect)


# Example 6 (hand computed likelihood / Pv3R computed posterior)
y <- list(list(m1="A"), list(m1="A"), list(m1="A"))
post <- compute_posterior(y, fs)
ll_liklihood <- fs$m1['A']^3
cc_liklihood <- fs$m1['A']
testthat::expect_equal(ll_liklihood/post$joint["II"],
                       cc_liklihood/post$joint["CC"])


# Example 7
y <- list(list(m1=c("A","T")), list(m1="T"))
n.m <- 1
fmat <- rdirichlet(n.m, alphas)
colnames(fmat) <- alleles
fs <- setNames(lapply(1:n.m, function(x) fmat[x,]), paste0("m", 1:n.m))
post <- compute_posterior(y, fs)
l.ntor <- 5/18*fs$m1['A']*fs$m1['T']^2 + 1/4*fs$m1['A']*fs$m1['T']
i.ntor <- 3/4*fs$m1['A']*fs$m1['T']^2
c.ntor <- 3/8*fs$m1['A']*fs$m1['T']
dtor <- sum(c(l.ntor, i.ntor, c.ntor))
expect <- matrix(c(c.ntor, l.ntor, i.ntor)/dtor,
                 ncol=3, byrow=T)
colnames(expect) <- causes
testthat::expect_equal(post$marg, expect)

# Example 9
y <- list(list(m1 = c("A","T"), m2 = "T", m3 = c("C", "G")),
          list(m1 = "T", m2 = "T", m3 = "C"))
n.m <- 3
fmat <- rdirichlet(n.m, alphas)
colnames(fmat) <- alleles
fs <- setNames(lapply(1:n.m, function(x) fmat[x,]), paste0("m", 1:n.m))
post <- compute_posterior(y, fs)
l_ntor <- (1/9)*  (2* fs$m1["T"]*fs$m2["T"]^2*fs$m3["C"] +
               (3/8)* fs$m1["T"]*(fs$m2["T"]^2 + fs$m2["T"])*fs$m3["C"] +
               (1/8)* fs$m1["T"]*(fs$m2["T"]^2 + fs$m2["T"])*(fs$m3["C"] + 1) +
               (1/8)*(fs$m1["T"] + 1)*(fs$m2["T"]^2 + fs$m2["T"])*fs$m3["C"] +
               (1/8)*(fs$m1["T"] + 1)*(fs$m2["T"]^2 + fs$m2["T"])*(fs$m3["C"] + 1) +
              (1/32)*(3*fs$m2["T"] + 1) +
               (1/8)*(9*fs$m2["T"] + 1))
i_ntor <- (1/8)*fs$m1["T"]*fs$m2["T"]*fs$m3["C"]*(9*fs$m2["T"] + 1)
c_ntor <- (1/32)*(9*fs$m2["T"] + 1)

dtor <- sum(c(l_ntor, i_ntor, c_ntor))

expect <- matrix(c(c_ntor, l_ntor, i_ntor)/dtor,
                 ncol=3, byrow=T)
colnames(expect) <- causes
testthat::expect_equal(post$marg, expect)




