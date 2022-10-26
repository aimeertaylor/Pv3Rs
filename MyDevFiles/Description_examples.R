library(Pv3Rs)
library(gtools)

set.seed(1)

causes <- c("C", "L", "I")
alleles <- c("A", "C", "G", "T")
n.as <- length(alleles)
alphas <- rep(1, n.as)


# Example 1

y <- list(list(m1="A"), list(m1="T"))
n.m <- 1
fmat <- rdirichlet(n.m, alphas)
colnames(fmat) <- alleles
fs <- setNames(lapply(1:n.m, function(x) fmat[x,]), paste0("m", 1:n.m))

print(compute_posterior(y, fs)$marg)
writeLines("Manual calculation:")
print(setNames(c(0, 1/3, 2/3), causes))


# Example 2

y <- list(list(m1="A"), list(m1="A"))
n.m <- 1
fmat <- rdirichlet(n.m, alphas)
colnames(fmat) <- alleles
fs <- setNames(lapply(1:n.m, function(x) fmat[x,]), paste0("m", 1:n.m))

print(compute_posterior(y, fs)$marg)
writeLines("Manual calculation:")
print(setNames(c(2/3/(1+fs$m1['A']),
                 1/3,
                 2/3*fs$m1['A']/(1+fs$m1['A'])),
               causes))


# Example 3
y <- list(list(m1="A", m2="T", m3="G"),
          list(m1="T", m2="T", m3="G"))
n.m <- 3
fmat <- rdirichlet(n.m, alphas)
colnames(fmat) <- alleles
fs <- setNames(lapply(1:n.m, function(x) fmat[x,]), paste0("m", 1:n.m))

print(compute_posterior(y, fs)$marg)
writeLines("Manual calculation:")
l.ntor <- (fs$m2['T']*fs$m3['G']+(fs$m2['T']+1)*(fs$m3['G']+1)/8)/3
i.ntor <- fs$m2['T']*fs$m3['G']
print(setNames(c(0,
                 l.ntor/(l.ntor+i.ntor),
                 i.ntor/(l.ntor+i.ntor)),
               causes))


# Example 4
y <- list(list(m1="A", m2="T", m3="G"),
          list(m1="A", m2="T", m3="G"))
n.m <- 3
fmat <- rdirichlet(n.m, alphas)
colnames(fmat) <- alleles
fs <- setNames(lapply(1:n.m, function(x) fmat[x,]), paste0("m", 1:n.m))

print(compute_posterior(y, fs)$marg)
writeLines("Manual calculation:")
l.ntor <- (fs$m1['A']+1)*(fs$m2['T']+1)*(fs$m3['G']+1)/8 + fs$m1['A']*fs$m2['T']*fs$m3['G'] + 1
i.ntor <- 3*fs$m1['A']*fs$m2['T']*fs$m3['G']
c.ntor <- 3
dtor <- l.ntor + i.ntor + c.ntor
print(setNames(c(c.ntor/dtor,
                 l.ntor/dtor,
                 i.ntor/dtor),
               causes))


# Example 5

y <- list(list(m1="A"), list(m1="T"), list(m1="T"))
n.m <- 1
fmat <- rdirichlet(n.m, alphas)
colnames(fmat) <- alleles
fs <- setNames(lapply(1:n.m, function(x) fmat[x,]), paste0("m", 1:n.m))

print(compute_posterior(y, fs)$marg)
writeLines("Manual calculation:")
l.ntor <- 17/24*fs$m1['A']*fs$m1['T']^2 + 11/16*fs$m1['A']*fs$m1['T']
i.ntor <- 7/5*fs$m1['A']*fs$m1['T']^2 + 13/10*fs$m1['A']*fs$m1['T']
c.ntor <- 0
dtor <- l.ntor + i.ntor + c.ntor
print(setNames(c(c.ntor/dtor,
                 l.ntor/dtor,
                 i.ntor/dtor),
               causes))
l.ntor <- 73/120*fs$m1['A']*fs$m1['T']^2 + 39/80*fs$m1['A']*fs$m1['T']
i.ntor <- 3/2*fs$m1['A']*fs$m1['T']^2
c.ntor <- 3/2*fs$m1['A']*fs$m1['T']
dtor <- l.ntor + i.ntor + c.ntor
print(setNames(c(c.ntor/dtor,
                 l.ntor/dtor,
                 i.ntor/dtor),
               causes))


# Example 7

y <- list(list(m1=c("A","T")), list(m1="T"))
n.m <- 1
fmat <- rdirichlet(n.m, alphas)
colnames(fmat) <- alleles
fs <- setNames(lapply(1:n.m, function(x) fmat[x,]), paste0("m", 1:n.m))

print(compute_posterior(y, fs)$marg)
writeLines("Manual calculation:")
l.ntor <- 5/18*fs$m1['A']*fs$m1['T']^2 + 1/4*fs$m1['A']*fs$m1['T']
i.ntor <- 3/4*fs$m1['A']*fs$m1['T']^2
c.ntor <- 3/8*fs$m1['A']*fs$m1['T']
dtor <- l.ntor + i.ntor + c.ntor
print(setNames(c(c.ntor/dtor,
                 l.ntor/dtor,
                 i.ntor/dtor),
               causes))