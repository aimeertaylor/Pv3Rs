library(Pv3Rs)

# =====================================================================
# Some allele frequencies
fs <- list(m1=c('1'=0.78, '2'=0.14, '3'=0.07, '4'=0.005, '5' = 0.005),
           m2=c('1'=0.27, '2'=0.35, '3'=0.38),
           m3=c('1'=0.55, '2'=0.45))


# =====================================================================
# Using a dummy marker to encode MOIs
y0 <- list(enrol = list(m1=c('1')),
           recur = list(m1=c('1','1')))

y1 <- list(enrol = list(m1=c('1'), m2=NA),
           recur = list(m1=c('1'), m2=c(NA,NA)))

compute_posterior(y0, fs)
compute_posterior(y1, fs)



# =====================================================================
# Multiple homoallelic observations at homoallelic alleles make no difference
# when the MOI is otherwise encoded:
y0 <- list(enrol = list(m1=c('1'), m2=c('1','2')),
           recur = list(m1=c('1'), m2=c('1')))

y1 <- list(enrol = list(m1=c('1','1'), m2=c('1','2')),
           recur = list(m1=c('1'), m2=c('1')))

compute_posterior(y0, fs)
compute_posterior(y1, fs)

# =====================================================================
# Multiple homoallelic observations at heteroallelic alleles provide a
# rough way to encode within-sample allele frequencies
y0 <- list(enrol = list(m1=c('1','1','2'), m2=c('1','2','3')),
           recur = list(m1=c('1'), m2=c('1')))

y1 <- list(enrol = list(m1=c('1','2'), m2=c('1','2','3')),
           recur = list(m1=c('1'), m2=c('1')))

y2 <- list(enrol = list(m1=c('1','2','2'), m2=c('1','2','3')),
           recur = list(m1=c('1'), m2=c('1')))

compute_posterior(y0, fs)
compute_posterior(y1, fs)
compute_posterior(y2, fs)

