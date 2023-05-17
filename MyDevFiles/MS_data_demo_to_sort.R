library(Pv3Rs)
library(plyr)
library(gtools) # for rdirichlet
library(tictoc)

# downloaded data from GitHub
load('data/MS_data_PooledAnalysis.RData')

# Update so data are either available as R package or this script is stored elsewhere
load('../RecurrentVivax/RData/GeneticModel/MS_data_PooledAnalysis.RData')

# copied from Pooled_Analysis.Rmd
MSs_all <- c("PV.3.502","PV.3.27","PV.ms8",
             "PV.1.501","PV.ms1","PV.ms5",
             "PV.ms6","PV.ms7","PV.ms16")

D_weight_Prior <- 1

Alpha_Posteriors <- apply(MS_pooled[,MSs_all], 2, function(x, Ind_Primary){
  # Extract xmax
  xmax <-  max(x,na.rm=T)
  # prior parameter vector (iterpolates unobserved repeat lengths < xmax)
  param_vector <- array(D_weight_Prior, dim=xmax, dimnames=list(1:xmax))
  # observed data summarised as counts
  obs_counts <-  table(x[Ind_Primary])
  # posterior parameter vector
  param_vector[names(obs_counts)] <-  param_vector[names(obs_counts)] + obs_counts
  return(param_vector)
})

Fs_Combined <- sapply(Alpha_Posteriors, function(x){x/sum(x)})

Fs_random <- lapply(Alpha_Posteriors, function(alpha) as.vector(rdirichlet(n=1, alpha)))

# this dlply is taken from post_prob_CLI
yns <- dlply(MS_pooled, 'ID')

#===============================================
# Pushing the limits
# part.list is what breaks it
# By trial an error, 10 seems to be the maximum allowed: does 10 work for monoclonal only?
# What is the maximum number of genotypes? 18!
genotype_count_store <- sapply(1:length(yns), function(i) {
  y.df <- yns[[i]]
  nrow(y.df)})

# genotype counts above 8:
genotype_count_store[genotype_count_store > 10]
which(genotype_count_store == 10)
# Can compute_posterior deal with 10 genotypes?

tic()
y.df <- yns[[101]]
# no recurrence, skip
if(length(unique(y.df$Episode_Identifier)) == 1) stop("no recurrence")

# only use markers for which there is at least one NA
ms <- MSs_all[apply(!is.na(y.df[MSs_all]), 2, any)]
if(length(ms) == 0) stop("no data")

ynts <- dlply(y.df, 'Episode_Identifier')
# ensure that Episode is non-decreasing
if(!all(order(y.df$Episode) == 1:nrow(y.df))) {
  stop(paste("Individual", i, "needs rows to be reordered according to infection number"))
}

# transform data frame format to format taken by 'compute_posterior'
y <- lapply(ynts, function(ynt) {
  setNames(lapply(ms, function(m) {
    alleles <- ynt[m][!is.na(ynt[m])] # extract non-NAs
    if(length(alleles) > 0) return(alleles)
    return(NA) # if all are NAs, change empty vector to NA
  }), ms)
})

post <- compute_posterior(y, Fs_Combined) # can change to Fs_random
post$marg
toc()
#===============================================




tic()
# this is currently a matrix, is data frame better?
Post_probs <- do.call(rbind, lapply(1:length(yns), function(i) {
  writeLines(paste("Individual", i, "out of", length(yns)))
  y.df <- yns[[i]]

  # no recurrence, skip
  if(length(unique(y.df$Episode_Identifier)) == 1) return(NULL)
  # >8 genotypes considered too complex, can change this criteria
  if(nrow(y.df) > 8) return(NULL)

  # only use markers for which there is at least one NA
  ms <- MSs_all[apply(!is.na(y.df[MSs_all]), 2, any)]
  if(length(ms) == 0) return(NULL)

  ynts <- dlply(y.df, 'Episode_Identifier')
  # ensure that Episode is non-decreasing
  if(!all(order(y.df$Episode) == 1:nrow(y.df))) {
    writeLines(paste("Individual", i, "needs rows to be reordered according to infection number"))
    return(NULL)
  }

  # transform data frame format to format taken by 'compute_posterior'
  y <- lapply(ynts, function(ynt) {
    setNames(lapply(ms, function(m) {
      alleles <- ynt[m][!is.na(ynt[m])] # extract non-NAs
      if(length(alleles) > 0) return(alleles)
      return(NA) # if all are NAs, change empty vector to NA
      }), ms)
    })

  post <- compute_posterior(y, Fs_Combined) # can change to Fs_random
  post$marg
}))
toc()

saveRDS(Post_probs, "data/Post_probs.rds")

postmat <- readRDS("data/Post_probs.rds")

hist(postmat[,1], breaks=100) # clonal
hist(postmat[,2], breaks=100) # relapse
hist(postmat[,3], breaks=100) # reinfection
