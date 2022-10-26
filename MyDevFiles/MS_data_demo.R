library(Pv3Rs)
library(plyr)
library(gtools) # for rdirichlet
library(tictoc)

# downloaded data from GitHub
load('data/MS_data_PooledAnalysis.RData')

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

tic()
# this is currently a matrix, is data frame better?
Post_probs <- do.call(rbind, lapply(1:length(yns), tmp <- function(i) {
  writeLines(paste("Individual", i, "out of", length(yns)))
  y.df <- yns[[i]]

  # no recurrence, skip
  if(length(unique(y.df$Episode_Identifier)) == 1) return(NULL)
  # cases that are too complex, can change this criteria
  if(nrow(y.df) > 8) return(NULL)

  ynts <- dlply(y.df, 'Episode_Identifier')

  # (?) only use markers for which every episode has at least 1 non-NA
  ms <- NULL
  for(m in MSs_all) {
    if(all(sapply(ynts, function(ynt) !all(is.na(ynt[m]))))) {
      ms <- c(ms, m)
    }
  }

  if(length(ms) == 0) return(NULL)

  # ensure that Episode is non-decreasing
  stopifnot(all(order(y.df$Episode) == 1:nrow(y.df)))

  # transform data frame format to format taken by 'compute_posterior'
  y <- lapply(ynts, function(ynt) {
    setNames(lapply(ms, function(m) ynt[m][!is.na(ynt[m])]), ms)
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
