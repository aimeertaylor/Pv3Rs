################################################################################
# Run benchmarking results and save, to be reported in ../benchmark.Rmd
# Warning: takes a few hours to run
################################################################################

library(Pv3Rs)
library(gtools)
library(tictoc)

# simulate allele frequencies for one marker
sim_fs_single <- function(n_a) {
  alleles <- letters[1:n_a] # n_a <= 26
  fs_unnamed <- as.vector(rdirichlet(1, alpha = rep(1, n_a)))
  setNames(fs_unnamed, alleles)
}

# simulate allele frequencies for multiple markers
# assumes each marker has the same number of alleles
sim_fs <- function(n_m, n_a) {
  markers <- paste0("m", 1:n_m) # marker names
  n_a_vec <- setNames(rep(n_a, n_m), markers)
  lapply(n_a_vec, sim_fs_single)
}

# simulate data for one marker, one episode
# sample without replacement to match desired MOI, so must have length(fs) >= MOI
sim_data_single <- function(fs_single, MOI) {
  sample(names(fs_single), MOI, prob=fs_single)
}

# simulate data for all markers, all episodes
sim_data <- function(MOIs, n_m, n_a) {
  fs <- sim_fs(n_m, n_a)
  y <- lapply(
    MOIs,
    function(MOI) lapply(fs, sim_data_single, MOI) # sample data for one episode
  )
  return(list(y=y, fs=fs))
}

# benchmark for increasing number of markers
set.seed(1)
n_markers <- 50*(1:8)
n_a <- 8
MOIs <- c(2, 2)
n_rep <- 30

data_by_m <- list()
for (n_m in n_markers) {
  data_by_m[[as.character(n_m)]] <- lapply(
    1:n_rep,
    function(x) sim_data(MOIs, n_m, n_a)
  )
}

time_by_m <- list()
for (n_m in n_markers) {
  tic(msg = paste0("Running n_m=", n_m))
  time_by_m[[as.character(n_m)]] <- sapply(
    data_by_m[[as.character(n_m)]],
    function(data) {
      tic("")
      suppressMessages(compute_posterior(data$y, data$fs))
      res <- toc(quiet=TRUE)
      res$toc - res$tic
    }
  )
  toc()
}

# benchmark for various MOI
set.seed(1)
n_m <- 1
n_a <- 8
MOIs_comb <- list()
for(MOItot in 3:7) {
  for(MOI1 in (MOItot-1):1) {
    MOI2 <- MOItot-MOI1
    if(MOI1 >= MOI2) MOIs_comb <- c(MOIs_comb, list(c(MOI1, MOI2)))
  }
}
n_rep <- 30

data_by_MOI <- list()
for (MOIs in MOIs_comb) {
  data_by_MOI[[paste(MOIs, collapse=",")]] <- lapply(
    1:n_rep,
    function(x) sim_data(MOIs, n_m, n_a)
  )
}

time_by_MOI <- list()
for (MOIs in MOIs_comb) {
  MOIstr <- paste(MOIs, collapse=",")
  tic(msg = paste0("Running MOI=", MOIstr))
  time_by_MOI[[MOIstr]] <- sapply(
    data_by_MOI[[MOIstr]],
    function(data) {
      tic("")
      suppressMessages(compute_posterior(data$y, data$fs))
      res <- toc(quiet=TRUE)
      res$toc - res$tic
    }
  )
  toc()
}

# benchmark for increasing number of monoclonal episodes
set.seed(1)
n_m <- 1
n_a <- 8
n_epis <- 2:7
n_rep <- 30

data_by_epi <- list()
for (epi in n_epis) {
  data_by_epi[[as.character(epi)]] <- lapply(
    1:n_rep,
    function(x) sim_data(rep(1, epi), n_m, n_a)
  )
}

time_by_epi <- list()
for (epi in n_epis) {
  tic(msg = paste0("Running epi=", epi))
  time_by_epi[[as.character(epi)]] <- sapply(
    data_by_epi[[as.character(epi)]],
    function(data) {
      tic("")
      suppressMessages(compute_posterior(data$y, data$fs))
      res <- toc(quiet=TRUE)
      res$toc - res$tic
    }
  )
  toc()
}

# benchmark for increasing number of episodes for fixed MOI = 7
set.seed(1)
n_m <- 1
n_a <- 8
maxepi <- 7 # fixed total MOI
n_rep <- 30

MOI_fix_comb <- list()
for(epi in 2:maxepi) {
  m <- maxepi %/% epi
  rem <- maxepi %% epi
  MOI_fix_comb <- c(MOI_fix_comb, list(c(rep(m+1, rem), rep(m, epi-rem))))
}

data_by_MOI_fix <- list()
for (MOIs in MOI_fix_comb) {
  data_by_MOI_fix[[paste(MOIs, collapse=",")]] <- lapply(
    1:n_rep,
    function(x) sim_data(MOIs, n_m, n_a)
  )
}

time_by_MOI_fix <- list()
for (MOIs in MOI_fix_comb) {
  MOIstr <- paste(MOIs, collapse=",")
  tic(msg = paste0("Running MOI=", MOIstr))
  time_by_MOI_fix[[MOIstr]] <- sapply(
    data_by_MOI_fix[[MOIstr]],
    function(data) {
      tic("")
      suppressMessages(compute_posterior(data$y, data$fs))
      res <- toc(quiet=TRUE)
      res$toc - res$tic
    }
  )
  toc()
}

# save all timing results
output_benchmark = list(
  time_by_m=time_by_m,
  time_by_MOI=time_by_MOI,
  time_by_epi=time_by_epi,
  time_by_MOI_fix=time_by_MOI_fix
)
save(output_benchmark, file = "output_benchmark.rda")
