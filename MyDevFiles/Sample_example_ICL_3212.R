library(tictoc)
library(Pv3Rs)
set.seed(3)

MOIs <- c(3,2,1,2)
infection_count <- length(MOIs)
gs_count <- sum(MOIs)
gs <- paste0("g", 1:gs_count) # Genotype names
# List time point indices per genotype
ts_per_gs <- rep(1:infection_count, MOIs)
# List genotypes per time point and time point pair
gs_per_ts <- split(gs, ts_per_gs)

# 4 possible alleles at each of 3 markers
n_a <- 4
n_m <- 3
alleles <- c('A', 'C', 'G', 'T')
ms <- paste0("m", 1:n_m) # Marker names

# generate an interesting relationship graph
if(file.exists("data/RG_ICL_3212.rds")) {
  RG0 <- readRDS("data/RG_ICL_3212.rds")
} else {
  tic()
  RGs <- enumerate_RGs_alt(MOIs, igraph=F)
  toc()
  # look for an RG where first two genotypes are strangers
  # infection 2 is reinfection
  # infection 3 is recrudescence
  # infection 4 is relapse
  for(RG0 in RGs) {
    if(RG0$sib.vec[RG0$clone.vec[1]] != RG0$sib.vec[RG0$clone.vec[2]]
       & length(intersect(RG0$sib.vec[RG0$clone.vec[1:3]],
                          RG0$sib.vec[RG0$clone.vec[4:5]]))==0
       & RG0$clone.vec[6] %in% RG0$clone.vec[4:5]) {
      break
    }
  }
  RG0 <- RG_to_igraph(RG0, gs, ts_per_gs)
  saveRDS(RG0, "data/RG_ICL_3212.rds")
}

plot_RG(RG0, edge.curved=0.2)

MOI_from_geno <- function(g) {
  if(is.vector(g)) 1
  else max(apply(g, 2, function(x) length(unique(x))))
}

set.seed(4)
# keep sampling data until MOIs are indeed the desired MOIs
while(TRUE) {
  fs <- setNames(lapply(ms, function(m) setNames(as.vector(MCMCpack::rdirichlet(1,
                                                                                rep(1, n_a))), alleles)),
                 ms)
  n.lines <- 2*sum(MOIs)
  g.lines <- sapply(fs, function(f) sample(alleles, n.lines,
                                           replace=TRUE, prob=f))
  rownames(g.lines) <- letters[1:n.lines]
  lines <- sample_lineages(RG0, n_m)
  genos <- sapply(1:n_m, function(i) g.lines[lines[,i],i])
  rownames(genos) <- gs
  colnames(genos) <- ms

  genos_MOI <- sapply(gs_per_ts, function(g) MOI_from_geno(genos[g,]))
  if(any(MOIs != genos_MOI)) next

  y <- lapply(gs_per_ts,
              function(g) {
                if(length(g)==1) as.list(genos[g,])
                else apply(genos[g,], 2, unique)
              })

  # find example that is harder to phase
  if(sum(sapply(y[[1]], length))==8 & sum(sapply(y[[2]], length))==5) break
}

saveRDS(y, "data/y_ICL_3212.rds")
saveRDS(fs, "data/fs_ICL_3212.rds")
