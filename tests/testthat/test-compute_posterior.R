# WORK in progress - commit today

test_that("Manually computed examples work", {
  expect_equal(2 * 2, 4)
})

fs <- list(m1=setNames(c(0.4, 0.6), c("A","B")))
prior <- matrix(rep(1/3,3), nrow=1)


# Single heterologous marker - works
y <- list(list(m1=c("A")), list(m1=c("B")))
expect_equal(compute_posterior(y, fs, prior)$marg[1,c("C", "L", "I")],
             c(C = 0, L = 1/3, I = 2/3))

# Single homologous marker - error in RG_inference
y <- list(list(m1=c("A")), list(m1=c("A")))
posterior_L <- 1/3
posterior_I <- (2*fs$m1["A"])/(3*(fs$m1["A"] + 1))
posterior_C <- 2/(3*(fs$m1["A"] + 1))
expect_equal(compute_posterior(y, fs)$marg[1,c("C", "L", "I")],
             c(C = posterior_C, L =  posterior_L, I = posterior_I))

# Multiple partially heterologous marker - same error, dropping a dimension
fs <- list(m1=setNames(c(0.4, 0.6), c("A","B")),
           m2=setNames(c(0.2, 0.8), c("C","D")),
           m3=setNames(c(0.3, 0.7), c("E","F")))
y <- list(list(m1="A", m2="C", m3="E"), list(m1="B", m2="C", m3="E"))
evidence <- (1/3)*fs$m1["A"]*fs$m1["B"]*fs$m2["C"]*fs$m3["E"]*
  ((1/3)*(fs$m2["C"]*fs$m3["E"] + (1/8)*(fs$m2["B"] + 1)*(fs$m3["G"] + 1)) + fs$m2["C"]*fs$m3["E"])
posterior_L <- ((1/3)*(fs$m2["C"]*fs$m3["E"] + (1/8)*(fs$m2["B"] + 1)*(fs$m3["G"] + 1))) /
  ((1/3)*(fs$m2["C"]*fs$m3["E"] + (1/8)*(fs$m2["B"] + 1)*(fs$m3["G"] + 1)) + fs$m2["C"]*fs$m3["E"])
posterior_I <- (fs$m2["C"]*fs$m3["E"])/
  ((1/3)*(fs$m2["C"]*fs$m3["E"] + (1/8)*(fs$m2["B"] + 1)*(fs$m3["G"] + 1)) + fs$m2["C"]*fs$m3["E"])
posterior_C <- 0
expect_equal(compute_posterior(y, fs, prior)$marg[1,c("C", "L", "I")],
             c(C = posterior_C, L =  posterior_L, I = posterior_I))

