# WORK in progress - commit today

test_that("Manually computed examples work", {
  expect_equal(2 * 2, 4)
})

fs <- list(m1 = setNames(c(0.4, 0.6), c("A", "T")))
prior <- matrix(rep(1 / 3, 3), nrow = 1, dimnames = list(NA, c("C", "L", "I")))


# Single heterologous marker
y <- list(list(m1 = c("A")), list(m1 = c("T")))
expect_equal(
  suppressMessages(
    compute_posterior(y, fs, prior)$marg[1, c("C", "L", "I")]),
  c(C = 0, L = 1 / 3, I = 2 / 3)
)

# Single homologous marker
y <- list(list(m1 = c("A")), list(m1 = c("A")))
posterior_L <- 1 / 3
posterior_I <- (2 * fs$m1[["A"]]) / (3 * (fs$m1[["A"]] + 1))
posterior_C <- 2 / (3 * (fs$m1[["A"]] + 1))
expect_equal(
  suppressMessages(compute_posterior(y, fs)$marg[1, c("C", "L", "I")]),
  c(C = posterior_C, L = posterior_L, I = posterior_I)
)

# Multiple partially heterologous marker
fs <- list(
  m1 = setNames(c(0.4, 0.6), c("A", "T")),
  m2 = setNames(c(0.2, 0.8), c("T", "other")),
  m3 = setNames(c(0.3, 0.7), c("G", "other"))
)
y <- list(list(m1 = "A", m2 = "T", m3 = "G"), list(m1 = "T", m2 = "T", m3 = "G"))
evidence <- (1 / 3) * fs$m1[["A"]] * fs$m1[["T"]] * fs$m2[["T"]] * fs$m3[["G"]] *
  ((1 / 3) * (fs$m2[["T"]] * fs$m3[["G"]] + (1 / 8) * (fs$m2[["T"]] + 1) * (fs$m3[["G"]] + 1)) + fs$m2[["T"]] * fs$m3[["G"]])
posterior_L <- ((1 / 3) * (fs$m2[["T"]] * fs$m3[["G"]] + (1 / 8) * (fs$m2[["T"]] + 1) * (fs$m3[["G"]] + 1))) /
  ((1 / 3) * (fs$m2[["T"]] * fs$m3[["G"]] + (1 / 8) * (fs$m2[["T"]] + 1) * (fs$m3[["G"]] + 1)) + fs$m2[["T"]] * fs$m3[["G"]])
posterior_I <- (fs$m2[["T"]] * fs$m3[["G"]]) /
  ((1 / 3) * (fs$m2[["T"]] * fs$m3[["G"]] + (1 / 8) * (fs$m2[["T"]] + 1) * (fs$m3[["G"]] + 1)) + fs$m2[["T"]] * fs$m3[["G"]])
posterior_C <- 0
expect_equal(
  suppressMessages(compute_posterior(y, fs, prior)$marg[1, c("C", "L", "I")]),
  c(C = posterior_C, L = posterior_L, I = posterior_I)
)
