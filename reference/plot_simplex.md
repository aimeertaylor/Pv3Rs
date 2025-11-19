# Plots a 2D simplex

Plots a 2D simplex (a triangle with unit sides centered at the origin)
onto which per-recurrence posterior probabilities of recrudescence,
relapse, reinfection (or any other probability triplet summing to one)
can be projected; see
[`project2D()`](https://aimeertaylor.github.io/Pv3Rs/reference/project2D.md)
and **Examples** below.

## Usage

``` r
plot_simplex(
  v.labels = c("Recrudescence", "Relapse", "Reinfection"),
  v.cutoff = 0.5,
  v.colours = c("yellow", "purple", "red"),
  plot.tri = TRUE,
  lim.mar = 0.1,
  p.coords = NULL,
  p.labels = rownames(p.coords),
  p.labels.pos = 3,
  p.labels.cex = 1,
  ...
)
```

## Arguments

- v.labels:

  Vertex labels anticlockwise from top (default: "Recrudescence",
  "Relapse", "Reinfection"). If NULL, vertices are not labelled.

- v.cutoff:

  Number between 0.5 and 1 that separates lower vs higher probability
  regions. Use with caution for recrudescence and reinfection
  classification; see [Understand posterior
  probabilities](https://aimeertaylor.github.io/Pv3Rs/articles/posterior-probabilities.html).

- v.colours:

  Vertex colours anticlockwise from top.

- plot.tri:

  Logical; draws the triangular boundary if `TRUE` (default).

- lim.mar:

  Margin away from simplex for axes limits; defaults to 0.1.

- p.coords:

  Matrix of 3D simplex coordinates (e.g., per-recurrence probabilities
  of recrudescence, relapse and reinfection), one vector of 3D
  coordinates per row, each row is projected onto 2D coordinates using
  [`project2D()`](https://aimeertaylor.github.io/Pv3Rs/reference/project2D.md)
  and plotted as a single simplex point using
  [`graphics::points()`](https://rdrr.io/r/graphics/points.html). If the
  user provides a vector encoding a probability triplet summing to one,
  it is converted to a matrix with one row.

- p.labels:

  Labels of points in `p.coords` (default row names of `p.coords`) No
  labels if `NA`.

- p.labels.pos:

  Position of `p.labels`: `1` = below, `2` = left, `3` = above (default)
  and `4` = right. Can be a single value or a vector.

- p.labels.cex:

  Size expansion of `p.labels` passed to
  [`text`](https://rdrr.io/r/graphics/text.html).

- ...:

  Additional parameters passed to
  [`graphics::points()`](https://rdrr.io/r/graphics/points.html).

## Value

None

## Examples

``` r
# Running example (runs across compute_posterior, plot_data and plot_simplex)
# based on real data from chloroquine-treated participant 52 of the Vivax
# History Trial (Chu et al. 2018a, https://doi.org/10.1093/cid/ciy319)
y <- ys_VHX_BPD[["VHX_52"]] # y is a list of length two (two episodes)
post <- compute_posterior(y, fs_VHX_BPD, progress.bar = FALSE)
#> Number of valid relationship graphs (RGs) is 971
#> Computing log p(Y|RG) for 971 RGs
#> Finding log-likelihood of each vector of recurrence states
#> 
plot_simplex(p.coords = post$marg, p.labels = "", pch = 20, cex = 2)


# Basic example
plot_simplex(p.coords = diag(3),
             p.labels = c("(1,0,0)", "(0,1,0)", "(0,0,1)"),
             p.labels.pos = c(1,3,3))


# ==============================================================================
# Given data on an enrollment episode and a recurrence,
# compute the posterior probabilities of the 3Rs and plot the deviation of the
# posterior from the prior
# ==============================================================================

# Some data:
y <- list(list(m1 = c('a', 'b'), m2 = c('c', 'd')), # Enrollment episode
          list(m1 = c('a'), m2 = c('c'))) # Recurrent episode

# Some allele frequencies:
fs <- list(m1 = setNames(c(0.4, 0.6), c('a', 'b')),
           m2 = setNames(c(0.2, 0.8), c('c', 'd')))

# A vector of prior probabilities:
prior <- array(c(0.2, 0.3, 0.5), dim = c(1,3),
               dimnames = list(NULL, c("C", "L", "I")))

# Compute posterior probabilities
post <- compute_posterior(y, fs, prior, progress.bar = FALSE)
#> Number of valid relationship graphs (RGs) is 9
#> Computing log p(Y|RG) for 9 RGs
#> Finding log-likelihood of each vector of recurrence states
#> 

# Plot simplex with the prior and posterior
plot_simplex(p.coords = rbind(prior, post$marg),
             p.labels = c("Prior", "Posterior"),
             pch = 20)

# Add the deviation between the prior and posterior: requires obtaining 2D
# coordinates manually
xy_prior <- project2D(as.vector(prior))
xy_post <- project2D(as.vector(post$marg))
arrows(x0 = xy_prior["x"], x1 = xy_post["x"],
       y0 = xy_prior["y"], y1 = xy_post["y"], length = 0.1)
```
