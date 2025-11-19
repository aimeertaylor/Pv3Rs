# Determine multiplicities of infection (MOIs)

Returns a MOI estimate for each episode based on allelic diversity
across markers.

## Usage

``` r
determine_MOIs(y, return.names = FALSE)
```

## Arguments

- y:

  List of lists encoding allelic data; see
  [`compute_posterior`](https://aimeertaylor.github.io/Pv3Rs/reference/compute_posterior.md)
  for more details. The outer list contains episodes in chronological
  order. The inner list contains named markers per episode. For each
  marker, one must specify an allelic vector: a set of distinct alleles
  detected at that marker; or `NA` if marker data are missing.

- return.names:

  Logical; if TRUE and `y` has named episodes, episode names are
  returned.

## Value

Numeric vector containing one MOI estimate per episode, each estimate
representing the maximum number of distinct alleles observed at any
marker per episode.

## Details

A true MOI is a number of genetically distinct groups of clonal
parasites within an infection. Give or take *de novo* mutations, all
parasites within a clonal group share the same DNA sequence, which we
call a genotype. As such, MOIs are distinct parasite genotype counts.
Under the Pv3Rs model assumption that there are no genotyping errors,
the true MOI of an episode is greater than or equal to the maximum
distinct allele count for any marker in the data on that episode. In
other words, under the assumption of no genotyping errors, maximum
distinct allelic counts are the most parsimonious MOI estimates
compatible with the data. By default, these MOI estimates are used by
[`compute_posterior`](https://aimeertaylor.github.io/Pv3Rs/reference/compute_posterior.md).

## Examples

``` r
y <- list(enrol = list(m1 = c("A", "B"), m2 = c("A"), m3 = c("C")),
          recur = list(m1 = c("B"), m2 = c("B", "C"), m3 = c("A", "B", "C")))
determine_MOIs(y) # returns c(2, 3)
#> [1] 2 3
```
