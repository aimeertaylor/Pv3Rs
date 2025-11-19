# Sample a relationship graph (RG)

Uses the techniques in
[`enumerate_RGs`](https://aimeertaylor.github.io/Pv3Rs/reference/enumerate_RGs.md)
to sample a single RG uniformly. All clonal partitions are generated,
each weighted by its number of consistent sibling partitions. A clonal
partition is sampled proportional to its weight, then a consistent
sibling partition is drawn uniformly. The resulting nested partition
represents the RG; see
[`enumerate_RGs`](https://aimeertaylor.github.io/Pv3Rs/reference/enumerate_RGs.md)
for details.

## Usage

``` r
sample_RG(MOIs, igraph = TRUE)
```

## Arguments

- MOIs:

  Vector of per-episode multiplicities of infection (MOIs), i.e.,
  numbers of per-episode genotypes / vertices.

- igraph:

  Logical; if `TRUE` (default), returns the RG as an `igraph` object.

## Value

An RG encoded either as an `igraph` object (default), or as a list; see
[`enumerate_RGs`](https://aimeertaylor.github.io/Pv3Rs/reference/enumerate_RGs.md)
for details.

## Examples

``` r
set.seed(1)
RG <- sample_RG(c(3, 2))
plot_RG(RG, vertex.label = NA)

```
