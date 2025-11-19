# Allele frequencies computed using example *Plasmodium vivax* data

The posterior mean of a multinomial-Dirichlet model with uniform prior
fit to data on allele prevalence in initial episodes of
[ys_VHX_BPD](https://aimeertaylor.github.io/Pv3Rs/reference/ys_VHX_BPD.md).
Because the model is fit to allele prevalence (observed) not allele
frequency ( requires integrating-out unknown multiplicities of
infection) it is liable to underestimate the frequencies of common
alleles and overestimate those of rare but detected alleles.

## Usage

``` r
fs_VHX_BPD
```

## Format

A list of nine markers; for each marker a named vector of allele
frequencies that sum to one.

## Source

- MS_data_PooledAnalysis.RData downloaded from
  <https://zenodo.org/records/3368828>

- <https://github.com/aimeertaylor/Pv3Rs/blob/main/data-raw/fs_VHX_BPD.R>
