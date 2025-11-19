# Example *Plasmodium vivax* data

Previously-published microsatellite data on *P. vivax* parasites
extracted from study participants enrolled in the Best Primaquine Dose
(BPD) and Vivax History (VHX) trials; see Taylor & Watson et al. 2019
([doi:10.1038/s41467-019-13412-x](https://doi.org/10.1038/s41467-019-13412-x)
) for more details of the genetic data; for more details of the VHX and
BPD trials, see Chu et al. 2018a
([doi:10.1093/cid/ciy319](https://doi.org/10.1093/cid/ciy319) ) and Chu
et al. 2018b
([doi:10.1093/cid/ciy735](https://doi.org/10.1093/cid/ciy735) ).

## Usage

``` r
ys_VHX_BPD
```

## Format

A list of 217 study participants; for each study participant, a list of
one or more episodes; for each episode, a list of three or more
microsatellite markers; for each marker, a vector of observed alleles
(repeat lengths). For example:

- BPD_103:

  Study participant identifier: study participant 103 in the BPD trial

- BPD_103_1:

  Episode identifier: episode one of study participant 103 in the BPD
  trial

- PV.3.27:

  Marker identifier: *P. vivax* 3.27

- 18:

  Allele identifier: 18 repeat lengths

## Source

- MS_data_PooledAnalysis.RData downloaded from
  <https://zenodo.org/records/3368828>

- <https://github.com/aimeertaylor/Pv3Rs/blob/main/data-raw/ys_VHX_BPD.R>
