# Pv3Rs <img src="man/figures/logo.png" align="right" alt="" width="120" />

An R package for *Plasmodium vivax* molecular correction: statistical genetic
inference of *P. vivax*

[//]: # (use same order as software note abstract)

- Relapse
- Recrudescence
- Reinfection

The core function `compute_posterior()` can be used to obtain per-person
posterior probabilities of relapse, recrudescence, and reinfection (recurrence
states) using *P. vivax* genetic data across multiple episodes to update a prior 
over recurrence states, which is ideally informative (e.g., based on 
time-to-event information).

Two other important features are more general:

- `plot_data()` can be used to visualise genetic data for molecular
correction regardless of the analytical method; for example, *Plasmodium
falciparum* data intended for analysis using a WHO match-counting algorithm.

- `plot_simplex()` can be used to visualise per-recurrence posterior 
probabilities of relapse, recrudescence, and reinfection, or a vector of any three
numbers in zero to one that sum to one.

## Please be aware of the following points!

### Pre-release state: 

The R package is still under development and may contain some errors. Moreover,
the statistical model around which the Pv3Rs R package is being developed is not
yet peer-reviewed and thus liable to modification. The model is documented in
the preprint [1]. The model builds upon a prototype documented in [2] (see [1]
to better understand how the model underpinning Pv3Rs differs from the
prototype).

[1] [Taylor, Foo & White, 2022](https://www.medrxiv.org/content/10.1101/2022.11.23.22282669v1)

[2] [Taylor & Watson et al. 2019](https://www.nature.com/articles/s41467-019-13412-x)

### Prior considerations: 

- Genetic data are modelled using a Bayesian model, whose prior is ideally
informative (in [2] prior estimates were generated by a time-to-event model
built by James Watson) because the cause of recurrent *P. vivax* malaria is not
always identifiable from genetic data alone: when the data are consistent with
recurrent parasites that are relatively unrelated to those in all preceding
infections, both reinfection and relapse are plausible; meanwhile, when the data
are compatible with recurrent parasites that are clones of those in the
preceding infection, both recrudescence and relapse are plausible.

- The main Pv3Rs function, `compute_posterior()`, could be used to estimate the
probable cause of recurrent *P. falciparum* malaria by setting the prior
probability of relapse to zero. However, the current version of Pv3Rs is
suboptimal for *P. falciparum* molecular correction: genotyping errors, which
are not accounted for under the current Pv3Rs model, are liable to lead to the
misclassification of recrudescence as reinfection when the probability of
relapse is zero *a priori* (and of recrudescence as relapse when the prior
probability of relapse exceeds zero). 

### Notable assumptions and limitations: 

As with any model, Pv3Rs makes various assumptions that limit its capabilities in
some settings. They are explained in some detail below and summarised in a
table.

#### Mutually exclusive recurrent states
We model recurrent states (relapse, recrudescence, and reinfection) as mutually
exclusive. Pv3Rs was designed with a view towards clinical trials where study 
participants are actively followed up frequently and where all detected 
infections are treated to the extent that post-treatment parasitaemia drops 
below some detectable level before recurrence, if it occurs. In studies where 
infections persist untreated and/or where events have time to accumulate, the 
concept of recurrence is ill-defined, and the output of Pv3Rs might not be 
meaningful.

#### Complexities of molecular correction that exceed data sampled
We do not model all the complexities around molecular correction. For example, we do 
not model any population structure (e.g., household effects) or the failure to 
capture a low-density clone in the body by a blood sample of limited volume [[Snounou
and Beck, 1998]](https://doi.org/10.1016/S0169-4758(98)01340-4), nor hidden biomass in 
the spleen and bone marrow (an alternative source of P. vivax recurrence) [[Markus, 2019]](https://doi.org/10.1016/j.pt.2019.08.009). 
As with study design above, we urge  the user to interpret Pv3Rs’ output conditional 
on the input. For example, we expect Pv3Rs to output probable relapse if a person is 
reinfected by a new mosquito with parasites that are recently related to those that 
caused a previous infection, as might happen in household transmission chains.

#### Sibling misspecification
Relapsing parasites that are siblings of parasites in previous infections can be
meiotic, parent-child-like, regular or half siblings, but we model all sibling
parasites as regular siblings via the following assumptions: 
    - allele inheritance is independent (not true of meiotic siblings)
    - sibling relationships are transitive (not true of parent-child-like trios or some half-sibling trios)
    - alleles of a sibling cluster are drawn from at most two parental alleles (not true of half siblings).

In our experience, half sibling misspecification leads to some misclassification
of relapses as reinfections; see vignette ["Understand half-sibling misspecification"](https://aimeertaylor.github.io/Pv3Rs/articles/understand-half-sibs.html). 
A descriptive study to explore the extent of half-sibling misspecification is 
recommended (an example will be provided in an upcoming manuscript).

#### Observation errors and *de novo* mutations
We do not model undetected alleles, genotyping errors, or *de novo* mutations. 
Recall that recrudescent parasites are modelled as perfect clones under Pv3Rs. As 
such, the posterior probability of recrudescence is rendered zero by errors and 
mutations. This becomes more likely when there are data on more markers. Sensitivity 
analyses that explore the impact of errors and mutations on recurrence state 
estimates are merited.

#### Interpreting probable reinfection and recrudescence
When data are not sufficiently informative to distinguish between recrudescence and 
relapse (or reinfection and relapse), the posterior probabilities of recrudescence and 
relapse (or reinfection and relapse) are heavily influenced by assumption that the distribution
over graphs is uniform *a priori* (see vignette ["Understand graph prior ramifications"](https://aimeertaylor.github.io/Pv3Rs/articles/understand-graph-prior.html)). 
The development of a more biologically-principled generative model on parasite 
relationships is merited.

Limitation | Reason
----------- | ------
Possible misclassification of persistent and/or accumulated states as relapse | Modelling recurrent states as mutually exclusive
Possible inconsistency with data on more-and-more markers | Not modelling errors
Possible misclassification of relapse as reinfection | Half-sibling misspecification and not modelling errors
Possible misclassification of recrudescence | Not modelling errors
Possible misclassification of reinfection | Not modelling population structure
Strong prior impact on posterior probabilities | Recurrent states are not always identifiable from genetic data alone


### Computational limits:

- Pv3Rs scales to hundreds of markers but not whole-genome sequence (WGS) data.  

- We do not recommend running `compute_posterior()` for data whose total 
genotype count (sum of per-episode multiplicities of infection) exceeds eight. 
If the total genotype counts exceeds eight but there are multiple recurrences,
it might be possible to generate recurrent state estimates for recurrences
analysing episodes pairwise (this approach was used in [2]).

- We have not tested the per-marker allele limit of `compute_posterior()`. 
Very high marker cardinalities could lead to very small allele frequencies and 
thus underflow problems. 


### Population-level allele frequencies: 

In addition to paired data, `compute_posterior()` requires as input
population-level allele frequencies. To minimise bias due to within-host
selection of recrudescent parasites, we recommend using only enrolment episodes
to estimate population-level allele frequencies, and ideally enrolment episodes
from study participants selected at random, not only study participants who
experience recurrence. That said, if most recurrences are either reinfections or
relapses, both of which are draws from the mosquito population (albeit a delayed
draw in the case of a relapse), assuming there is no systematic within-patient
selection (as might occur when infections encounter lingering drug pressure),
estimates based on all episodes should be unbiased and more precise than those
based on enrolment episodes only.

### Read-count data: 

Unfortunately, the Pv3Rs model does not exploit read count data at present.
However, read count data could be used to compute population-level allele
frequencies, assuming they are not biased by experimental artefacts.


## Installation 

```r
# Install or update latest stable version of devtools from CRAN
install.packages("devtools")

# Install R.rsp, a vignette builder used to build LaTeX vignettes in Pv3Rs
install.packages("R.rsp")

# Install Pv3Rs from GitHub 
# I highly recommend doing this in RStudio as RStudio installs pandoc needed to build vignettes.
# If you're working in R outside of RStudio you might need to install pandoc and check its path; 
# otherwise set build_vignettes = FALSE
devtools::install_github("aimeertaylor/Pv3Rs", build_vignettes = TRUE)

# Load and attach Pv3Rs
library(Pv3Rs)

# List links to all available documentation
help(package = "Pv3Rs")

# List links to vignettes
vignette(package = "Pv3Rs")

# View function documentation, e.g., 
?compute_posterior
```
