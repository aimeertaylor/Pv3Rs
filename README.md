# Pv3Rs

An R package developed to enable statistical genetic inference of 
*Plasmodium vivax* 

- Relapse
- Reinfection
- Recrudescence

The main function `Pv3Rs::compute_posterior` can be used to estimate posterior 
probabibilities of relapse, reinfection and recrudescence for one or more 
recurrence using genetic data to update a prior, which is ideally informative 
(e.g., based on time-to-event information). 

Other features are more general:

- `Pv3Rs::plot_data` can be used to visualise genetic data for molecular
correction regardless of the analytical method (e.g., it could be used to
visualise *Plasmodium falciparum* data intended for analysis using a WHO match
counting algorithm).

- Together, `Pv3Rs::plot_simplex` and `Pv3Rs::project2D` can be used to plot a 
vector of any three numbers between 0 and 1 that sum to one.

## Please be aware of the following points!

### Unstable state: 

The R package is currently under development and thus liable to contain many
errors. The model around which the Pv3Rs R package is being developed is
documented in the preprint [1]. The model builds upon a prototype documented in
[2] (see [1] to better understand how the model underpinning Pv3Rs differs from
the prototype).

[1] [Taylor, Foo & White, 2022](https://www.medrxiv.org/content/10.1101/2022.11.23.22282669v1)

[2] [Taylor & Watson et al. 2019](https://www.nature.com/articles/s41467-019-13412-x)

### Prior considerations: 

- Genetic data are modelled using a Bayesian model, whose prior is ideally
informative (in [2] prior estimates were generated by a time-to-event model
built by James Watson) because the cause of recurrent *P. vivax* malaria is not
always identifiable from genetic data alone (when the data are consistent with
recurrent parasites that are relatively unrelated to those in all preceding
infections, both reinfection and relapse are plausible; meanwhile, when the data
are compatible with recurrent parasites that are clones of those in the
preceding infection, both recrudescence and relapse are plausible).

- The main Pv3Rs function, `compute_posterior()`, could be used to estimate the
probable cause of recurrent *P. falciparum* malaria by setting the prior
probability of relapse to zero. However, the current version of Pv3Rs is 
suboptimal for *P. falciparum* recurrent state inference: genotyping errors, which
are not accounted for under the current Pv3Rs model, are liable to lead to the
misclassification of recrudescence as reinfection when the probability of relapse
is zero *a priori*.

### Notable assumptions: 

1) Relationship graphs compatible with a given recurrent state sequence are
equally likely *a priori*
2) Perfect detection of alleles (no genotyping error)
3) Perfect detection of parasites (problematic for low density clones)
4) Perfect detection of episodes (requires active follow up) 
5) Mutually exclusive recurrent states (requires frequent follow up with treatment)
6) No within-host *de novo* mutations 
7) Parasites are outbred (implies the parasite population is infinitely large and panmitic)
8) Genetic markers are independent
    - Ignores long-range linkage disequilibrium due to background relatedness
    - Ignores short-range linkage disequilbrium between between sibling parasites 
10) All siblings are regular siblings
    - Siblings are transitive (not true of some parent child-like sibling trios)
    - Siblings are independent (not true of meiotic siblings)
    - Siblings draw from at most two parents (not true of half-siblings)

The first assumption listed above has a small artefactual effect on posterior 
estimates when relationship graphs grow in size. This is demonstrated in one of 
the examples of `Pv3Rs::compute_posterior()` and will be explained in more 
detail in an upcoming vignette. 

For studies with possibly high rates of recrudescence, we recommend a
sensitivity analysis to explore the impact of genotyping errors on
recrudescence. An example will be provided in an upcoming vignette. 

The marker independence assumption is liable to lead to overconfident posterior 
probabilities when markers are linked. 

### Computational limits:

- Pv3Rs scales to 100s of markers but not whole-genome sequence (WGS) data.  

- We do not recommend running `Pv3Rs::compute_posterior()` for data whose total 
genotype count (sum of per-episode multiplicities of infection) exceeds eight.
If the total genotype counts exceeds eight but there are multiple recurrences,
it might be possible to generate recurrent state estimates for recurrences
one-by-one (this approach was used in [2]).

- We have not tested the per-marker allele limit of `Pv3Rs::compute_posterior()`. 
Very high marker cardinalities could lead to very small allele frequencies and 
thus underflow problems. 


### Population-level allele frequencies 

`compute_posterior()` requires population-level alelle frequencies. To avoid
bias due to within-host selection of recrudescent parasites, we recommended
using only enrolment episodes to estimate population-level allele frequencies.
That said, if most recurrences are either reinfections or relapses, both of
which are draws from the mosquito population (albeit a delayed draw in the
case of a relapse), in the absence of systematic within-patient selection (as
might occur when break-through infections encounter lingering drug pressure),
estimates based on all episodes should be unbiased and more precise than those
based on enrolment episodes only.

### Other points: 

Unfortunately, the Pv3Rs model does not exploit read count data at present.
However, read count data could be used to compute population-level allele
frequencies, assuming they are not biased by experimental artefacts.


## Installation 

```r
# Install or update latest stable version of devtools from CRAN
# I highly recommend doing this if you've recently updated R and RStudio to versions 4.3.0 and 2023.3.1.446, respectively; 
# otherwise, you might encounter problems rendering documentation
install.packages("devtools")

# Install paneljudge from GitHub 
# I highly recommend doing this in RStudio as RStudio installs pandoc needed to build vignettes
# If you're working in R outside of RStudio you might need to install pandoc and check its path; 
# otherwise set build_vignettes = FALSE
devtools::install_github("aimeertaylor/Pv3Rs", build_vignettes = TRUE)

# Load and attach the package
library(Pv3Rs)

# View documentation and examples for main function
?compute_posterior

# Load the demo vignette
vignette("demo", package = "Pv3Rs")

# Lists available functions, as well as example data sets and their documentation [check]
help(package = "Pv3Rs")
```
