# Pv3Rs

An R package developed to support model-based classification, using genetic data, of the cause of recurrent *Plasmodium vivax* malaria, i.e. the Pv3Rs:  

- Relapse 
- Reinfection
- Recrudescence

The R package is currently under development; however, a model around which it is being developed is already 
published: [Taylor & Watson et al. 2019](https://www.nature.com/articles/s41467-019-13412-x). In fact, we might re-name this R-package and use it as a basis for more streamlined versions of Pv3R and Pf2R.

Genetic data are modelled using a Bayesian model, whose prior is ideally informative, because the cause of recurrent *P. vivax* malaria is not always 
identifiable from genetic data alone (when the data suggest that recurrent parasites are relatively unrelated to those in all preceding infections, 
both reinfection and relapse are plausible; meanwhile, when the data suggest that recurrent parasites are clones of those in the 
preceding infection, both recrudescence and relapse are plausible). 

## Installation 

```r
# Install or update latest stable version of devtools from CRAN
install.packages("devtools")

# Install paneljudge from GitHub 
devtools::install_github("aimeertaylor/Pv3Rs", build_vignettes = TRUE)

# Load and attach the package
library(Pv3Rs)

# View documentation inc. examples for main function
?compute_posterior

# Load the demo vignette
vignette("demo", package = "Pv3Rs")

# Lists available functions, as well as example data sets and their documentation [check]
help(package = "Pv3Rs")
```
