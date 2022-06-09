# Pv3Rs

An R package developed to support model-based classification of the cause of recurrent *Plasmodium vivax* malaria using genetic data, i.e. the Pv3Rs:  

- Relapse 
- Recrudescence
- Reinfection

The R package is currently under development; however, the model around which it is being developed is already 
published: [Taylor & Watson et al. 2019](https://www.nature.com/articles/s41467-019-13412-x). 

Genetic data are modelled using a Bayesian model, whose prior is ideally informative, because the cause of recurrent *P. vivax* malaria is not always 
identifiable from genetic data alone (when the data suggest that recurrent parasites are unrelated to those in the initial infection, 
both reinfection and relapse are plausible; meanwhile, when the data suggest that recurrent parasites are clones of those in the 
initial infection, both recrudescence and relapse are plausible). 
