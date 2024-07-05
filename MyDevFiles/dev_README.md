## Current files

Verify the examples in https://doi.org/10.1101/2022.11.23.22282669
- `Description_examples.R`

Example for Fig1 of the MedRxiv 2022 preprint: 
- `Example_PreprintFig.R`

Examples of unwanted effects due to uniform distribution over graphs
- `UnderstandingGraphBias.R`

Illustrating how bulk data from meiotic sibs is identical to bulk data from their parents:
- `UnderstandingMeioticSibBulkData.R'

Exploring repeat values and how they could encode MOIs (possibly) and intra-sample allele frequency imbalance (to-avoid):
- `UnderstandingRepeatValues.R'

Example of half siblings (not yet described)
- `Half_siblings_YF.R` (`Half_siblings.R` renamed)
- `Half_siblings_AT.R` (Aimee's attempt to understand how the frequency of a relapse classification tends to one with the number of markers despite the inter-to-intra match ratio converging to some fixed value with the number of markers). 

Study of model behaviour given different relationships of recurrent parasites (not yet described)
RelationshipStudy/

Benchmarking study
Benchmark/

## Files from developing current implementation (removed)

Demo a synthetic example with MOI = 8, now in vignettes/demo.Rmd
- `Inference_example_ICL_3212.R`
- `Sample_example_ICL_3212.R`

Demo a microsatellite example, now in vignettes/demo.Rmd
- `MS_data_demo.R`

Compare current implementation against previous implementation
- `Verify_enumerate_IP_RG.R`

## Files from previous implementation (to remove)

- `Check_match_freq_against_prob.R`
- `Check_meiotic_sib_relatedness.R`
- `Check_pr_IP_inbred_RG3.R`
- `Check_pr_IP_inbred_RG4.R`
- `Check_pr_IP_inbred_allstranger.R`
- `Check_pr_IP_outbredRG.R`
- `Graphs_2_2.pdf`
- `Pr_IG_RG_inbred.pdf`