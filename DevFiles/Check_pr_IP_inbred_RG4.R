# NOTE: from deprecated model

################################################################################
# This script could be converted into the equivalent of Check_pr_IP_inbred_RG3.R
# It contains probabilities that were validated in IG_RG4_model_checking.R,
# which has been deleted as no-longer runs.
# The derivation of the probabilities can be found in Pr_IG_RG_inbred.pdf
# Example RGs can be found in Graphs_2_2.pdf
# RG names refer to one specific RG graph that type
# IGs are summarised as counts of elements per set, e.g. "31" means three nodes
# are IBD and one is not.
# As such, probabilities sum over different partitions,
# e.g. "31" sums over (1,2,3)(4), (1,2,4)(3), (1,3,4)(2) and (2,3,4)(1) whose
# individual probabilities may differ depending on RG (see Pr_IG_RG_inbred).
################################################################################

# RG names:
# RG_s2_11 refers to a specific RG graph where two nodes are siblings
# RG_s2_s2 refers to a specific RG graph in which there are two pairs of siblings
# etc.
#
# Examples of these graphs can be found in Graphs_2_2.pdf. Note that RG_c3_1 is not available because
# Graphs_2_2.pdf was constructed assuming MOIs of 2 and 2, respectively, and clonal relationships
# within infections are not modelled.
RG_is <- c(RG_s2_11 = 2,
           RG_s2_s2 = 13,
           RG_s4 = 33,
           RG_s3_1 = 22,
           RG_c2_11 = 3,
           RG_c2_c2 = 11,
           #RG_c3_1 = NA, # Not available because within infection clones are indistinguishable
           RG_s2_c2 = 8,
           RG_cs3_1 = 21,
           RG_cs24 = 35,
           RG_cs44 = 37)

Pr_IG_RGs <- cbind(RG_s2_11 = c("1111" = 0.5*(1-f)*(1-2*f)*(1-3*f), # Not written down
                                "211" = 2*f*(1-f)*(1-2*f) + 0.5*(1-f)*(1-2*f)*(1-3*f) + # when IBD and sibling edge align (page 1)
                                        0.5*f*(1-f)*(1-2*f) + # when IBD and sibling edge are opposite one another (page 2)
                                        4*0.5*f*(1-f)*(1-2*f), # when IBD and sibling edges touch one another (page 2)
                                "22" = 0.5*f*(1-f)*(1-2*f) + (5/2)*(1-f)*f^2, # pages 3 and 4
                                "31" = 2*((1-f)*f^2 + 0.5*f*(1-f)*(1-2*f)) +  # when IBD triangle does not encompasses sibling edge (two situations - page 4)
                                  2*(1-f)*f^2, # when IBD triangle encompasses sibling edge (two situations - page 4)
                                "4" = f^3 + 0.5*(1-f)*f^2), # page 4
                   RG_s2_s2 = c("1111" = 0.25*(1-f)*(1-2*f)*(1-3*f), # page 5
                                "211" = 2*f*(1-f)*(1-2*f) + 0.5*(1-f)*(1-2*f)*(1-3*f) + # page 5
                                  f*(1-f)*(1-2*f), # page 6
                                "22" = 0.25*(1-f)*(1-2*f)*(1-3*f) + (1+3/4)*f*(1-f)*(1-2*f) + (2+1/4)*(1-f)*f^2 + # page 6
                                  2*(1/4)*(1-f)*f^2, # page 7
                                "31" = 4*(0.25*f*(1-f)*(1-2*f) + 0.75*(1-f)*f^2), # page 8
                                "4" = f^3 + 1.25*(1-f)*f^2 + 0.25*f*(1-f)*(1-2*f)), # page 7
                   RG_4sib = c("1111" = 0,
                               "211" = 0,
                               "22" = (3/8)*(1-f), # page 9
                               "31" = 0.5*(1-f), # page 9
                               "4" = (1/8)*(1-f) + f), # page 9
                   RG_s3_1 = c("1111" = 0,
                               "211" = 3*(1/4)*(1-f)*(1-2*f),
                               "22" = 3*0.25*f*(1-f),
                               "31" = (1+1/4)*f*(1-f) + (1/4)*(1-f)*(1-2*f) + 3*0.25*f*(1-f),
                               "4" = (1/4)*f*(1-f) + f^2),
                   RG_c2_11 = c("1111" = 0,
                                "211" = (1-f)*(1-2*f),
                                "22" = f*(1-f),
                                "31" = 2*f*(1-f),
                                "4" = f^2),
                   RG_c2_c2 = c("1111" = 0,Pr
                                "211" = 0,
                                "22" = 1-f,
                                "31" = 0,
                                "4" = f),
                   RG_c3_1 = c("1111" = 0,
                               "211" = 0,
                               "22" = 0,
                               "31" = 1-f,
                               "4" = f),
                   RG_s2_c2 = c("1111" = 0,
                                "211" = 0.5*(1-f)*(1-2*f),
                                "22" = 0.5*(1-f)*(1-2*f) + 1.5*f*(1-f),
                                "31" = 2*0.5*f*(1-f),
                                "4" = 0.5*f*(1-f) + f^2),
                   RG_cs3_1 = c("1111" = 0,
                                "211" = 0.5*(1-f)*(1-2*f),
                                "22" = 0.5*f*(1-f),
                                "31" = 0.5*(1-f)*(1-2*f) + 1.5*f*(1-f) + 0.5*f*(1-f),
                                "4" = 0.5*f*(1-f) + f^2),
                   RG_cs24 = c("1111" = 0,
                               "211" = 0,
                               "22" = 0.25*(1-f),
                               "31" = 2*0.25*(1-f),
                               "4" = 0.25*(1-f) + f),
                   RG_cs44 = c("1111" = 0,
                               "211" = 0,
                               "22" = 0.5*(1-f),
                               "31" = 0,
                               "4" = 0.5*(1-f) + f))
