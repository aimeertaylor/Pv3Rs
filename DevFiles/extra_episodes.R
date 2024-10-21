#####
# Toy example involving a pair of episodes that is suggestive of recrudescence,
# along with an extra, uninformative episode with no genetic data.
#
# Posterior  probabilities are noticeably different when the uninformative
# episode occurs before vs after the pair of episodes corresponding to
# recrudescence. In particular, the probability of recrudescence is higher if
# the uninformative episode occurs before the recrudescence pair of episodes.
# The lack of symmetry is primarily due to the fact that out of the graph space
# of MOI = (1, 1, 1), there are 3 graphs that support (relapse, reinfection),
# but 5 graphs that support (reinfection, relapse).
#
# This analysis is also repeated with the recrudescence replaced by reinfection.
#####

library(Pv3Rs)

marker_count <- 100
ms <- paste0("m", 1:marker_count)

match <- as.list(sapply(ms, function(t) "A")) # O
mismatch <- as.list(sapply(ms, function(t) "B")) # X
no_data <- as.list(sapply(ms, function(t) NA)) # N
fs <- sapply(ms, FUN = function(m) c("A" = 0.01, "B" = 0.99), simplify = FALSE)

ys_list <- list("OO" = list(enroll = match,
                            recur1 = match),
                "OON" = list(enroll = match,
                             recur1 = match,
                             recur2 = no_data),
                "NOO" = list(enroll = no_data,
                             recur1 = match,
                             recur2 = match),
                "OX" = list(enroll = match,
                            recur1 = mismatch),
                "OXN" = list(enroll = match,
                             recur1 = mismatch,
                             recur2 = no_data),
                "NOX" = list(enroll = no_data,
                             recur1 = match,
                             recur2 = mismatch))

res_list <- suppressMessages(lapply(ys_list, compute_posterior, fs))

res_list[["OO"]]$marg
res_list[["OON"]]$marg
res_list[["NOO"]]$marg # Pr of C increases if you have a no-data episode prior to matching episodes

res_list[["OON"]]$joint # P(Cx) = 0.255, P(LC) = 0.085, P(LI) = 0.085, P(LL) = 0.064
res_list[["NOO"]]$joint # P(xC) = 0.264, P(CL) = 0.088, P(IL) = 0.053, P(LL) = 0.066

res_list[["OX"]]$marg
res_list[["OXN"]]$marg
res_list[["NOX"]]$marg # Pr of I increases if you have a no-data episode prior to mismatching episodes
