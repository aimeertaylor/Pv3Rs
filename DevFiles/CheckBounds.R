# TO MOVE TO CHECK
ms <- paste0("m", 1:100)
fA <- 0.14
fm <- c("A" = fA, "B" = 1-fA)
fs <- sapply(ms, function(m) fm, simplify = F)

# For seqs CC, CI, CL, LC, LI, LL when MOIs are 1, 1, 1,
# there are 1, 1,  3,  3,  3,  12 compatible graphs

y = list(enroll = all_As, recur1 = all_As, recur2 = no_data)
post <- sm(compute_posterior(y, MOIs = c(2,1,1), fs = fs))
post$joint # CC, CL, CI are equal
c("Estimate" =  post$marg["recur1","C"],
  "Maximum" = 3 / (3 + 1/3 + 1/3 + 3/12)) # Relies on knowledge of the data

# CC should be at its maximum
post <- sm(compute_posterior(y = list(enroll = all_As,
                                      recur1 = all_As,
                                      recur2 = all_As), fs = fs))
c("Estimate" = as.numeric(post$joint["CC"]),
  "Maximum" = 1 / (1 + 1/3 + 1/3 + 1/12)) # Not reliant on knowledge of data

# CL should be at its maximum
post <- sm(compute_posterior(y = list(enroll = all_As,
                                      recur1 = all_As,
                                      recur2 = c(all_As[1:50], all_Bs[51:100])),
                             fs = fs))
c("Estimate" = as.numeric(post$joint["CL"]), # Not reliant on knowledge of data
  "Maximum" = 1 / (1 + 3/12))

