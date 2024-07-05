################################################################################
# Model behaviour given repeat values It is possible to encode elevated MOIs
# using repeat values (explore further) It is not possible to encode
# within-sample allele values using repeat values.
#
# To-do: check readme, check analytical solution for MOI encoding (ask YS?)
################################################################################
rm(list = ls())
library(Pv3Rs)

# =====================================================================
# Some allele frequencies
fs <- list(m1=c('1'=0.78, '2'=0.14, '3'=0.07),
           m2=c('1'=0.78, '2'=0.14, '3'=0.07),
           m3=c('1'=0.78, '2'=0.14, '3'=0.07))


# =====================================================================
# Repeat values (inc. NAs) encode MOIs:
y0 <- list(enrol = list(m1=c('1')),
           recur = list(m1=c('1')))

y1 <- list(enrol = list(m1=c('1')),
           recur = list(m1=c('1','1')))

y2 <- list(enrol = list(m1=c('1'), m2=NA),
           recur = list(m1=c('1'), m2=c(NA,NA)))

determine_MOIs(y0)
determine_MOIs(y1)
determine_MOIs(y2)

suppressMessages(compute_posterior(y0, fs))$marg
suppressMessages(compute_posterior(y1, fs))$marg
suppressMessages(compute_posterior(y2, fs))$marg


# =====================================================================
# Repeat observations at homoallelic markers make no difference
# when the MOI is otherwise encoded:
y0 <- list(enrol = list(m1=c('1'), m2=c('1','2')),
           recur = list(m1=c('1'), m2=c('1')))

y1 <- list(enrol = list(m1=c('1','1'), m2=c('1','2')),
           recur = list(m1=c('1'), m2=c('1')))

determine_MOIs(y0)
determine_MOIs(y1)

suppressMessages(compute_posterior(y0, fs))$marg
suppressMessages(compute_posterior(y1, fs))$marg

# =====================================================================
# Repeat observations at heteroallelic markers do make a difference. The
# temptation to encode imbalanced within-sample allele frequencies is
# tantalizing but...

# Most fundamentally, the imbalance that can be enbcoded is constrained by the
# MOI. That is to say, only imbalances that can be represented within MOI
# repeats can be encoded (otherwise the imbalance inf inflates the MOI - see
# below) only imbalances in the rare; For example, the following encodes MOIs of
# 2 and 3; not MOIs of 2 and 2 with with an imbalanced recurrence.
y <- list(enrol = list(m1=c('1','2'), m2=c('1')),
          recur = list(m1=c('1','2','2'), m2=c('1')))
determine_MOIs(y)
suppressMessages(compute_posterior(y, fs))$marg

# There are technical issues with the ordering of markers:
y0 <- list(enrol = list(m1=c('1','1','2'), m2=c('1','2','3')), # runs
           recur = list(m1=c('1'), m2=c('1')))
y1 <- list(enrol = list(m1=c('1','2','3'), m2=c('1','1','2')), # throws an error
           recur = list(m1=c('1'), m2=c('1')))

determine_MOIs(y0)
determine_MOIs(y1)
try(suppressMessages(compute_posterior(y0, fs))$marg)
try(suppressMessages(compute_posterior(y1, fs))$marg)

# There are technical issues when multiple markers have different repeat value
# combinations:
y2 <- list(enrol = list(m1=c('1','1','1'), m2=c('1','1','2'), m3=c('1','2','3')),
           recur = list(m1='1', m2='1', m3=NA))
determine_MOIs(y2)
try(suppressMessages(compute_posterior(y2, fs))$marg)



