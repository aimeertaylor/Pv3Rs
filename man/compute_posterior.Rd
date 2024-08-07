% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_posterior.R
\name{compute_posterior}
\alias{compute_posterior}
\title{Compute the posterior distribution of recurrence states}
\usage{
compute_posterior(
  y,
  fs,
  prior = NULL,
  MOIs = NULL,
  return.RG = FALSE,
  return.logp = FALSE
)
}
\arguments{
\item{y}{Observed data in the form of a list of lists. Alleles should be
provided as a set of distinct alleles for each infection (per marker).
The number of entries is the number of episodes in chronological order.
Episode names can be specified, but they are not used. Each episode is in
turn a list of distinct observed alleles for each marker, which must be
named, or \code{NA} if not observed. If a vector of observed alleles
contains \code{NA} for a given marker, the non-\code{NA} entries will be
ignored due to ambiguity. For a given marker, alleles are modeled as
categorical random variables. As such, allele names are arbitrary, but must
correspond with frequency names (see examples below). The same names can be
used for alleles belonging to different markers. As such, frequencies must
be specified per named allele per named marker.}

\item{fs}{List of allele frequencies as vectors. Names of the list must match
with the marker names in \code{y}. Within lists (i.e., for each marker),
frequencies must be specified per allele name.}

\item{prior}{Matrix of prior probabilities of the recurrence states for each
recurrent episode. Each row corresponds to an episode in chronological
order. The column names must be C, L, and I for recrudescence, relapse and
reinfection respectively. Row names can be specified by they are not used.
If prior is not provided, a uniform prior will be used.}

\item{MOIs}{Multiplicity of infection for each episode. If MOIs are not
provided, the most parsimonious MOIs will be used; see
\code{\link{determine_MOIs}}.}

\item{return.RG}{Boolean for whether to return the relationship graphs,
defaults to \code{FALSE}.}

\item{return.logp}{Boolean for whether to return the log-likelihood for each
relationship graph, defaults to \code{FALSE}. Setting this to \code{FALSE} allows
for permutation symmetries to be exploited to save computational time, see
\code{\link{enumerate_alleles}}. Setting this to \code{TRUE} will result in
longer runtimes, especially in the case of a larger multiplicity of
infection.}
}
\value{
List containing:
\describe{
\item{marg}{Matrix of marginal posterior probabilities of the possible
recurrence states for each recurrent episode, one row per reinfection.}
\item{joint}{Vector of joint posterior probabilities of each possible string
of recurrence states.}
\item{RGs}{List of relationship graphs with their log-likelihoods stored.
Only returned if \code{return.RG} is \code{TRUE}. See
\code{\link{enumerate_RGs}}.}
}
}
\description{
Entry point to Bayesian inference for \emph{P. vivax} recurrence states based
on genetic data. Specifically, this function finds the posterior
probabilities of relapse, reinfection and recrudescence from genetic data.
Please note that the progress bar does not necessarily increment at a uniform
rate, and may sometimes appear to be stuck while the code is still running.
}
\details{
We enumerate all possible relationship graphs between genotypes, where each
pair of genotypes may be clones, siblings, or strangers, each with a
different level of expected genetic relatedness. The likelihood of a
sequence of recurrence states can be determined from the likelihood of all
relationship graphs compatible with said sequence. More details on the
enumeration and likelihood calculation of relationship graphs can be found in
\code{\link{enumerate_RGs}} and \code{\link{RG_inference}} respectively.

Model assumptions:
\itemize{
\item{No within-host mutations, genotyping errors, undetected alleles}
\item{Parasites are outbred}
\item{Relationship graphs are equally likely given recurrence states}
}

\code{compute_posterior} supports missing data, which should be encoded as
\code{NA}s. However, to avoid estimates based entirely on the weakly
informative nature of multiple per-marker allele calls (see example below and
\code{vignette("missing_data", "Pv3Rs")} for more details), we recommend against
generating estimates for recurrences that have no paired data due to
missingness (see Microsatellite data example in \code{vignette("demo", "Pv3Rs")}).

The data input expects each list of alleles (for an infection) to consist of
a set of distinct alleles. Note that \code{Pv3Rs} supports prevalence data
but not quantitative (proportional abundance) data. If repeated alleles are
provided, they will be collapsed to one occurrence.
}
\examples{

# ===========================================================================
# Example where alleles are named arbitrarily
# ===========================================================================
# Data on an enrolment episode and a recurrence:
y <- list(episode0 = list(marker1 = c("Tinky Winky", "Dipsy"),
                          marker2 = c("Laa-Laa", "Po")),
          episode1 = list(marker1 = "Tinky Winky",
                          marker2 = "Laa-Laa"))

# Allele frequencies:
fs <- list(
  marker1 = setNames(c(0.4, 0.6), c("Tinky Winky", "Dipsy")),
  marker2 = setNames(c(0.2, 0.8), c("Laa-Laa", "Po"))
)

# Compute posterior probabilities using default uniform prior, note that
# since there is only one recurrence, the marginal probabilities are the same
# as the joint probabilities:
compute_posterior(y, fs)



# ===========================================================================
# Example where alleles are given numeric names: 1, 2, 3, 4, 5
# ===========================================================================
# Allele frequencies with numeric names
fs <- list(m1=c('1'=0.78, '2'=0.14, '3'=0.07, '4'=0.005, '5' = 0.005),
           m2=c('1'=0.27, '2'=0.35, '3'=0.38),
           m3=c('1'=0.55, '2'=0.45))

# Data:
y <- list(enrol = list(m1=c('3','2'), m2=c('1','3'), m3=c('1','2')),
          recur1 = list(m1=c('1','4'), m2=c('1','2'), m3=c('1','2')),
          recur2 = list(m1=c('1','5'), m2=c('2','3'), m3=c('1')))

compute_posterior(y, fs)



# ===========================================================================
# Example demonstrating the cosmetic-only nature of episode names. Due to the
# cosmetic-only nature of episode names, input info (episode data and prior)
# and output results should be chronologically ordered and interpreted.
# ===========================================================================
y <- list(enrol = list(m1 = NA),
          recur2 = list(m1 = NA),
          recur1 = list(m1 = NA))
prior <- array(c(0.2,0.7,0.2,0.3,0.6,0), dim = c(2,3),
               dimnames = list(c("recur1", "recur2"), c("C", "L", "I")))

# The first prior row named "recur1" is returned for the first recurrence
# despite it being named "recur2" and the the second prior row named "recur2"
# is returned for the second recurrence despite it being named "recur1":
compute_posterior(y, fs, prior)
prior



#============================================================================
# compute_posterior() returns the prior when there are no data. However, the
# prior will be re-weighted if NAs encode MOIs that are incompatible with
# recrudescence. (Recrudescing parasites are clones of parasites in the
# preceding blood-stage infection. The Pv3R model assumes no within-host
# mutations, genotyping errors or undetected alleles. As such, recrudescence
# is incompatible with an MOI increase relative to the preceding infection.)
#============================================================================
# Allele frequencies:
fs <- list(m1 = setNames(c(0.25, 0.75), c("A", "Other")))

# Data on an enrolment episode and two recurrences:
y_no_data_MOIs111 <- list(enrol = list(m1 = NA),
                          recur1 = list(m1 = NA),
                          recur2 = list(m1 = NA))

y_no_data_MOIs121 <- list(enrol = list(m1 = NA),
                          recur1 = list(m1 = c(NA, NA)),
                          recur2 = list(m1 = NA))

# Compute posterior probabilities:
post3Rs_no_data_MOIs111 <- compute_posterior(y_no_data_MOIs111, fs)
post3Rs_no_data_MOIs121 <- compute_posterior(y_no_data_MOIs121, fs)

# Returns the default prior since MOIs 111 are compatible with all 3R
# sequences:
post3Rs_no_data_MOIs111

# Returns the default prior re-weighted to the exclusion of 3R sequences with
# a recrudescence at the first recurrence:
post3Rs_no_data_MOIs121



#============================================================================
# Example demonstrating the weakly informative nature of multiple per-marker
# allele calls; see vignette("missing_data", "Pv3Rs") for more details.
#============================================================================
fs = list(m1 = c('1' = 0.5, '2' = 0.5))
y <- list(enrol = list(m1 = c('1', '2')), recur = list(m1 = NA))

# compute_posterior() does not return the prior despite there being no
# recurrent data:
compute_posterior(y, fs)



#============================================================================
# Example of the small but undesirable effect on the posterior of prior on
# graphs: the marginal probability that the first recurrence is a
# recrudescence increases as the number of recurrences increases even though
# only the first recurrence has data (also see vignette [to-do - base on
# MyDevFiles/Graph_prior_bias_examples.R])
#============================================================================
# Allele frequencies:
fs <- list(m1 = setNames(c(0.25, 1-0.25), c("A", "Other")))

# Data for different scenarios; scenarios where the number of recurrences
# increases but only the first recurrence has data
ys <- list(y1 = list(enrol = list(m1 = "A"), recur1 = list(m1 = "A")),
           y2 = list(enrol = list(m1 = "A"), recur1 = list(m1 = "A"), recur2 = list(m1 = NA)),
           y3 = list(enrol = list(m1 = "A"), recur1 = list(m1 = "A"), recur2 = list(m1 = NA), recur3 = list(m1 = NA)),
           y4 = list(enrol = list(m1 = "A"), recur1 = list(m1 = "A"), recur2 = list(m1 = NA), recur3 = list(m1 = NA), recur4 = list(m1 = NA)))

# Compute posterior probabilities and extract marginal probabilities:
results <- lapply(ys, function(y) compute_posterior(y, fs)$marg)

# Extract results for the first recurrence only:
first_recur <- sapply(results, function(result) result[1,])

# Plot 2D simplex
n_recur <- max(sapply(ys, length)-1)
pardefault <- par()
par(mar = c(0,0,0,0))
plot_simplex(v_labels = rownames(first_recur))

# Project probabilities onto 2D simplex coordinates
xy <- apply(first_recur, 2, project2D)

# Plot divergence from one recurrence to four:
arrows(x0 = xy["x", 1], x1 = xy["x", n_recur],
       y0 = xy["y", 1], y1 = xy["y", n_recur],
       length = 0.05, col = "red")

# Plot a point for each recurrence from one to four:
points(x = xy["x", ], y = xy["y", ], pch = ".")

# Restore plotting margins
par(mar = pardefault$mar)

}
