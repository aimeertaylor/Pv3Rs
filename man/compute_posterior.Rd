% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_posterior.R
\name{compute_posterior}
\alias{compute_posterior}
\title{Compute posterior probabilities of \emph{P. vivax} recurrence states}
\usage{
compute_posterior(
  y,
  fs,
  prior = NULL,
  MOIs = NULL,
  return.RG = FALSE,
  return.logp = FALSE,
  progress.bar = TRUE
)
}
\arguments{
\item{y}{List of lists encoding allelic data. The outer list contains
episodes in chronological order. The inner list contains named
markers per episode. Marker names must be consistent across episodes. \code{NA}
indicates missing marker data; otherwise, specify a per-marker
vector of distinct alleles detected (presently, \code{compute_posterior()} does
not support data on the proportional abundance of detected alleles).
Repeat alleles and \code{NA} entries within allelic vectors are ignored.
Allele names are arbitrary, allowing for different data types, but must
correspond with frequency names.}

\item{fs}{List of per-marker allele frequency vectors, with names matching
marker names in \code{y}. Per-marker alleles frequencies mut contain one
frequency per named allele, with names matching alleles in \code{y}. Per-marker
frequencies must sum to one.}

\item{prior}{Matrix of prior probabilities of recurrence states per episode,
with rows as episodes in chronological order, and columns named "C", "L",
and "I" for recrudescence, relapse and reinfection, respectively. Row names
are ignored. If \code{NULL} (default), per-episode recurrence states are assumed
equally likely \emph{a priori}.}

\item{MOIs}{Vector of per-episode multiplicities of infection (MOIs); because
the Pv3Rs model assumes no genotyping errors, \code{MOIs} must be greater than
or equal to the most parsimonious MOI estimates compatible with the data;
see \code{determine_MOIs(y)}. These are the estimates used when
\code{MOIs = NULL} (default).}

\item{return.RG}{Logical; returns the relationship graphs
(default \code{FALSE}). Automatically set to \code{TRUE} if \code{return.logp = TRUE}.}

\item{return.logp}{Logical; returns the log-likelihood for each relationship
graph (default \code{FALSE}). Setting \code{TRUE} disables permutation symmetry
optimisation and thus increases runtime, especially when MOIs are large.
Does not affect the output of the posterior probabilities; for an an
example of permutation symmetry, see \strong{Exploration of relationship graphs} in
\HTMLVignette{demonstrate-usage}{exploration-of-relationship-graphs}{Demonstrate Pv3Rs usage}.}

\item{progress.bar}{Logical; show progress bars (default \code{TRUE}).
Note that the progress bar may update non-uniformly.}
}
\value{
List containing:
\describe{
\item{\code{marg}}{Matrix of marginal posterior probabilities for each recurrence,
with rows as recurrences and columns as "C" (recrudescence), "L"
(relapse), and "I" (reinfection). Each marginal probability sums over a
subset of joint probabilities. For example, the marginal probability of
"C" at the first of two recurrences sums over the joint probabilities
of "CC", "CL", and "CI".}
\item{\code{joint}}{Vector of joint posterior probabilities for each recurrence
state sequence; within a sequence "C", "L", and "I" are used as above.}
\item{\code{RGs}}{List of lists encoding relationship graphs; returned only if
\code{return.RG = TRUE} (default \code{FALSE}), and with log-likelihoods if
\code{return.logp = TRUE} (default \code{FALSE}). A relationship graph
encoded as a list can be converted into a \code{igraph} object using
\code{\link{RG_to_igraph}} and thus plotted using
\code{\link{plot_RG}}. For more details on relationship graphs, see
\code{\link{enumerate_RGs}}.}
}
}
\description{
Computes per-person posterior probabilities of \emph{P. vivax} recurrence states —
recrudescence, relapse, reinfection — using per-person genetic data on two
or more episodes. For usage, see \strong{Examples} below and
\HTMLVignette{demonstrate-usage}{}{Demonstrate Pv3Rs usage}.
For a more complete understanding of \code{compute_posterior} output, see
\href{https://aimeertaylor.github.io/Pv3Rs/articles/posterior-probabilities.html}{Understand posterior probabilities}.

Note: The progress bar may increment non-uniformly (see
\strong{Details}); it may appear stuck when computations are ongoing.
}
\details{
\code{compute_posterior()} computes posterior probabilities proportional to
the likelihood multiplied by the prior. The likelihood sums over:
\itemize{
\item ways to phase allelic data onto haploid genotypes
\item graphs of relationships between haploid genotypes
\item ways to partition alleles into clusters of identity-by-descent
}

We enumerate all possible relationship graphs between haploid genotypes,
where pairs of genotypes can either be clones, siblings, or strangers. The
likelihood of a sequence of recurrence states can be determined from the
likelihood of all relationship graphs compatible with said sequence. More
details on the enumeration of relationship graphs can be found in
\code{\link{enumerate_RGs}}. For each relationship graph, the model sums over
all possible identity-by-descent partitions. Because some graphs are
compatible with more partitions than others, the log p(Y|RG) progress bar may
advance non-uniformly. We do not recommend running `compute_posterior() when
the total genotype count (sum of MOIs) exceeds eight because there are too
many relationship graphs.

Notable model assumptions and limitations:
\itemize{
\item{All siblings are regular siblings}
\item{Recrudescent parasites derive only from the immediately preceding
episode}
\item{Recrudescence, relapse and reinfection are mutually exclusive}
\item{Undetected alleles, genotyping errors, and \emph{de novo} mutations
are not modelled}
\item{Population structure and various other complexities that confound
molecular correction are not modelled}
}
}
\examples{
# Numerically named alleles
y <- list(enrol = list(m1 = c('3','2'), m2 = c('1','2')),
          recur1 = list(m1 = c('1','4'), m2 = c('1')),
          recur2 = list(m1 = c('1'), m2 = NA))
fs <- list(m1 = c('1' = 0.78, '2' = 0.14, '3' = 0.07, '4' = 0.01),
           m2 = c('1' = 0.27, '2' = 0.73))
compute_posterior(y, fs, progress.bar = FALSE)


# Arbitrarily named alleles, plotting per-recurrence posteriors
y <- list(enrolment = list(marker1 = c("Tinky Winky", "Dipsy"),
                          marker2 = c("Tinky Winky", "Laa-Laa", "Po")),
          recurrence = list(marker1 = "Tinky Winky",
                          marker2 = "Laa-Laa"))
fs <- list(marker1 = c("Tinky Winky" = 0.4, "Dipsy" = 0.6),
           marker2 = c("Tinky Winky" = 0.1, "Laa-Laa" = 0.1, "Po" = 0.8))
plot_simplex(p.coords = compute_posterior(y, fs, progress.bar = FALSE)$marg)


# Episode names are cosmetic: "r1_prior" is returned for "r2"
y <- list(enrol = list(m1 = NA), r2 = list(m1 = NA), r1 = list(m1 = NA))
prior <- matrix(c(0.6,0.7,0.2,0.3,0.2,0), ncol = 3,
                dimnames = list(c("r1_prior", "r2_prior"), c("C", "L", "I")))
suppressMessages(compute_posterior(y, fs = list(m1 = c(a = 1)), prior))$marg
prior


# Prior is returned when all data are missing
y_missing <- list(enrol = list(m1 = NA), recur = list(m1 = NA))
suppressMessages(compute_posterior(y_missing, fs = list(m1 = c("A" = 1))))


# Return of the prior re-weighted to the exclusion of recrudescence:
suppressMessages(compute_posterior(y_missing, fs = list(m1 = c("A" = 1)),
                 MOIs = c(1,2)))
# (Recrudescing parasites are clones of previous blood-stage parasites. The
# Pv3R model assumes no within-host de-novo mutations and perfect allele
# detection. As such, recrudescence is incompatible with an MOI increase on
# the preceding infection.)


# Beware provision of unpaired data: the prior is not necessarily returned;
# for more details, see link above to "Understand posterior estimates"
y <- list(list(m1 = c('1', '2')), list(m1 = NA))
fs <- list(m1 = c('1' = 0.5, '2' = 0.5))
suppressMessages(compute_posterior(y, fs))$marg

}
