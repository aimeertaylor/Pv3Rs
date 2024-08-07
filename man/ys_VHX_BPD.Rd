% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ys_VHX_BPD}
\alias{ys_VHX_BPD}
\title{Example Plasmodium vivax data}
\format{
A list of 217 patients; for each patient, a list of one or more
episodes; for each episode, a list of three or more microsatellite markers;
for each marker, a vector of observed alleles (repeat lengths). For example:
\describe{
\item{BPD_103}{Patient identifier: patient 103 in the BPD trial}
\item{BPD_103_1}{Episode identifier: episode one of patient 103 in the BPD trial}
\item{PV.3.27}{Marker identifier: Plasmodium vivax 3.27}
\item{18}{Repeat length: 18}
}
}
\source{
\itemize{
\item \url{https://github.com/jwatowatson/RecurrentVivax/blob/master/RData/GeneticModel/MS_data_PooledAnalysis.RData}
\item \url{https://github.com/aimeertaylor/Pv3Rs/blob/main/data-raw/ys_VHX_BPD.R}
}
}
\usage{
ys_VHX_BPD
}
\description{
Previously-published microsatellite data on Plasmodium vivax parasites
extracted from patients enrolled in the Best Primaquine Dose (BPD) and Vivax
History (VHX) trials; see
\url{https://www.nature.com/articles/s41467-019-13412-x}.
}
\keyword{datasets}
