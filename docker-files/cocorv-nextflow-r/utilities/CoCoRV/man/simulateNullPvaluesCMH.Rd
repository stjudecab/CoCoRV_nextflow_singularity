% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulateNullPvalues.R
\name{simulateNullPvaluesCMH}
\alias{simulateNullPvaluesCMH}
\title{simulated null P values for CMH exact test using the CDF of the 
p-values under the null directly. This is much faster than simulating the 
2x2xK tables as a intermediate step}
\usage{
simulateNullPvaluesCMH(
  data,
  alternative = "two.sided",
  nReplication = 1000,
  seed = NULL
)
}
\arguments{
\item{data}{a  matrix or data frame with the number of columns being 4xK:
The columns of each 2x2 contingency table should be arranged as follows
   counts with rare mutations in case, 
   counts without rare mutations in case, 
   counts with rare mutations in controls
   counts without rare mutations in controls}

\item{alternative}{see alternative in fisher.test}

\item{nReplication}{the number of replications of empirical null P values}

\item{seed}{set the random seed}
}
\value{
A list of the following components. 
pvalues: observed p values
nullP: the matrix of null P values from the simulated counts
}
\description{
simulated null P values for CMH exact test using the CDF of the 
p-values under the null directly. This is much faster than simulating the 
2x2xK tables as a intermediate step
}
