% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/countFunctions.R
\name{extractACANCounts}
\alias{extractACANCounts}
\title{This function calculate AC AN from a gds File stratified in groups}
\usage{
extractACANCounts(
  gds,
  sampleID,
  bedGRange,
  caseGroupInfo,
  groupIDSpecified,
  nVariant,
  overlapType = "within",
  sexID = NULL,
  reference = NULL
)
}
\arguments{
\item{gds}{a gds file handle}

\item{sampleID}{restrict to the samples specified}

\item{bedGRange}{If set, it is a GRanges or GRangesList to specify a region}

\item{caseGroupInfo}{A vector of case group ID}

\item{groupIDSpecified}{A vector of group ID used}

\item{nVariant}{number of variants}

\item{overlapType}{the overlap type between the bedGRange and the variants
in the gdsFile. See type in the findOverlaps function in GenomicRanges}

\item{sexID}{the sex ID, 1 for male, 2 for female, used for counting the
total alleles in chromosome X and Y, default is NULL}

\item{reference}{the reference build, used to identify the PAR region in 
chromosome X and Y, default is NULL}
}
\value{
a list of two components. The first is the AC matrix of different
groups, the second is the AN matrix of different groups
}
\description{
This function calculate AC AN from a gds File stratified in groups
}
