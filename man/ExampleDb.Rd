% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Data.R
\docType{data}
\name{ExampleDb}
\alias{ExampleDb}
\alias{ClonedExampleDb}
\title{Example Change-O database}
\format{A data.frame with the following Change-O style columns:
  \itemize{
    \item  \code{SEQUENCE_ID}:           Sequence identifier
    \item  \code{SEQUENCE_IMGT}:         IMGT-gapped observed sequence.
    \item  \code{GERMLINE_IMGT_D_MASK}:  IMGT-gapped germline sequence with N, P and
                                         D regions masked.
    \item  \code{V_CALL}:                V region allele assignments.
    \item  \code{V_CALL_GENOTYPED}:      TIgGER corrected V region allele assignment.
    \item  \code{D_CALL}:                D region allele assignments.
    \item  \code{J_CALL}:                J region allele assignments.
    \item  \code{JUNCTION}:              Junction region sequence.
    \item  \code{JUNCTION_LENGTH}:       Length of the junction region in nucleotides.
    \item  \code{NP1_LENGTH}:            Combined length of the N and P regions proximal
                                         to the V region.
    \item  \code{NP2_LENGTH}:            Combined length of the N and P regions proximal
                                         to the J region.
    \item  \code{SAMPLE}:                Sample identifier. Time in relation to vaccination.
    \item  \code{ISOTYPE}:               Isotype assignment.
    \item  \code{DUPCOUNT}:              Copy count (number of duplicates) of the sequence.
}}
\usage{
ExampleDb

ClonedExampleDb
}
\description{
A small example database subset from Laserson and Vigneault et al, 2014.

Includes the \code{CLONE} column.
}
\references{
\enumerate{
  \item  Laserson U and Vigneault F, et al. High-resolution antibody dynamics of
           vaccine-induced immune responses.
           Proc Natl Acad Sci USA. 2014 111:4928-33.
}
}
\keyword{datasets}