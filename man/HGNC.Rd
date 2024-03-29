% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{HGNC}
\alias{HGNC}
\title{HUGO Gene Nomenclature Committee (HGNC) data}
\format{
## `HGNC`
A data frame with 43667 rows and 2 columns:
\describe{
 \item{hgnc_id}{character HGNC ID, a unique ID created by the HGNC for every approved symbol}
 \item{symbol}{character The HGNC approved gene symbol}
 \item{cd}{character Symbol used within the Human Cell Differentiation Molecule (HCDM) database}
 }
}
\source{
<https://www.genenames.org/download/archive/>
}
\usage{
HGNC
}
\description{
Data frame that contains three columns from the current (August 2023)
hgnc_complete_set database. These data are the second reference for
naming the cellular markers
}
\keyword{datasets}
