% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{LUNC}
\alias{LUNC}
\title{LUCID Nomenclature (LUNC) data}
\format{
## `LUNC`
A data frame with 35 rows and 3 columns:
\describe{
 \item{fcs_parameters_desc}{character the parameters description we get from fcs}
 \item{source}{character some information on how we got this symbol}
 \item{symbol}{character symbol we've chosen as standard label, this should preferably be from `HCDM` or `HGNC`}
 \item{symbol_database}{character name of the database where we got the symbol from,
 department means that we kind of made up the name ourselves. Note that e.g. CD3 is not in HCDM but CD8A is}
 }
}
\usage{
LUNC
}
\description{
Data frame with nomenclature used at our department. This is basically a list
of parameter descriptions that come from .fcs files that not have a CD label
or HGNC label. We use these data to convert our nomenclature
to (if able) `HCDM` nomenclature, or else to `HGNC`.
}
\keyword{datasets}
