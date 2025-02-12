% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/st_assistantFunctions.R
\name{convertMatrixMouseGene}
\alias{convertMatrixMouseGene}
\title{convertMatrixMouseGene}
\usage{
convertMatrixMouseGene(counts, host = "https://dec2021.archive.ensembl.org/")
}
\arguments{
\item{counts}{A matrix whose rownames are mouse genes}

\item{host}{Host to connect to in `useMart` function}
}
\description{
Convert mouse genes to human genes in a matrix with mouse genes as rownames
}
