% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataset.building.R
\name{dataset.building}
\alias{dataset.building}
\title{dataset construction for the application of the CCC method}
\usage{
dataset.building(x)
}
\arguments{
\item{x}{must be a data.frame containing the compounds' names ("compound"), the monoisotopic mass ("MM"), the retention time ("RT"), the chemical formula ("Formula") and the smiles ("smiles")}
}
\value{
the basic dataset to build a CCC_method
}
\description{
the function builds the dataset to apply the the CCC_method to an in-house dataset of standards. The new-built must be saved as model
}
\examples{
STD_RP <- read.csv(system.file("extdata", "STD_RP.csv", package = "CCC"), row.names = 1, stringsAsFactors = FALSE)
X1Y <- dataset.building(STD_RP)
}
\author{
Luca Narduzzi "nardluca@gmail.com"
}
