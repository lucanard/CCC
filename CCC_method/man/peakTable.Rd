\name{peakTable}
\alias{peakTable}
\docType{data}
\title{Garniga exudate extract peakTable
}
\description{It contains a peakTable obtained from metaMS. The data has been obtained extracting and injecting the exudate obtained from grass fermentation in a LC-Q-TOF-MS system.
}
\usage{data("peakTable")}
\format{
  A data frame with 2000 observations on the following 31 variables.
  \describe{
    \item{\code{ChemSpiderID}}{a character vector containing the chemspider ID of the identified compounds}
    \item{\code{compound}}{a character vector containing the name of the compounds identified}
    \item{\code{pcgroup}}{a numeric vector containing the CAMERA pseudo-spectra grouping}
    \item{\code{adduct}}{a character vector containing the adducts detected by CAMERA}
    \item{\code{isotopes}}{a character vector containing the isotopic pattern detected by CAMERA}
    \item{\code{mz}}{a numeric vector of the detected m/z detected in the sampleset}
    \item{\code{rt}}{a numeric vector of the detected RT detected in the sampleset}
    \item{\code{A002_S1D2_RP_neg01}}{a numeric vector}
    \item{\code{A002_S1D2_RP_neg02}}{a numeric vector}
    \item{\code{A003_S2D2_RP_neg01}}{a numeric vector}
    \item{\code{A003_S2D2_RP_neg02}}{a numeric vector}
    \item{\code{A004_S3D2_RP_neg01}}{a numeric vector}
    \item{\code{A004_S3D2_RP_neg02}}{a numeric vector}
    \item{\code{A005_S4D2_RP_neg01}}{a numeric vector}
    \item{\code{A005_S4D2_RP_neg02}}{a numeric vector}
    \item{\code{A006_S1D3_RP_neg01}}{a numeric vector}
    \item{\code{A006_S1D3_RP_neg02}}{a numeric vector}
    \item{\code{A007_S2D3_RP_neg01}}{a numeric vector}
    \item{\code{A007_S2D3_RP_neg02}}{a numeric vector}
    \item{\code{A008_S3D3_RP_neg01}}{a numeric vector}
    \item{\code{A008_S3D3_RP_neg02}}{a numeric vector}
    \item{\code{A009_S4D3_RP_neg01}}{a numeric vector}
    \item{\code{A009_S4D3_RP_neg02}}{a numeric vector}
    \item{\code{A010_S1D4_RP_neg01}}{a numeric vector}
    \item{\code{A010_S1D4_RP_neg02}}{a numeric vector}
    \item{\code{A011_S2D4_RP_neg01}}{a numeric vector}
    \item{\code{A011_S2D4_RP_neg02}}{a numeric vector}
    \item{\code{A012_S3D4_RP_neg01}}{a numeric vector}
    \item{\code{A012_S3D4_RP_neg02}}{a numeric vector}
    \item{\code{A013_S4D4_RP_neg01}}{a numeric vector}
    \item{\code{A013_S4D4_RP_neg02}}{a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(peakTable)
}
\keyword{datasets}
