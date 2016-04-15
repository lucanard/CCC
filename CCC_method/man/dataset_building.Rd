\name{dataset.building}
\alias{dataset.building}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Builds the standards dataset for the model testing in the CCC method
%%  ~~function to do ... ~~
}
\description{The function builds the standards dataset in the CCC_code necessary to test the models in the CCC_method
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dataset.building(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{x must be a data.frame containing the compounds' names ("compound"), the monoisotopic mass ("MM"), the retention time ("RT"), the chemical formula ("Formula") and the smiles ("smiles")
  \item{x}{
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{a matrix of 18 columns containing the X-Y variables to build the CCC_model. The Y variables are: "Sulphur", "Phosphorous", "phenolics", "Acidic", "NN", "het", "aliph", "CO", and "bs". The X variables are: Monoisotopic mass ("MM"), the retention time ("RT"), the mass defect ("md"), the relative mass defect ("RMD"), the number of carbons ("nC"), the percentage of carbons ("pC"), the residual relative mass defect ("rRMD"), the presence/absence of odd value for mass ("odd") and the isotopic mass defect and pattern ("IMD").
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Luca Narduzzi nardluca@gmail.com
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
dataset <- dataset.building(STD_RP)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
