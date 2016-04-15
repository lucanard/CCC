\name{apply.model}
\alias{apply.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A wrapper function to apply the CCC_model to xsAnnotate CAMERA objects
%%  ~~function to do ... ~~
}
\description{the function gets the outcome from CAMERA, groups the isotopic patterns, groups the various isotopic patterns in pseudospectras, gives an estimation of their amount of Carbons and builds the CCC_model on each of the isotopic pattern.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
apply.model(xsaFA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{xsaFA = an xsAnnotate object from CAMERA
  \item{x}{
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{a data.frame grouped in spectras, each containg the isotopic clusters, their estimated amount of Carbons and the CCC_model on each isotopic cluster.
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
compi <- apply.model(xsaFA)


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
