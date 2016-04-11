\name{permutize}
\alias{permutize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{performms permutation of a numeric vector
%%  ~~function to do ... ~~
}
\description{When per.test = TRUE, permutes the numeric vector, to generate a permutation test. per.test default is FALSE.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
permutize(x, per.test=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{x <- the numeric vector to be permuted
per.test <- if TRUE, it permutes the Y variable. default is FALSE.
  \item{x}{
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A numeric vector.
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
Yvar <- permutize(as.numeric(dataset[,1]))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
