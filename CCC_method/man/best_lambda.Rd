\name{best_lambda}
\alias{best_lambda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Chooses the best lambda for your Lasso and Ridge models
%%  ~~function to do ... ~~
}
\description{The functions selects the best lambda for your Lasso and/or Ridge model through 1000 iterative training-test set splitting and testing for the best lambda.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
best_lambda(X1Y, ny, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{X1Y <- the dataset given by the function dataset.building
ny <- the Y variable under evaluation (values from 1 to 9) according to the order given by the function dataset.building
alpha <- select 0 for RIDGE and 1 for LASSO model
  \item{x}{
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{It gives a single value to be used for your LASSO and RIDGE model testing
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
lambda <- best_lambda(dataset, ny = 1, alpha = 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
