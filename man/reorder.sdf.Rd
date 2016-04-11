\name{reorder.sdf}
\alias{reorder.sdf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{reordering function to get the best candidate according to the CCC_method from the sdf files given by MetFrag
%%  ~~function to do ... ~~
}
\description{The function reorders the candidates given by Metfrag according to the values predicted by the CCC method. it handles all the .sdf files given in the directory folder
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
reorder.sdf(directory, compi, selection)
}
%- maybe also 'usage' for other objects documented here.
\arguments{directory <- the directory where the sdf files from metfrag are stored
compi <- the data.frame given by the apply.model function
selection <- the selected mzs by the statistical analysis performed by the user. The values contained here should correspond to the one interesting from the compi data.frame. 
  \item{x}{
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{list of SDFset files similar to the one given by metFrag. Their elements reordered according to the CCC_method
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
candidates <- reorder.csv("H:/R/metFrag_files", compi = compi, selection = selection)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
