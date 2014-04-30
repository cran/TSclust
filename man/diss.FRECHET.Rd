\name{diss.FRECHET}
\alias{diss.FRECHET}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Frechet Distance%%  ~~function to do ... ~~
}
\description{
Computes the Frechet distance between time series.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
diss.FRECHET(x, y, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  Numeric vector containing the first of the two time series.
%%     ~~Describe \code{x} here~~
}
  \item{y}{
  Numeric vector containing the second of the two time series.
%%     ~~Describe \code{y} here~~
}
\item{\dots}{
  Additional parameters for the function. See \code{link[longitudinalData]{distFrechet}}.
}

}

\details{
 This is basically a wrapper for \code{distFrechet} of the \code{\link{longitudinalData}} package, intended for an easier discovery of the functionalities used in TSclust.
}

\value{
The computed distance.

}



\seealso{
 \code{\link{diss}}, \code{link[longitudinalData]{distFrechet}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Create three sample time series
x <- cumsum(rnorm(100))
y <- cumsum(rnorm(100))
diss.FRECHET(x, y)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
