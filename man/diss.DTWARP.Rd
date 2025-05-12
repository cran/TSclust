\name{diss.DTWARP}
\alias{diss.DTWARP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dynamic Time Warping Distance%%  ~~function to do ... ~~
}
\description{
Computes Dynamic Time Warping distance between time series.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
diss.DTWARP(x, y, ...)
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
  Additional parameters for the function. See \code{link[dtw]{dtw}}.
}

}

\details{
 This is basically a wrapper for \code{dtw} of the \code{\link[dtw]{dtw}} package, intended for an easier discovery of the functionalities used in TSclust.
}

\value{
The computed distance.

}



\seealso{
 \code{\link{diss}}, \code{link[dtw]{dtw}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Create three sample time series
x <- cumsum(rnorm(100))
y <- cumsum(rnorm(100))
diss.DTWARP(x, y)

}

