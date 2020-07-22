\name{diss.EUCL}
\alias{diss.EUCL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Euclidean Distance%%  ~~function to do ... ~~
}
\description{
Computes Euclidean distance between time series.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
diss.EUCL(x, y)
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

}

\value{
The computed distance.

}



\seealso{
 \code{\link{diss}}, \code{\link{dist}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Create three sample time series
x <- cumsum(rnorm(100))
y <- cumsum(rnorm(100))
diss.EUCL(x, y)
}
