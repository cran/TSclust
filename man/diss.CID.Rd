\name{diss.CID}
\alias{diss.CID}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Complexity-Invariant Distance Measure For Time Series
}
\description{
Computes the distance based on the Euclidean distance corrected by the complexity estimation of the series.
}
\usage{
diss.CID(x, y)
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
\details{
	This distance is defined \deqn{CID(x,y) = ED(x,y) \times CF(x,y)} where \eqn{CF(x,y)} is a complexity correction factor defined as: \deqn{ max(CE(x), CE(y)) / min(CE(x), CE(y)) } and \eqn{CE(x)} is a compexity estimate of a time series \eqn{x}. \code{diss.CID} therefore increases the distance between series with different complexities. If the series have the same complexity estimate, the distance defenerates Euclidean distance. The complexity is defined in \code{diss.CID} as:
		 \deqn{ CE(x) = \sqrt{ \sum_{t=1} (x_{t+1} - x_t)^2  } }{CE(x) = sqrt( sum ( diff(x)^2 ) )}
}
\value{
The computed dissimilarity.
}
\references{

Batista, G. E., Wang, X., & Keogh, E. J. (2011). A Complexity-Invariant Distance Measure for Time Series. In SDM (Vol. 31, p. 32).

}

\author{
Pablo Montero Manso, José Antonio Vilar. 
}

\seealso{
 \code{\link{diss}}, \code{\link{diss.CORT}}
}

\examples{
n = 100
x <- rnorm(n)  #generate sample series, white noise and a wiener process
y <- cumsum(rnorm(n))

diss.CID(x, y)

z <- rnorm(n)
w <- cumsum(rnorm(n))
series = rbind(x, y, z, w)
diss(series, "CID")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
