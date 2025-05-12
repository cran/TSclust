\name{diss.ACF}
\alias{diss.ACF}
\alias{diss.PACF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Autocorrelation-based Dissimilarity%%  ~~function to do ... ~~
}
\description{
Computes the dissimilarity between two time series as the distance between their estimated simple (ACF) or partial (PACF) autocorrelation coefficients.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
diss.ACF(x, y, p = NULL, omega=NULL, lag.max=50)
diss.PACF(x, y, p = NULL, omega=NULL, lag.max=50)
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
  \item{p}{
  	If not NULL, sets the weight for the geometric decaying of the autocorrelation coefficients. Ranging  from \code{0} to \code{1}.
  }
  \item{lag.max}{
  	Maximum number of simple or partial autocorrelation coefficients to be considered.
  }
  \item{omega}{
  	If not NULL, completely specifies the weighting matrix for the autocorrelation coefficients. \code{p} is ignored if \code{omega} is used.
  }
}
\details{
Performs the weighted Euclidean distance between the simple autocorrelation ( \code{dist.ACF}) or partial autocorrelation ( \code{dist.PACF} ) coefficients.
If neither \code{p} nor \code{omega} are specified, uniform weighting is used. If \code{p} is specified, geometric wights decaying with the lag in the form \eqn{ p(1-p)^i} are applied. If \code{omega} (\eqn{\Omega}) is specified, \deqn{ d(x,y) = {\{ ( \hat{\rho}_{x} - \hat{\rho}_{y} )^t \bm{\Omega} (\hat{\rho}_{x} - \hat{\rho}_{y} ) \}}^\frac{1}{2} }{ d(x,y) = ((\rho_x - \rho_y)^t \Omega (\rho_x - \rho_y) )^0.5} with \eqn{\hat{\rho}_{x}}{\rho_x} and \eqn{\hat{\rho}_{y}}{\rho_y} the respective (partial) autocorrelation coefficient vectors.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The computed distance.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Galeano, P. and Peña, D. (2000). Multivariate analysis in vector time series. \emph{Resenhas}, \bold{4 (4)}, 383--403.

Montero, P and Vilar, J.A. (2014) \emph{TSclust: An R Package for Time Series Clustering.}  Journal of Statistical Software, 62(1), 1-43. \doi{doi:10.18637/jss.v062.i01}

}
\author{
Pablo Montero Manso, José Antonio Vilar.
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{diss.COR}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Create three sample time series
x <- cumsum(rnorm(100))
y <- cumsum(rnorm(100))
z <- sin(seq(0, pi, length.out=100))
## Compute the distance and check for coherent results
diss.PACF(x, y)
diss.ACF(x, z)
diss.PACF(y, z)
#create a dist object for its use with clustering functions like pam or hclust
diss( rbind(x,y,z), "ACF", p=0.05)

}
