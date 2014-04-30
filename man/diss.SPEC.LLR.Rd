\name{diss.SPEC.LLR}
\alias{diss.SPEC.LLR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
General Spectral Dissimilarity Measure Using Local-Linear Estimation of the Log-Spectra%%  ~~function to do ... ~~
}
\description{
Computes a general dissimilarity measure based on the ratio of local linear spectral estimators.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
diss.SPEC.LLR(x, y, alpha=0.5, method="DLS", plot=FALSE, n=length(x))
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
 \item{alpha}{ Power for the ratio of densities in the Chernoff information measure. Between 0 and 1.}
  \item{method}{
        \code{"DLS"} for least squares estimation of the spectral density and \code{"LK"} for maximum likelihood estimation.
%%     ~~Describe \code{y} here~~
}
  \item{plot}{
		if \code{TRUE}, plot the smoothed spectral densities of the two series.	
}

\item{n}{
	The number of points to use for the linear interpolation. A value of n=0 uses numerical integration instead of linear interpolation. See details.
}

}
\details{
\deqn{ d_W = \int_{-\pi}^{\pi} W^{\prime}\Bigg( \frac{f_x(\lambda) } { f_y(\lambda) } \Bigg) d\lambda }{d_W = INT( W^\prime( f_x(\lambda) / f_y(\lambda) ) )d\lambda | \lambda=-\pi to \pi}
where:
\itemize{
    \item  \eqn{f_x} and \eqn{f_y} are nonparametric approximations of spectral densities of \code{x} and \code{y} respectively.
    \item \eqn{W^{\prime}(x) = W(x) + W(1/x)} with \eqn{W(x) = \log( \alpha x + (1- \alpha) x) - \alpha \log(x)}, so that \eqn{W(.)} is a divergence function depending on \eqn{\alpha}. 
}
This dissimilarity measure corresponds to the limiting spectral approximation of the Chernoff information measure in the time domain (see Kakizawa et al., 1998). The spectral densities are approximated by using local linear fitting by generalized least squared if \code{method=”DLS”} or by maximum likelihood if \code{method=”LK”} (in this case, higher computational cost is required).

By default, for performance reasons, the spectral densities are estimated using linear interpolation using \code{n} points. If \code{n} is 0, no linear interpolation is performed, and \code{integrate} is used to calculate the integral, using as many points as \code{integrate} sees fit.
 If the dissimilarity will be calculated for more than two series, calling SPEC.LLR from the \code{diss} wrapper function is preferred, since it saves some computations.

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
Vilar, J.A. and Pértega, S. (2004) Discriminant and cluster analysis for gaussian stationary processes:
  local linear fitting approach. \emph{J. Nonparametr. Stat.}, \bold{16(3-4)} 443--462.\cr
  
Kakizawa, Y.,  Shumway, R. H. and Taniguchi M. (1998) Discrimination and clustering for multivariate time series. \emph{J. Amer. Statist. Assoc.}, \bold{93(441)}, 328--
340.

Montero, P and Vilar, J.A. (2014) \emph{TSclust: An R Package for Time Series Clustering.}  Journal of Statistical Software, 62(1), 1-43. \url{http://www.jstatsoft.org/v62/i01/.}
}
\author{
Pablo Montero Manso, José Antonio Vilar.
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{diss.SPEC.GLK}}, \code{\link{diss.SPEC.ISD}}, \code{\link{integrate}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Create three sample time series
x <- cumsum(rnorm(50))
y <- cumsum(rnorm(50))
z <- sin(seq(0, pi, length.out=50))
## Compute the distance and check for coherent results
diss.SPEC.LLR(x, y, plot=TRUE)
diss.SPEC.LLR(x, z, n=0) #try integrate instead of interpolation
diss.SPEC.LLR(y, z, method="LK", n=0) #maximum likelihood with integration
#create a dist object for its use with clustering functions like pam or hclust
diss(rbind(x,y,z), METHOD="SPEC.LLR", method="DLS", alpha=0.5, n=50)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
