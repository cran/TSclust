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
diss.SPEC.LLR(x, y, alpha=0.5, method="LS", plot=FALSE)
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
        \code{"LS"} for least squares estimation of the spectral density and \code{"LK"} for maximum likelihood estimation.
%%     ~~Describe \code{y} here~~
}
  \item{plot}{
		if \code{TRUE}, plot the smoothed spectral densities of the two series.	
}

}
\details{
\deqn{ d_W = \int_{-\pi}^{\pi} W^{\prime}\Bigg( \frac{f_x(\lambda) } { f_y(\lambda) } \Bigg) d\lambda }{d_W = INT( W^\prime( f_x(\lambda) / f_y(\lambda) ) )d\lambda | \lambda=-\pi to \pi}
where:
\itemize{
    \item  \eqn{f_x} and \eqn{f_y} are nonparametric approximations of spectral densities of \code{x} and \code{y} respectively.
    \item \eqn{W^{\prime}(x) = W(x) + W(1/x)} with \eqn{W(x) = \log( \alpha x + (1- \alpha) x) - \alpha \log(x)}, so that \eqn{W(.)} is a divergence function depending on \eqn{\alpha}. 
}
This dissimilarity measure corresponds to the limiting spectral approximation of the Chernoff information measure in the time domain (see Kakizawa et al., 1998). The spectral densities are approximated by using local linear fitting by generalized least squared if \code{method=”LS”} or by maximum likelihood if \code{method=”LK”} (in this case, higher computational cost is required).

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
}
\author{
Pablo Montero Manso, José Antonio Vilar.
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{diss.SPEC.GLK}, \link{diss.SPEC.ISD}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Create three sample time series
x <- cumsum(rnorm(50))
y <- cumsum(rnorm(50))
z <- sin(seq(0, pi, length.out=50))
## Compute the distance and check for coherent results
diss.SPEC.LLR(x, y, plot=TRUE)
diss.SPEC.LLR(x, z)
diss.SPEC.LLR(y, z)
#create a dist object for its use with clustering functions like pam or hclust
\dontrun{
proxy::dist( rbind(x,y,z), method=distance.W, meth="LK")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
