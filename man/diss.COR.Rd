\name{diss.COR}
\alias{diss.COR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Correlation-based Dissimilarity%%  ~~function to do ... ~~
}
\description{
Computes dissimilarities based on the estimated Pearson's correlation of two given time series.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
diss.COR(x, y, beta = NULL)
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
  \item{beta}{
  	If not NULL, specifies the regulation of the convergence in the second method.
  }

}
\details{
 Two different measures of dissimilarity between two time series based on the estimated Pearson's correlation can be computed.
 If \code{beta} is not specified, the value \eqn{ d_1 = \sqrt{ 2 ( 1 - \rho) } } is computed, where \eqn{(\rho)} denotes the Pearson's correlation between series \code{x} and \code{y}.
 If \code{beta} is specified, the function \eqn{ d_2 = \sqrt{ (\frac{ 1 - \rho}{ 1 + \rho})^\beta } }{ d_2 = \sqrt( ((1 - \rho) / (1 + \rho))^\beta )} is used, where \eqn{\beta} is \code{beta} .
 

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
	
Golay, X., Kollias, S., Stoll, G., Meier, D., Valavanis, A., and Boesiger, P. (2005) A new correlation-based fuzzy logic clustering algorithm for FMRI. \emph{Magnetic Resonance in Medicine}, \bold{40.2}, 249--260.

Montero, P and Vilar, J.A. (2014) \emph{TSclust: An R Package for Time Series Clustering.}  Journal of Statistical Software, 62(1), 1-43. \doi{doi:10.18637/jss.v062.i01}

}
\author{
Pablo Montero Manso, José Antonio Vilar.  
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{diss.PACF}}, \code{\link{diss.ACF}}, \code{\link{diss}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Create three sample time series
x <- cumsum(rnorm(100))
y <- cumsum(rnorm(100))
z <- sin(seq(0, pi, length.out=100))
## Compute the distance and check for coherent results
diss.COR(x, y)
diss.COR(x, z)
#create a dist object for its use with clustering functions like pam or hclust
\donttest{
diss( rbind(x,y,z), "COR")
}
}

