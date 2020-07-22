\name{diss.SPEC.GLK}
\alias{diss.SPEC.GLK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dissimilarity based on the Generalized Likelihood Ratio Test%%  ~~function to do ... ~~
}
\description{
The dissimilarity between two time series is computed by using an adaptation of the generalized likelihood ratio test to check the equality of two log-spectra.
}
\usage{
diss.SPEC.GLK(x, y, plot=FALSE)
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
	\item{plot}{
		If \code{TRUE}, the smoothed spectral densities of the two series are plotted.	
}

}
\details{
The dissimilarity between two series \code{x} and \code{y} is measured in terms of the vaue of a test statistic to check the equality of their log-spectra, \eqn{m_X(\lambda)} and \eqn{m_Y(\lambda)} respectivelty. The test statistic is constructed by using the generalized likelihood ratio test criterion (Fan and Zhang, 2004). Specifically, the test statistic takes the form:
\deqn{ d(x,y) = \sum_{k=1}^T [ Z_k - \hat{\mu}( \lambda_k) - 2 \log( 1 + e^{ \{ Z_k - \hat{\mu}(\lambda_k)\}})] - \sum_{k=1}^T[Z_k - 2 \log(1 + e^{Z_k})],}{ d(x,y) = \sum[Z_k - hat\mu(\lambda_k) - 2 \log( 1 + exp(Z_k - hat\mu(\lambda_k))] - \sum[ Z_k -2 \log(1+ exp(Z_k))] | k=1 to T,}

where \eqn{I_x(\lambda_k)} and \eqn{I_y(\lambda_k)} are the periodograms of \code{x} and \code{y}, \eqn{ Z_k = \log(I_x(\lambda_k)) - \log( I_y(\lambda_k))}, and \eqn{\hat{\mu}(\lambda_k)}{hat\mu(\lambda_k)} is the local maximum log-likelihood estimator of \eqn{\mu(\lambda_k)= m_x(\lambda_k) - m_y(\lambda_k)} computed by local linear fitting.


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
Fan, J. and Zhang, W. (2004) Generalised likelihood ratio tests for spectral density. \emph{Biometrika}, 195--209.\cr

Pértega, S. and Vilar, J.A. (2010) Comparing several parametric and nonparametric approaches to time series clustering: A simulation study. \emph{J. Classification}, \bold{27(3)}, 333--362.

Montero, P and Vilar, J.A. (2014) \emph{TSclust: An R Package for Time Series Clustering.}  Journal of Statistical Software, 62(1), 1-43. \url{http://www.jstatsoft.org/v62/i01/.}

}
\author{
Pablo Montero Manso, José Antonio Vilar.
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{diss.SPEC.ISD}}, \code{\link{diss.SPEC.LLR}} 
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Create two sample time series
x <- cumsum(rnorm(50))
y <- cumsum(rnorm(50))
z <- sin(seq(0, pi, length.out=50))
## Compute the distance and check for coherent results
diss.SPEC.GLK(x, y, plot=TRUE)
#create a dist object for its use with clustering functions like pam or hclust
\donttest{
diss( rbind(x,y,z), "SPEC.GLK" )
}
}

