\name{diss.INT.PER}
\alias{diss.INT.PER}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Integrated Periodogram Based Dissimilarity%%  ~~function to do ... ~~
}
\description{
Computes the dissimilarity between two time series in terms of the distance between their integrated periodograms.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
diss.INT.PER(x, y, normalize=TRUE)
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
 \item{normalize}{
 If \code{TRUE} the normalized version is computed.
 }
}
\details{

The distance is computed as:
\deqn{ d(x,y) = \int_{-\pi}^{\pi} | F_x(\lambda) - F_y(\lambda) | \, d\lambda, }{ d(x,y) = INT( F_x(\lambda) - F_y(\lambda))d\lambda | \lambda = -\pi to \pi}
 where \eqn{ F_x(\lambda_j) = C_x^{-1} \sum_{i=1}^{j} I_x(\lambda_i)}{ F_x(\lambda_j)=(sum F_x(\lambda_i)|i=1 to j)/C_x} and \eqn{F_y(\lambda_j) = C_y^{-1} \sum_{i=1}^{j} I_y(\lambda_i)}{ F_y(\lambda_j)=(sum F_y(\lambda_i)|i=1 to j)/C_xy}, with \eqn{C_x = \sum_i I_x(\lambda_i)}  and \eqn{C_y = \sum_i I_y(\lambda_i)} in the normalized version. \eqn{C_x = 1}  and \eqn{C_y = 1} in the non-normalized version. \eqn{I_x(\lambda_k)} and \eqn{I_y(\lambda_k)} denote the periodograms of \code{x} and \code{y}, respectively.
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
Casado de Lucas, D. (2010) Classification techniques for time series and functional data.

Montero, P and Vilar, J.A. (2014) \emph{TSclust: An R Package for Time Series Clustering.}  Journal of Statistical Software, 62(1), 1-43. \url{http://www.jstatsoft.org/v62/i01/.}

}
\author{
Pablo Montero Manso, José Antonio Vilar. 
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{diss.PER}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Create three sample time series
x <- cumsum(rnorm(100))
y <- cumsum(rnorm(100))
z <- sin(seq(0, pi, length.out=100))
## Compute the distance and check for coherent results
diss.INT.PER(x, y, normalize=TRUE)
diss.INT.PER(x, y, normalize=TRUE)
diss.INT.PER(x, y, normalize=TRUE)
\donttest{
diss( rbind(x,y,z), "INT.PER", normalize=FALSE )
}
}

