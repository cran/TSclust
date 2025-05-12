\name{diss.PER}
\alias{diss.PER}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Periodogram Based Dissimilarity%%  ~~function to do ... ~~
}
\description{
Computes the distance between two time series based on their periodograms.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
diss.PER(x, y, logarithm=FALSE, normalize=FALSE)
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
  \item{logarithm}{
  Boolean. If \code{TRUE} logarithm of the periodogram coefficients will be taken.
  }
  \item{normalize}{
  Boolean. If \code{TRUE}, the periodograms will be normalized by the variance of their respective series.
  }
}
\details{
	Computes the Euclidean distance between the periodogram coefficients of the series \code{x} and \code{y}. Additional transformations can be performed on the coefficients depending on the values of \code{logarithm} and \code{normalize}.
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
Caiado, J., Crato, N. and Peña, D. (2006) A periodogram-based metric for time series classification. \emph{Comput. Statist. Data Anal.}, \bold{50(10)}, 2668--2684.

Montero, P and Vilar, J.A. (2014) \emph{TSclust: An R Package for Time Series Clustering.}  Journal of Statistical Software, 62(1), 1-43. \doi{doi:10.18637/jss.v062.i01}

}
\author{
Pablo Montero Manso, José Antonio Vilar.
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{link{diss.INT.PER}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Create three sample time series
x <- cumsum(rnorm(100))
y <- cumsum(rnorm(100))
z <- sin(seq(0, pi, length.out=100))
## Compute the distance and check for coherent results
diss.PER(x, y)
diss.PER(x, z)
diss.PER(y, z)
diss.PER(x, y, TRUE, TRUE)
diss.PER(x, z, TRUE, TRUE)
diss.PER(y, z, TRUE, TRUE)
#create a dist object for its use with clustering functions like pam or hclust
diss( rbind(x,y,z), "PER", logarithm=TRUE, normalize=TRUE)
}
