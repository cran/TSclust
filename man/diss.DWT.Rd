
\name{diss.DWT}

\alias{diss.DWT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dissimilarity for Time Series Based on Wavelet Feature Extraction%%  ~~function to do ... ~~
}
\description{
 Performs an unsupervised feature extration using orthogonal wavelets on the series and returns the Euclidean distance between the wavelet approximations in an appropriate scale.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
diss.DWT(series)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{series}{
  Numeric matrix with row order time series
%%     ~~Describe \code{x} here~~
}

}
\details{
This method differs from other dissimilarities in that pairwise dissimilaries depend on the whole dataset that is given to \code{diss.DWT}, hence, there is no pairwise version of the function defined, only accepts whole datasets.
The set of original series is replaced by their wavelet approximation coefficients in an appropriate scale, and the dissimilarity between two series is computed as the Euclidean distance between these coefficients.  The appropriate scale is automatically determined by using an algorithm addressed to obtain an efficient reduction of the dimensionality but preserving as much information from the original data as possible. The algorithm is introduced by Zhang, Ho, Zhang, and Lin (2006).
 
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns an object of type \code{dist} with the pairwise distances.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\seealso{
\code{\link[wmtsa]{wavDWT}}
}
\references{
Zhang, H., Ho, T. B., Zhang, Y., and Lin, M. (2006) Unsupervised feature extraction for time series clustering using orthogonal wavelet transform. \emph{INFORMATICA-LJUBLJANA-}, \bold{30(3)}, 305.
}
\author{
Pablo Montero Manso, José Antonio Vilar.
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## Create three sample time series
x <- cumsum(rnorm(100))
y <- cumsum(rnorm(100))
z <- sin(seq(0, pi, length.out=100))

#compute the distance
diss.DWT(rbind(x, y, z))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
