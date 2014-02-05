\name{diss.NCD}
\alias{diss.NCD}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normalized Compression Distance
}
\description{
Computes the distance based on the sizes of the compressed time series.
}
\usage{
diss.NCD(x, y, type = "min")

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
 \item{type}{
 Character string, the type of compression.  May be abbreviated to a single letter, defaults to the first of the alternatives.
 }
}
\details{

The compression based  dissimilarity is calculated: \deqn{ d(x,y) = C(xy) - max(C(x),C(y))/ min(C(x),C(y)) } where \eqn{C(x)}, \eqn{C(y)} are the sizes in bytes of the compressed series \eqn{x} and \eqn{y}.
\eqn{C(xy)} is the size in bytes of the series \eqn{x} and \eqn{y} concatenated. The algorithm used for compressing the series is chosen with \code{type}.
\code{type} can be "gzip", "bzip2" or "xz", see \code{\link{memCompress}}. "min" selects the best separately for \code{x}, \code{y} and the concatenation.
Since the compression methods are character-based, a symbolic representation can be used, see details for an example using SAX as the symbolic representation.
}
\value{
The computed distance.
}
\references{
Cilibrasi, R., & Vitányi, P. M. (2005). Clustering by compression. \emph{Information Theory, IEEE Transactions on}, \bold{51(4)}, 1523-1545.

Keogh, E., Lonardi, S., & Ratanamahatana, C. A. (2004). Towards parameter-free data mining. Proceedings of the tenth ACM SIGKDD international conference on Knowledge discovery and data mining (pp. 206-215).

}

\author{
Pablo Montero Manso, José Antonio Vilar. 
}

\seealso{
 \code{\link{memCompress}}, \code{\link{diss}}
}

\examples{
n = 50
x <- rnorm(n)  #generate sample series, white noise and a wiener process
y <- cumsum(rnorm(n))

diss.NCD(x, y)

z <- rnorm(n)
w <- cumsum(rnorm(n))
series = rbind(x, y, z, w)
diss(series, "NCD", type="bzip2")

################################################################
#####symbolic representation prior to compression, using SAX####
################################################################
#normalization function, required for SAX
z.normalize = function(x) {
    (x - mean(x)) / sd(x)
}

sx <- convert.to.SAX.symbol( z.normalize(x), alpha=4 )
sy <- convert.to.SAX.symbol( z.normalize(y), alpha=4 )
sz <- convert.to.SAX.symbol( z.normalize(z), alpha=4 )
sw <- convert.to.SAX.symbol( z.normalize(w), alpha=4 )

diss(rbind(sx, sy, sz, sw), "NCD", type="bzip2")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
