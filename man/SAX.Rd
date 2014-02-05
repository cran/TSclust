\name{diss.MINDIST.SAX}
\alias{diss.MINDIST.SAX}
\alias{diss.SAX}
\alias{PAA}
\alias{convert.to.SAX.symbol}
\alias{MINDIST.SAX}
\alias{SAX.plot}
\concept{PAA}

\title{
Symbolic Aggregate Aproximation related functions
}
\description{
\code{diss.MINDIST.SAX} computes a dissimilarity that lower bounds the Euclidean on the discretized, dimensionality reduced series. Function \code{PAA} produces the dimension reduction. Function \code{convert.to.SAX.symbol} produces the discretization.
}
\usage{
diss.MINDIST.SAX(x, y, w, alpha=4, plot=FALSE)
PAA(x, w)
convert.to.SAX.symbol(x, alpha)
MINDIST.SAX(x, y, alpha, n)
SAX.plot(series, w, alpha, col.ser=rainbow(ncol(series)))
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
  \item{w}{
  	The amount of equal sized frames that the series will be reduced to.
%%     ~~Describe \code{k} here~~
}
  \item{alpha}{
  	The size of the alphabet, the amount of symbols used to represents the values of the series.
%%     ~~Describe \code{k} here~~
}

\item{plot}{
	If \code{TRUE}, plot a graphic of the reduced series, with their corresponding symbols.
}

\item{n}{
	The original size of the series.
}

\item{series}{
    A \code{ts} or \code{mts} object with the series to plot.
}
\item{col.ser}{
    Colors for the series. One per series.
}

}
\details{

SAX is a symbolic representation of continuous time series.


\code{w} must be an integer but it does not need to divide the length of the series. If \code{w} divides the length of the series, the \code{diss.MINDIST.SAX} plot uses this to show the size of the frames.

\code{PAA} performs the Piecewise Aggregate Approximation of the series, reducing it to \code{w} elements, called frames. Each frame is composed by \eqn{n/w} observations of the original series, averaged. Observations are weighted when \code{w} does not divide \code{n}.

\code{convert.to.SAX.symbol} performs SAX discretization: Discretizes the series \code{x} to an alphabet of size \code{alpha}, \code{x} should be z-normalized in this case. The \eqn{N(0,1)} distribution is divided in \code{alpha} equal probability parts, if an observation falls into the \eqn{i}th part (starting from minus infinity), it is assigned the \eqn{i} symbol.

\code{MINDIST.SAX} calculates the MINDIST dissimilarity between symbolic representations.

\code{diss.MINDIST.SAX} combines the previous procedures to compute a dissimilarity between series. The series are z-normalized at first. Then the dimensionality is reduced uusin \code{PAA} to produce series of length \code{w}. The series are discretized to an alphabet of size \code{alpha} using \code{convert.to.SAX.symbol}. Finally the dissimilarity value is produced using \code{MINDIST.SAX}.

\code{SAX.plot} produces a plot of the SAX representation of the given \code{series}.


}
\value{
The computed dissimilarity.
}
\references{

Lin, J., Keogh, E., Lonardi, S. & Chiu, B. (2003) A Symbolic Representation of Time Series, with Implications for Streaming Algorithms. In Proceedings of the 8th ACM SIGMOD Workshop on Research Issues in Data Mining and Knowledge Discovery.

Keogh, E., Chakrabarti, K., Pazzani, M., & Mehrotra, S. (2001). Dimensionality reduction for fast similarity search in large time series databases. \emph{Knowledge and information Systems}, \bold{3(3)}, 263-286.

}

\author{
Pablo Montero Manso, José Antonio Vilar. 
}

\seealso{
 \code{\link{diss}}
}

\examples{
set.seed(12349)
n = 100
x <- rnorm(n)  #generate sample series, white noise and a wiener process
y <- cumsum(rnorm(n))
w <- 20 #amount of equal-sized frames to divide the series, parameters for PAA
alpha <- 4 #size of the alphabet, parameter for SAX

#normalize
x <- (x - mean(x)) /sd(x)
y <- (y - mean(y)) /sd(y)

paax <- PAA(x, w) #generate PAA reductions
paay <- PAA(y, w)

plot(x, type="l", main="PAA reduction of series x") #plot an example of PAA reduction
p <- rep(paax,each=length(x)/length(paax)) #just for plotting the PAA
lines(p, col="red")

#repeat the example with y
plot(y, type="l", main="PAA reduction of series y") 
py <- rep(paay,each=length(y)/length(paay))
lines(py, col="blue")

#convert to SAX representation
SAXx <- convert.to.SAX.symbol( paax, alpha)
SAXy <- convert.to.SAX.symbol( paay, alpha)

#CALC THE SAX DISTANCE
MINDIST.SAX(SAXx, SAXy, alpha, n)

#this whole process can be computed using diss.MINDIST.SAX
diss.MINDIST.SAX(x, y, w, alpha, plot=TRUE)

z <- rnorm(n)^2

diss(rbind(x,y,z), "MINDIST.SAX", w, alpha)

SAX.plot( as.ts(cbind(x,y,z)), w=w, alpha=alpha)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
