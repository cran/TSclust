\name{diss.AR.LPC.CEPS}
\alias{diss.AR.LPC.CEPS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dissimilarity Based on LPC Cepstral Coefficients%%  ~~function to do ... ~~
}
\description{
Computes the dissimilarity between two time series in terms of their Linear Predicitive Coding (LPC) ARIMA processes.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
diss.AR.LPC.CEPS(x, y, k = 50, order.x=NULL, order.y=NULL,
 seasonal.x=list(order=c(0, 0, 0), period=NA),
 seasonal.y=list(order=c(0, 0, 0), period=NA),
 permissive=TRUE)

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
  \item{k}{
  Number of cepstral coefficients to be considered.
%%     ~~Describe \code{k} here~~
}
\item{order.x}{ Numeric matrix. Specifies the ARIMA models to be fitted for the series x. When using \code{diss} wrapper, use \code{order} argument instead. See details.}
\item{order.y}{ Numeric matrix. Specifies the ARIMA ARIMA models to be fitted for the series y. When using \code{diss} wrapper, use \code{order} argument instead. See details.}
\item{seasonal.x}{ A list of \code{arima} seasonal elements for series x. When using \code{diss} wrapper, use \code{seasonal} argument instead. See details.}
\item{seasonal.y}{ A list of \code{arima} seasonal elements for series x. When using \code{diss} wrapper, use \code{seasonal} argument instead. See details.}
\item{permissive}{ Specifies whether to force an AR order of 1 if no order is found. Ignored if neither order.x or order.y are NULL}
}
\details{
If \code{order.x} or \code{order.y} are \code{NULL}, their respective series will be fitted automatically using a AR model.
\code{order.x} and \code{order.y} contain the three components of the ARIMA model: the AR order, the degree of differencing and the MA order, specified as in the function \code{\link{arima}}.

\code{seasonal.x} and \code{seasonal.y} are lists with two components: 'order' and 'period'. See \code{seasonal} parameter of \code{\link{arima}}, except that specification using a numeric \code{vector} of length 3 is not allowed.

If using \code{diss} function with "AR.LPC.CEPS" \code{method}, the argument \code{order} must be used instead of \code{order.x} and \code{order.y}. \code{order} is a matrix with one row per series, specified as in \code{\link{arima}}. If \code{order} is \code{NULL}, automatic fitting imposing a AR model is performed. The argument \code{seasonal} is used instead of \code{seasonal.x} and \code{seasonal.y}. \code{seasonal} is a list of elements, one per series in the same order that the series are input. Each element of \code{seasonal} must have the same format as the one in \code{\link{arima}}.


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
Kalpakis, K., Gada D. and Puttagunta, V. (2001) Distance measures for effective clustering of arima time-series. \emph{Proceedings 2001 IEEE International Conference on Data Mining}, 273--280. \cr

Montero, P and Vilar, J.A. (2014) \emph{TSclust: An R Package for Time Series Clustering.}  Journal of Statistical Software, 62(1), 1-43. \doi{doi:10.18637/jss.v062.i01}
}
\author{
Pablo Montero Manso, José Antonio Vilar. 
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{diss.AR.PIC}}, \code{\link{diss.AR.MAH}}, \code{\link{diss}}
}
\examples{
## Create three sample time series
x <- arima.sim(model=list(ar=c(0.4,-0.1)), n =100, n.start=100)
y <- arima.sim(model=list(ar=c(0.9)), n =100, n.start=100)
z <- arima.sim(model=list(ar=c(0.5, 0.2)), n =100, n.start=100)
## Compute the distance and check for coherent results
diss.AR.LPC.CEPS(x, y, 25) #impose an AR automatically selected for both series
#impose an ARIMA(2,0,0) for series x and an AR automatically selected for z
diss.AR.LPC.CEPS(x, z, 25, order.x = c(2,0,0), order.y = NULL ) 
diss.AR.LPC.CEPS(y, z, 25)
#create a dist object for its use with clustering functions like pam or hclust

diss( rbind(x,y,z), METHOD="AR.LPC.CEPS", k=20, order=rbind(c(2,0,0), c(1,0,0), c(2,0,0)),
 seasonal=list( list(order=c(1,0,0), period=1), list(order=c(2,0,0), period=3),
  list(order=c(1,0,0), period=1)) )

}
