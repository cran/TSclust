\name{diss.AR.PIC}
\alias{diss.AR.PIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model-based Dissimilarity Measure Proposed by Piccolo (1990)%%  ~~function to do ... ~~
}
\description{
Computes the distance between two time series as the Euclidean distance between the truncated AR operators approximating their ARMA structures.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
diss.AR.PIC(x, y, order.x=NULL, order.y=NULL, permissive=TRUE)

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

\item{order.x}{ Specifies the ARIMA model to be fitted for the series x. When using \code{diss} wrapper, use \code{order} argument instead. See details.}
\item{order.y}{ Specifies the ARIMA model to be fitted for the series y. When using \code{diss} wrapper, use \code{order} argument instead. See details.}
\item{permissive}{ Specifies whether to force an AR order of 1 if no order is found. Ignored if neither order.x or order.y are NULL}
}
\details{
If \code{order.x} or \code{order.y} are \code{NULL}, their respective series will be fitted automatically using a AR model. If \code{permissive} is \code{TRUE} and no AR order is found automatically, an AR order of 1 will be imposed, if this case fails, then no order can be found and the function produces an error.
\code{order.x} and \code{order.y} contain the three components of the ARIMA model: the AR order, the degree of differencing and the MA order, specified as in the function \code{\link{arima}}.

If using \code{diss} function with "AR.PIC" \code{method}, the argument \code{order} must be used instead of \code{order.x} and \code{order.y}. \code{orders} is a matrix with one row per ARIMA, specified as in \code{\link{arima}}. If \code{order} is \code{NULL}, automatic fitting imposing a AR model is performed.

}

\value{
The computed distance.
}

\references{
Piccolo, D. (1990) A distance measure for classifying arima models. \emph{J. Time Series
Anal.}, \bold{11(2)}, 153--164.

Montero, P and Vilar, J.A. (2014) \emph{TSclust: An R Package for Time Series Clustering.}  Journal of Statistical Software, 62(1), 1-43. \doi{doi:10.18637/jss.v062.i01}

}

\author{
Pablo Montero Manso, José Antonio Vilar.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{diss.AR.MAH}}, \code{\link{diss.AR.LPC.CEPS}}, \code{\link{diss}}, \code{\link{arima}}, \code{\link{ar}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Create three sample time series
x <- arima.sim(model=list(ar=c(0.4,-0.1)), n =100, n.start=100)
y <- arima.sim(model=list(ar=c(0.9)), n =100, n.start=100)
z <- arima.sim(model=list(ar=c(0.5, 0.2)), n =100, n.start=100)
## Compute the distance and check for coherent results
#ARIMA(2,0,0) for x and ARIMA(1,0,0) for y
diss.AR.PIC( x, y, order.x = c(2,0,0), order.y = c(1,0,0) )
diss.AR.PIC( x, z, order.x = c(2,0,0), order.y = c(2,0,0) )
# AR for y (automatically selected) and ARIMA(2,0,0) for z
diss.AR.PIC( y, z, order.x=NULL, order.y=c(2,0,0) ) 
#create a dist object for its use with clustering functions like pam or hclust
diss( rbind(x,y,z), METHOD="AR.PIC", order=rbind(c(2,0,0), c(1,0,0), c(2,0,0)) )

}
