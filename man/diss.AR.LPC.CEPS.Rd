\name{diss.AR.LPC.CEPS}
\alias{diss.AR.LPC.CEPS}
\alias{multidiss.AR.LPC.CEPS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dissimilarity Based on LPC Cepstral Coefficients%%  ~~function to do ... ~~
}
\description{
Computes the dissimilarity between two time series in terms of their Linear Predicitive Coding (LPC) ARMA processes.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
diss.AR.LPC.CEPS(x, y, k = 50, order=NULL,
 seasonal=list(list(order=c(0, 0, 0), period=NA), list(order=c(0, 0, 0), period=NA)))

multidiss.AR.LPC.CEPS(series, k = 50, order=NULL, seasonal=NULL )
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
\item{seasonal}{ A list of \code{arima} seasonal elements. }
\item{order}{ Numeric matrix. Each row specifies the ARIMA models to be fitted for the series.}
\item{series}{ Numeric matrix. Each row specifies one time series.}
}
\details{
Each row of \code{order} contains the three components of the ARIMA model: the AR order, the degree of differencing and the MA order, specified as in the function \code{\link{arima}}. The first two rows correspond to series \code{x} and \code{y} respectively.
If \code{order} is NULL, or any of the orders contains \code{NA} values, the series are automatically fitted by imposing an AR model, selected by AIC.
Each element of \code{seasonal}  is a list with two components: 'order' and 'period'. See \code{seasonal} parameter of \code{\link{arima}}.
The \code{seasonal} component of a series is ignored if the \code{order} component is not specified.
\code{multidiss.AR.LPC.CEPST} computes the similarity matrix for more that two series having different orders.

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
Kalpakis, K., Gada D. and Puttagunta, V. (2001) Distance measures for effective clustering of arima time-series. \emph{Proceedings 2001 IEEE International Conference on Data Mining}, 273--280. \cr
}
\author{
Pablo Montero Manso, José Antonio Vilar. 
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{link{diss.AR.PIC}}, \code{link{diss.AR.MAH}}
}
\examples{
## Create three sample time series
x <- arima.sim(model=list(ar=c(0.4,-0.1)), n =100, n.start=100)
y <- arima.sim(model=list(ar=c(0.9)), n =100, n.start=100)
z <- arima.sim(model=list(ar=c(0.5, 0.2)), n =100, n.start=100)
## Compute the distance and check for coherent results
diss.AR.LPC.CEPS(x, y, 25) #impose an AR automatically selected for both series
#impose an ARIMA(2,0,0) for series x and an AR automatically selected for z
diss.AR.LPC.CEPS(x, z, 25, order = rbind( c(2,0,0) , c(NA,NA,NA)) ) 
diss.AR.LPC.CEPS(y, z, 25)
#create a dist object for its use with clustering functions like pam or hclust

multidiss.AR.LPC.CEPS( rbind(x,y,z), k=20, order=rbind(c(2,0,0), c(1,0,0), c(2,0,0)),
 list( list(order=c(1,0,0), period=1), list(order=c(2,0,0), period=3),
  list(order=c(1,0,0), period=1)) )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
