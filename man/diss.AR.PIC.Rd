\name{diss.AR.PIC}
\alias{diss.AR.PIC}
\alias{multidiss.AR.PIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model-based Dissimilarity Measure Proposed by Piccolo (1990)%%  ~~function to do ... ~~
}
\description{
Computes the distance between two time series as the Euclidean distance between the truncated AR operators approximating their ARMA structures.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
diss.AR.PIC(x, y, order=NULL)
multidiss.AR.PIC( series, order=NULL)
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
\item{order}{ Numeric matrix. Each row specifies the ARIMA models to be fitted for the series.}
\item{series}{ Numeric matrix. Each row specifies one time series.}
}
\details{
Each row of \code{order} contains the three components of the ARIMA model: the AR order, the degree of differencing and the MA order, specified as in the function \code{\link{arima}}. The first two rows correspond to series \code{x} and \code{y} respectively.
If \code{order} is NULL, or any of the orders contains \code{NA} values, the corresponding series are fitted automatically imposing an AR model, selected automatically by AIC.
\code{multidiss.AR.PIC} computes the similarity matrix for more that two series having different orders.
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
Piccolo, D. (1990) A distance measure for classifying arima models. \emph{J. Time Series
Anal.}, \bold{11(2)}, 153--164.
}
\author{
Pablo Montero Manso, José Antonio Vilar.
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{diss.AR.MAH}} \code{\link{diss.AR.LPC.CEPS}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Create three sample time series
x <- arima.sim(model=list(ar=c(0.4,-0.1)), n =100, n.start=100)
y <- arima.sim(model=list(ar=c(0.9)), n =100, n.start=100)
z <- arima.sim(model=list(ar=c(0.5, 0.2)), n =100, n.start=100)
## Compute the distance and check for coherent results
diss.AR.PIC(x, y, rbind( c(2,0,0), c(1,0,0)) ) #ARIMA(2,0,0) for x and ARIMA(1,0,0) for y
diss.AR.PIC(x, z, rbind( c(2,0,0), c(2,0,0)) )
# AR for y (automatically selected) and ARIMA(2,0,0) for z
diss.AR.PIC(y, z, rbind( c(NA,NA,NA), c(2,0,0)) ) 
#create a dist object for its use with clustering functions like pam or hclust
multidiss.AR.PIC( rbind(x,y,z), order=rbind(c(2,0,0), c(1,0,0), c(2,0,0)) )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
