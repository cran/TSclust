\name{diss.PRED}
\alias{diss.PRED}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dissimilarity Measure Based on Nonparametric Forecast%%  ~~function to do ... ~~
}
\description{
Computes the dissimilarity between two time series as the L1 distance between the kernel estimators of their forecast densities at a pre-specified horizon.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
diss.PRED(x, y, h, B=500, logarithm.x=FALSE, logarithm.y=FALSE,
differences.x=0, differences.y=0, plot=FALSE)
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
  \item{h}{
  The horizon of interest, i.e the number of steps-ahead where the prediction is evaluated. 
%%     ~~Describe \code{k} here~~
}
  \item{B}{
  The amount of bootstrap resamples. 
%%     ~~Describe \code{k} here~~
}
 \item{logarithm.x}{
 	Boolean. Specifies whether to transform series x by taking logarithms or not. When using \code{diss} wrapper, use \code{logarithms} argument instead. See details.
}
\item{logarithm.y}{
     Boolean. Specifies whether to transform series y by taking logarithms or not.  When using \code{diss} wrapper, use \code{logarithms} argument instead. See details.
}
\item{differences.x}{
	Specifies the amount of differences to apply to series x.  When using \code{diss} wrapper, use \code{differences} argument instead. See details.
	}
\item{differences.y}{
    Specifies the amount of differences to apply to series y.  When using \code{diss} wrapper, use \code{differences} argument instead. See details.
}
\item{plot}{
	If \code{TRUE}, plot the resulting forecast densities.
}

}
\details{


The dissimilarity between the time series \code{x} and \code{y}  is given by \deqn{ d(x,y) = \int{ | f_{x,h}(u) - f_{y,h}(u) | du} } where \eqn{f_{x,h}}{ d(x,y) = INT( | f_{x,h}(u) - f_{y,h}(u) )du} and \eqn{f_{y,h}} are kernel density estimators of the forecast densities h-steps ahead of \code{x} and \code{y}, respectively. The horizon of interest h is pre-specified by the user.
The kernel density estimators are based on B bootstrap replicates obtained by using a resampling procedure that mimics the generating processes, which are assumed to follow an arbitrary autoregressive structure (parametric or non-parametric). The procedure is completely detailed in Vilar et al. (2010). This function has high computational cost due to the bootstrapping procedure. 

The procedure uses a bootstrap method that requires stationary time series. In order to support a wider range of time series, the method allows some transformations on the series before proceeding with the bootstrap resampling. This transformations are inverted before calculating the densities. The transformations allowed are logarithm and differenciation.
The parameters \code{logarithm.x}, \code{logarithm.y}, \code{differences.x}, \code{differences.y} can be specified with this purpose.

If using \code{diss} function with "PRED" \code{method}, the argument \code{logarithms} must be used instead of \code{logarithm.x} and \code{logarithm.y}. \code{logarithms} is a boolean vector specifying if the logarithm transform should be taken for each one of the \code{series}. The argument \code{differences}, a numeric vector specifying the amount of differences to apply the \code{series}, is used instead of \code{differences.x} and \code{differences.y}. The plot is also different, showing all the densities in the same plot.
}
\value{
\code{diss.PRED} returns a list with the following components.

\item{L1dist}{ The computed distance.}
\item{dens.x}{ A 2-column matrix with the density of predicion of series \code{x}. First column is the base (x) and the second column is the value (y) of the density.}
\item{dens.y}{ A 2-column matrix with the density of predicion of series \code{y}. First column is the base (x) and the second column is the value (y) of the density.}

When used from the \code{diss} wrapper function, it returns a list with the following components.
\item{dist}{A \code{dist} object with the pairwise L1 distances between series.}
\item{densities}{ A list of 2-column matrices containing the densities of each series, in the same format as 'dens.x' or 'dens.y' of \code{diss.PRED}.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

}
\references{

Alonso, A.M., Berrendero, J.R., Hernandez, A. and Justel, A. (2006) Time series clustering based on forecast densities. \emph{Comput. Statist. Data Anal.}, \bold{51},762--776.\cr

Vilar, J.A.,  Alonso, A. M. and Vilar, J.M. (2010) Non-linear time series clustering based on non-parametric forecast densities. \emph{Comput. Statist. Data Anal.}, \bold{54 (11)}, 2850--2865.
}
\author{
José Antonio Vilar, Pablo Montero Manso. 
%%  ~~who you are~~
}
\seealso{
 \code{\link{diss}}
}
\examples{
x <- (rnorm(100))
x <- x + abs(min(x)) + 1 #shift to produce values greater than 0, for a correct logarithm transform
y <- (rnorm(100))
z <- sin(seq(0, pi, length.out=100))
## Compute the distance and check for coherent results
diss.PRED(x, y, h=6, logarithm.x=FALSE, logarithm.y=FALSE, differences.x=1, differences.y=0)
#create a dist object for its use with clustering functions like pam or hclust
diss( rbind(x,y,z), METHOD="PRED", h=3, B=200,
 logarithms=c(TRUE,FALSE, FALSE), differences=c(1,1,2) )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
