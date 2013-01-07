\name{diss.PRED}
\alias{diss.PRED}
\alias{multidiss.PRED}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dissimilarity Measure Based on Nonparametric Forecast%%  ~~function to do ... ~~
}
\description{
Computes the dissimilarity between two time series as the L1 distance between the kernel estimators of their forecast densities at a pre-specified horizon.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
diss.PRED(x, y, h = 5, B=500, logarithms=c(FALSE, FALSE), differences=c(0,0), plot=FALSE)
multidiss.PRED( series, h=5, B=500, logarithms=NULL, differences=NULL, plot=FALSE)
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
 \item{logarithms}{
 	Boolean vector. Specifies whether to transform each series by taking logarithms or not.
}
\item{differences}{
	Numeric vector. Specifies the amount of differences to apply to each series.
	}
\item{plot}{
	If \code{TRUE}, plot the resulting forecast densities.
}
\item{series}{
	numeric matrix. Each row specifies one time series.
}
}
\details{
The dissimilarity between the time series \code{x} and \code{y}  is given by \deqn{ d(x,y) = \int{ | f_{x,h}(u) - f_{y,h}(u) | du} } where \eqn{f_{x,h}}{ d(x,y) = INT( | f_{x,h}(u) - f_{y,h}(u) )du} and \eqn{f_{y,h}} are kernel density estimators of the forecast densities h-steps ahead of \code{x} and \code{y}, respectively. The horizon of interest h is pre-specified by the user.
The kernel density estimators are based on B bootstrap replicates obtained by using a resampling procedure that mimics the generating processes, which are assumed to follow an arbitrary autoregressive structure (parametric or non-parametric). The procedure is completely detailed in Vilar et al. (2010). This function has high computational cost due to the bootstrapping procedure. 
\code{multidiss.PRED} computes the similarity matrix for more than two series requiring different a different logarithm transform or different amount of differences. In this case the default values of \code{logarithms} are \code{FALSE} and \code{0} for \code{differences}.

}
\value{
\code{diss.PRED}returns a list with the following components.

\item{L1dist}{ The computed distance.}
\item{dens.x}{ A 2-column matrix with the density of predicion of series \code{x}. First column is the base (x) and the second column is the value (y) of the density.}
\item{dens.y}{ A 2-column matrix with the density of predicion of series \code{y}. First column is the base (x) and the second column is the value (y) of the density.}

\code{multidiss.PRED} returns a list with the following components.
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

\examples{
x <- (rnorm(100))
x <- x + abs(min(x)) + 1 #shift to produce values greater than 0, for a correct logarithm transform
y <- (rnorm(100))
z <- sin(seq(0, pi, length.out=100))
## Compute the distance and check for coherent results
diss.PRED(x, y, 5, logarithms=c(FALSE,FALSE), differences=c(1,0))
#create a dist object for its use with clustering functions like pam or hclust
multidiss.PRED( rbind(x,y,z), h=5, B=500, logarithms=c(TRUE,FALSE, FALSE), differences=c(1,1,2) )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
