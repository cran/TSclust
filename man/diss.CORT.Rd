\name{diss.CORT}
\alias{diss.CORT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dissimilarity Index Combining Temporal Correlation and Raw Values Behaviors%%  ~~function to do ... ~~
}
\description{
Computes an adaptive dissimilarity index between two time series that covers both dissimilarity
on raw values and dissimilarity on temporal correlation behaviors.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
diss.CORT(x, y, k = 2, deltamethod="Euclid")
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
  Parameter controlling the weight of the dissimilarity between dynamic behaviors (See Details).
%%     ~~Describe \code{k} here~~
}
  \item{deltamethod}{
  	Defines the method for the raw data discrepancy. Either \code{"Euclid"}, \code{"Frechet"} or \code{"DTW"}.
  }
}
\details{
The dissimilarity between  time series \code{x}  and \code{y} is given by:
\deqn{ d(x,y) = \Phi[CORT(x,y)] \delta(x,y) }
where:

    CORT(x,y) measures the proximity between the dynamic behaviors of x and y by means of the first order temporal correlation coefficient defined by:
\deqn{ CORT(x,y) = \frac{ \sum_{t=1} (x_{t+1} - x_t) ( y_{t+1} - y_t) }{ \sqrt{ \sum_{t=1} (x_{t+1} - x_t)^2}  \sqrt{ \sum_{t=1} (y_{t+1} - y_t)^2  } } }{ CORT(x,y) = \sum ( ( x_(t+1) - x_t ) ( y_(t+1) - y_t ) ) /  (\sqrt( \sum  (x_(t+1) - x_t)^2) \sqrt( \sum  (y_(t+1) - y_t)^2))}
    \eqn{\Phi[u]} is an adaptive tuning function taking the form:
 \deqn{  \frac{2}{1+e^{ku}} }{2/(1+exp(ku))} with \eqn{k \geq 0} so that both \eqn{\Phi} and \code{k} modulate the weight that CORT(x,y) has on d(x,y).
 
 \eqn{\delta(x,y)} denotes a dissimilarity measure between the raw values of series \code{x} and \code{y}, such as the Euclidean distance, the Frechet distance or the Dynamic Time Warping distance. Note that \eqn{ d(x,y) = \delta(x,y)} if \code{k=0}.

More details of the procedure can be seen in Chouakria-Douzal and Nagabhushan (2007).
 
\code{deltamethod} (\eqn{\delta}) can be either Euclidean (\code{deltamethod = "Euclid"}), Frechet (\code{ deltamethod = "Frechet"}) or Dynamic Time Warping (\code{deltamethod ="DTW"}) distances. When calling from \code{dis.CORT}, DTW uses Manhattan as local distance.


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
Chouakria-Douzal, A. and Nagabhushan P. N. (2007) Adaptive dissimilarity index for measuring time series proximity. \emph{Adv. Data Anal. Classif.}, \bold{1(1)}, 5--21.
}
\author{
Pablo Montero Manso, José Antonio Vilar. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{diss.COR}}, \code{\link{diss.DTWARP}}, \code{\link{diss.FRECHET}}, \code{\link[longitudinalData]{distFrechet}}, \code{\link[dtw]{dtw}}.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Create three sample time series
x <- cumsum(rnorm(100))
y <- cumsum(rnorm(100))
z <- sin(seq(0, pi, length.out=100))
## Compute the distance and check for coherent results
diss.CORT(x, y, 2)
diss.CORT(x, z, 2)
diss.CORT(y, z, 2)
#create a dist object for its use with clustering functions like pam or hclust
\dontrun{
diss( rbind(x,y,z), "CORT", k=3, deltamethod="DTW")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
