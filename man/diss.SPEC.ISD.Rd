\name{diss.SPEC.ISD}
\alias{diss.SPEC.ISD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dissimilarity Based on the Integrated Squared Difference between the Log-Spectra%%  ~~function to do ... ~~
}
\description{
Computes the dissimilarity between two time series in terms of the integrated squared difference between non-parametric estimators of their log-spectra.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
diss.SPEC.ISD(x, y, plot=FALSE)
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
	\item{plot}{
		If \code{TRUE}, plot the smoothed spectral densities of the two series.
}
}
\details{

\deqn{ d(x,y) = \int ( \hat{m}_x(\lambda) - \hat{m}_y(\lambda))^2 \, d\lambda, }{ d(x,y) = INT( (m_x(\lambda) - m_y(\lambda))^2 )d\lambda }
where \eqn{ \hat{m}_x(\lambda) }{m_x(\lambda)} and \eqn{ \hat{m}_y(\lambda) }{m_y(\lambda)} are local linear smoothers of the log-periodograms, obtained using the maximum local likelihood criterion.
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
Pértega, S. and Vilar, J.A. (2010) Comparing several parametric and nonparametric approaches to time series clustering: A simulation study. \emph{J. Classification}, \bold{27(3)}, 333--362.
}
\author{
Pablo Montero Manso, José Antonio Vilar.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{diss.SPEC.GLK}}, \code{\link{diss.SPEC.LLR}} 
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Create two sample time series
x <- cumsum(rnorm(50))
y <- cumsum(rnorm(50))
z <- sin(seq(0, pi, length.out=100))
## Compute the distance and check for coherent results
diss.SPEC.ISD(x, y, plot=TRUE)
#create a dist object for its use with clustering functions like pam or hclust
\dontrun{
proxy::dist( rbind(x,y,z), diss.SPEC.ISD )
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
