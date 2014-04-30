\name{diss}
\alias{diss}

\title{
TSclust Dissimilarity Computation
}
\description{
Computes the dissimilarity matrix of the given numeric matrix, list, data.frame or \code{mts} object using the selected TSclust dissimilarity method.
}
\usage{
diss(SERIES, METHOD, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SERIES}{
       Numeric matrix, \code{list}, \code{data.frame} or \code{mts} object. Numeric matrices are interpreted row-wise (one series per row) meanwhile \code{data.frame} and \code{mts} objects are interpredted column-wise.
}
  \item{METHOD}{
      the dissimilarity measure to be used. This must be one of "ACF", "AR.LPC.CEPS", "AR.MAH", "AR.PIC", "CDM", "CID", "COR", "CORT", "DTWARP", "DWT", "EUCL", "FRECHET", INT.PER",
                "NCD", "PACF", "PDC",  PER", "PRED", "MINDIST.SAX", "SPEC.LLR", "SPEC.GLK" or "SPEC.ISD". Any unambiguous substring can be given.
                 See details for individual usage.
}
  \item{\dots}{
      Additional arguments for the selected method.
}

}
\details{
\code{SERIES} argument can be a numeric matrix, with one row per series, a \code{list} object with one numeric vector per element, a \code{data.frame} or a \code{mts} object.
Some methods can have additional arguments. See the individual help page for each dissimilarity method, detailed below.
Methods that have arguments that require one value per time series in \code{series} must provide so using a vector, a matrix (in the case of a multivalued argument) or a list when appropiate. In the case of a matrix, the values are conveyed row-wise. See the AR.LPC.CEPS example below.
\itemize{
\item "ACF" Autocorrelation-based method. See \code{\link{diss.ACF}}.
\item "AR.LPC.CEPS" Linear Predictive Coding ARIMA method. This method has two value-per-series arguments, the ARIMA order, and the seasonality.See \code{\link{diss.AR.LPC.CEPS}}.
\item "AR.MAH" Model-based ARMA method. See \code{\link{diss.AR.MAH}}.
\item "AR.PIC" Model-based ARMA method. This method has a value-per-series argument, the ARIMA order. See \code{\link{diss.AR.PIC}}.
\item "CDM" Compression-based dissimilarity method. See \code{\link{diss.CDM}}.
\item "CID" Complexity-Invariant distance. See \code{\link{diss.CID}}.
\item "COR" Correlation-based method. See \code{\link{diss.COR}}.
\item "CORT" Temporal Correlation and Raw values method. See \code{\link{diss.CORT}}.
\item "DTWARP" Dynamic Time Warping method. See \code{\link{diss.DTWARP}}.
\item "DWT" Discrete wavelet transform method. See \code{\link{diss.DWT}}.
\item "EUCL" Euclidean distance. See \code{\link{diss.EUCL}}. For many more convetional distances, see \code{link[stats]{dist}}, though you may need to transpose the dataset.
\item "FRECHET" Frechet distance. See \code{\link{diss.FRECHET}}.
\item "INT.PER" Integrate Periodogram-based method. See \code{\link{diss.INT.PER}}.
\item "NCD" Normalized Compression Distance. See \code{\link{diss.NCD}}.
\item "PACF" Partial Autocorrelation-based method. See \code{\link{diss.PACF}}.
\item "PDC" Permutation distribution divergence. Uses the \code{\link{pdc}} package. See \code{\link[pdc]{pdcDist}} for 
additional arguments and details. Note that series given by numeric matrices are interpreted row-wise and not column-wise, opposite as in \code{\link[pdc]{pdcDist}}.
\item "PER" Periodogram-based method. See \code{\link{diss.PER}}.
\item "PRED" Prediction Density-based method. This method has two value-per-series agument, the logarithm and difference transform. See \code{\link{diss.PRED}}.
\item "MINDIST.SAX" Distance that lower bounds the Euclidean, based on the Symbolic Aggregate approXimation measure. See \code{\link{diss.MINDIST.SAX}}.
\item "SPEC.LLR" Spectral Density by Local-Linear Estimation method. See \code{\link{diss.SPEC.LLR}}.
\item "SPEC.GLK" Log-Spectra Generalized Likelihood Ratio test method. See \code{\link{diss.SPEC.GLK}}.
\item "SPEC.ISD" Intregated Squared Differences between Log-Spectras method. See \code{\link{diss.SPEC.ISD}}.
}
}
\value{
\item{dist}{A \code{dist} object with the pairwise dissimilarities between \code{series}.}
Some methods produce additional output, see their respective documentation pages for more information.
}

\author{
Pablo Montero Manso, José Antonio Vilar. 
%%  ~~who you are~~
}
\seealso{
 \code{\link[pdc]{pdc}}, \code{\link[stats]{dist}}
}
\examples{
data(electricity)
diss(electricity, METHOD="INT.PER", normalize=FALSE)

## Example of multivalued, one per series argument
## The AR.LPC.CEPS dissimilarity allows the specification of the ARIMA model for each series
## Create three sample time series and a mts object
x <- arima.sim(model=list(ar=c(0.4,-0.1)), n =100, n.start=100)
y <- arima.sim(model=list(ar=c(0.9)), n =100, n.start=100)
z <- arima.sim(model=list(ar=c(0.5, 0.2)), n =100, n.start=100)
seriests <- rbind(x,y,z)

## If we want to provide the ARIMA order for each series
## and use it with AR.LPC.CEPS, we create a matrix with the row-wise orders
orderx <- c(2,0,0) 
ordery <- c(1,0,0)
orderz <- c(2,0,0)
orders = rbind(orderx, ordery, orderz)

diss( seriests, METHOD="AR.LPC.CEPS", k=30, order= orders )

##other examples
diss( seriests, METHOD="MINDIST.SAX", w=10, alpha=4 )
diss( seriests, METHOD="PDC" )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
