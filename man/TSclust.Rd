\name{TSclust}
\alias{TSclust-package}
\docType{package}
\title{
Package for Time Series Clustering.
}

\description{
This package contains several measures of dissimilarity between time series, some examples of time series datasets, specific clustering algorithms, and dimension reduction algorithms.
}


\examples{
### The most common use case begins with a set of time series we want to cluster.
### This package includes several example datasets.
### 
data(interest.rates)
###transformation of the interest rates
trans.inter.rates <- log(interest.rates[2:215,]) - log(interest.rates[1:214,])

##use the dist function of the proxy package to easily create the dist object
#applying ACF with geometric decaying to each pair of time series
if (require(proxy)) {
  tsdist <- proxy::dist( t(trans.inter.rates) , diss.ACF, p=0.05)

  names(tsdist) <- colnames(interest.rates)

  #perform hierachical clustering to the dist object
  hc <- hclust(tsdist)

  #show the results
  plot(hc)
}
}
