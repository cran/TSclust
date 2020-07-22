\name{TSclust}
\alias{TSclust-package}
\alias{TSclust}
\docType{package}
\title{
Package for Time Series Clustering.
}

\description{
This package contains several measures of dissimilarity between time series, some examples of time series datasets, specific clustering algorithms, and dimension reduction algorithms.
dissimilarities begin with diss.*, and a wrapper function \code{diss} is available. Cluster evaluation methods include \code{cluster.evaluation} and \code{loo1nn.cv}. A clustering algorithm based on pairwise p-values is implemented in \code{pvalues.clust}. The package should be used along with other existing clustering packages and function such as \code{hclust}, packages \code{cluster}, ...
}

\references{
Montero, P and Vilar, J.A. (2014) \emph{TSclust: An R Package for Time Series Clustering.}  Journal of Statistical Software, 62(1), 1-43. \url{http://www.jstatsoft.org/v62/i01/.}
}

\author{
Pablo Montero Manso, José Antonio Vilar. 
}

\seealso{
 \code{\link{diss}}
}

\examples{
#the available dissimilarities can be found in the diss help, page (?diss)
#and their individual pages from there.

### The most common use case begins with a set of time series we want to cluster.
### This package includes several example datasets.
### 
data(interest.rates)
###transformation of the interest rates
trans.inter.rates <- log(interest.rates[2:215,]) - log(interest.rates[1:214,])

##use the dist function of the proxy package to easily create the dist object
#applying ACF with geometric decaying to each pair of time series
tsdist <- diss( t(trans.inter.rates) , "ACF", p=0.05)

names(tsdist) <- colnames(interest.rates)

#perform hierachical clustering to the dist object
hc <- hclust(tsdist)

#show the results
plot(hc)

mahdist <- diss( t(trans.inter.rates) , "AR.MAH", p=0.05)$p_value

pvalues.clust(mahdist, 0.05)




}
