\name{paired.tseries}
\alias{paired.tseries}
\docType{data}
\title{
Pairs of Time Series from Different Domains
}
\description{
  Dataset formed by pairs of time series from different domains. Series were selected from the UCR Time Series Archive.
}
\details{
	Each pair of series in the dataset (Series 1 and 2, Series 3 and 4, etc.) comes from the same domain, so this pairing could constitute a possible ground truth solution.
}
\note{
 \code{abbreviate} can be used on the \code{colnames}.
}
\usage{data(paired.tseries)}
\format{
  A \code{mts} object with 36 series of length 1000.
}


\source{http://www.cs.ucr.edu/~eamonn/SIGKDD2004/All_datasets/}
\references{
Keogh, E., Lonardi, S., & Ratanamahatana, C. A. (2004). Towards parameter-free data mining. Proceedings of the tenth ACM SIGKDD international conference on Knowledge discovery and data mining (pp. 206-215).

Montero, P and Vilar, J.A. (2014) \emph{TSclust: An R Package for Time Series Clustering.}  Journal of Statistical Software, 62(1), 1-43. \url{http://www.jstatsoft.org/v62/i01/.}

}
\examples{
data(paired.tseries)
#Create the true solution, the pairs
true_cluster <- rep(1:18, each=2)
#test a dissimilarity metric and a cluster algorithm
intperdist <- diss( paired.tseries, "INT.PER") #create the distance matrix
#use hierarchical clustering and divide the tree in 18 clusters
intperclust <- cutree( hclust(intperdist), k=18 )
#use a cluster simmilarity index to rate the solution
cluster.evaluation( true_cluster, intperclust)

#### other evaluation criterion used in this dataset  consist in counting the correct pairs
#### formed during agglomerative hierarchical cluster (see references)
true_pairs = (-matrix(1:36, ncol=2, byrow=TRUE))
hcintper <- hclust(intperdist, "complete")
#count within the hierarchical cluster the pairs
sum( match(data.frame(t(true_pairs)), data.frame(t(hcintper$merge)), nomatch=0) > 0 ) / 18

}
\keyword{datasets}
