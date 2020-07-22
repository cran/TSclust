\name{synthetic.tseries}
\alias{synthetic.tseries}
\docType{data}
\title{
  Synthetic Time Series for Clustering Performace Comparisons.
}
\description{
	This dataset features three repetitions of several models of time series.
}
\usage{data(synthetic.tseries)}


\details{
The dataset is a \code{mts} object, formed by several repetitions of each of the following models.
\tabular{cll}{ 
M1 \tab AR \tab   \eqn{X_t = 0.6 X_{t-1} + \varepsilon_{t}}  \cr
M2 \tab Bilinear  \tab  \eqn{X_t = \left( 0.3 -0.2 \varepsilon_{t-1} \right) X_{t-1} + 1.0 +\varepsilon_{t}} \cr
M3 \tab EXPAR \tab \eqn{X_t =\left( 0.9 \exp \left( - X_{t-1}^2 \right) -0.6 \right) X_{t-1} + 1.0 + \varepsilon_{t}} \cr
M4 \tab SETAR \tab  \eqn{ X_t =\left( 0.3 X_{t-1} +1.0 \right) I \left( X_{t-1} \geq 0.2 \right) - } \cr
   \tab      \tab  \eqn{ \left( 0.3 X_{t-1} -1.0 \right) I \left( X_{t-1} < 0.2 \right) + \varepsilon_{t}} \cr
M5 \tab NLAR \tab  \eqn{ X_t = 0.7 \left| X_{t-1} \right| \left( 2 + \left| X_{t-1} \right| \right)^{-1} + \varepsilon_{t}} \cr
M6 \tab STAR \tab  \eqn{ X_t = 0.8 X_{t-1} -0.8 X_{t-1} \left( 1 + \exp \left( -10 X_{t-1} \right) \right)^{-1}  + \varepsilon_{t}} \cr
}

	Three simulations of each model are included. This dataset can be used for comparing the performance of different dissimilarity measures between time series  or clustering algorithms.
}

\references{
Montero, P and Vilar, J.A. (2014) \emph{TSclust: An R Package for Time Series Clustering.}  Journal of Statistical Software, 62(1), 1-43. \url{http://www.jstatsoft.org/v62/i01/.}
}

\examples{
data(synthetic.tseries)
#Create the true solution, for this dataset, there are three series of each model
true_cluster <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6)
#test a dissimilarity metric and a cluster algorithm
intperdist <- diss( synthetic.tseries, "INT.PER") #create the distance matrix
#use hierarchical clustering and divide the tree in 6 clusters
intperclust <- cutree( hclust(intperdist), 6 ) 
#use a cluster simmilarity index to rate the solution
cluster.evaluation( true_cluster, intperclust)

#test another dissimilarity metric and a cluster algorithm
acfdist <- diss( synthetic.tseries, "ACF", p=0.05) 
acfcluster <- pam( acfdist, 6 )$clustering #use pam clustering to form 6 clusters
cluster.evaluation( true_cluster, acfcluster)
\donttest{
#test another dissimilarity metric and a cluster algorithm
chernoffdist <- diss( synthetic.tseries, "SPEC.LLR")
chernoffclust <- pam( chernoffdist, 6 )$clustering 
cluster.evaluation( true_cluster, chernoffclust)
}

}
