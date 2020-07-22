\name{cluster.evaluation}
\alias{cluster.evaluation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Clustering Evaluation Index Based on Known Ground Truth%%  ~~function to do ... ~~
}
\description{
Computes the similarity between the true cluster solution and the one obtained with a method under evaluation.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cluster.evaluation(G, S)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{
    Integer vector with the labels of the true cluster solution. Each element of the vector specifies the cluster 'id' that the element belongs to.
%%     ~~Describe \code{x} here~~
}
  \item{S}{
    Integer vector with the labels of the cluster solution to be evaluated. Each element of the vector specifies the cluster 'id' that the element belongs to.
%%     ~~Describe \code{y} here~~
}

}
\details{
The measure of clustering evaluation is defined as \deqn{ Sim(G,C) = 1/k \sum_{i=1}^k \max_{1\leq j\leq k} Sim(G_i,C_j),  } where \deqn{Sim(G_i, C_j) = \frac{ 2 | G_i \cap C_j|}{ |G_i| + |C_j|}}

with |.| denoting the cardinality of the elements in the set. This measure has been used for comparing different clusterings, e.g. in Kalpakis et al. (2001) and Pértega and Vilar (2010). 

}
\value{
The computed index.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Larsen, B. and Aone, C. (1999) Fast and effective text mining using linear-time document clustering. \emph{Proc. KDD' 99}.16--22. \cr

Kalpakis, K., Gada D. and Puttagunta, V. (2001) Distance measures for effective clustering of arima time-series. \emph{Proceedings 2001 IEEE International Conference on Data Mining}, 273--280. \cr

Pértega S. and Vilar, J.A (2010) Comparing several parametric and nonparametric approaches to time series clustering: A simulation study. \emph{J. Classification}, \bold{27(3)}, 333-362.

Montero, P and Vilar, J.A. (2014) \emph{TSclust: An R Package for Time Series Clustering.}  Journal of Statistical Software, 62(1), 1-43. \url{http://www.jstatsoft.org/v62/i01/.}

}
\author{
Pablo Montero Manso, José Antonio Vilar.
%%  ~~who you are~~
}
\note{
This index is not simmetric.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
 \code{\link[fpc]{cluster.stats}},  \code{\link[clValid]{clValid}}, \code{\link[clv]{std.ext}}
}

\examples{

 #create a true cluster 
 #(first 4 elements belong to cluster '1', next 4 to cluster '2' and the last 4 to cluster '3'.
 true_cluster <- c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3)
 #the cluster to be tested
 new_cluster <- c( 2, 1, 2, 3, 3, 2, 2, 1, 3, 3, 3, 3)
 
 #get the index
 cluster.evaluation(true_cluster, new_cluster)
 
 #it can be seen that the index is not simmetric
 cluster.evaluation(new_cluster, true_cluster)
}
