\name{loo1nn.cv}
\alias{loo1nn.cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Clustering Evaluation Index Based on Leave-one-out One-nearest-neighbor Evaluation
}
\description{
 Computes the leave-one-out one-nearest-neighbor cross-validation of an arbitrary distance matrix.
}
\usage{
loo1nn.cv(d, G)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{d}{
    A \code{dist} object.
}
  \item{G}{
    Integer vector with the labels of the true cluster solution. Each element of the vector specifies the cluster 'id' that the element belongs to.
}
}
\details{
    Computes the proportion of succesful clusters that the given distance matrix produces using leave-one-out one-nearest-neighbor cross-validation. Distance ties are solved by majority vote. A tie while voting produces a warning and is solved by selecting a candidate cluster at random. 
}
\value{
The computed proportion.
}

\author{
Pablo Montero Manso, José Antonio Vilar.
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
 \code{\link{cluster.evaluation}}, \code{\link[pdc]{loo1nn}},  \code{\link[class]{knn.cv}},
}

\examples{
 data(synthetic.tseries)
 
 #create the ground thruth cluster 
 G <- rep(1:6, each = 3)
 
 #obtain candidate distance matrix (dist object)
 dACF <- diss(synthetic.tseries, "ACF")
 
 #calculate the cross-validation
 loo1nn.cv(dACF, G)
 
}

