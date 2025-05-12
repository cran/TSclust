\name{pvalues.clust}
\alias{pvalues.clust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Clustering Algorithm Based on p-values.%%  ~~function to do ... ~~
}
\description{
Clustering algorithm based on p-values. Each group in the cluster solution is formed by
series with associated p-values greater than a pre-specified level of significance.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pvalues.clust(pvalues, significance)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pvalues}{
  A \code{dist} object containing the p-values from testing the equality of each pair of time series under study.
%%     ~~Describe \code{x} here~~
}
  \item{significance}{
  The significance level.
%%     ~~Describe \code{y} here~~
}

}
\details{
Each element (i,j) in \code{pvalues} corresponds to the p-value obtained from checking whether or not the \eqn{i}-th and \eqn{j}-th series come from the same generating
model. The clustering algorithm will only group together those series whose associated p-values are greater than the pre-specified significance level. The algorithm was originally developed for its use with the p-values obtained with in \code{diss.AR.MAH} (see Maharaj, 2000), but it can be applied to any similar test.
}
\value{
 An integer vector of length n, the number of observations, giving for each observation the number (id) of the cluster to which it belongs.

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Maharaj E.A. (2000) Clusters of time series. \emph{J. Classification}, \bold{17(2)}, 297--314.

Montero, P and Vilar, J.A. (2014) \emph{TSclust: An R Package for Time Series Clustering.}  Journal of Statistical Software, 62(1), 1-43. \doi{doi:10.18637/jss.v062.i01}

}
\author{
Pablo Montero Manso, José Antonio Vilar. 
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{diss.AR.MAH}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Create three sample time series
x <- cumsum(rnorm(100))
y <- cumsum(rnorm(100))
z <- sin(seq(0, pi, length.out=100))
##

## Compute the distance and check for coherent results
dd <- diss( rbind(x,y,z), "AR.MAH")
pvalues.clust( dd$p_value, 0.05 )

}

