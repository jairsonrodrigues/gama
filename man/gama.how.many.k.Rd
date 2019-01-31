\name{gama.how.many.k}
\alias{gama.how.many.k}
\title{Estimates the optimal number of partitions.}
\description{
  This function estimates the best k value for the number of partitions the dataset should be segmented.
}
\usage{
gama.how.many.k(dataset = NULL, method = "minimal")
}
\arguments{
  \item{dataset}{the original dataset used for clustering.}
  \item{method}{the method used to estimate the number of partitions. If 'minimal' is used, the function will perform estimation based on finding the 'elbow' in the Within-cluster Sum of Squares Error graphic. It uses a second derivative approximation, in order to suggest k. If 'broad' is used, the function will proceed an estimation by majority voting of 24 indices, by using the \pkg{NbClust} package.}

}
\references{
Malika Charrad, Nadia Ghazzali, Veronique Boiteau, Azam Niknafs (2014). NbClust: An R
  Package for Determining the Relevant Number of Clusters in a Data Set. Journal of
  Statistical Software, 61(6), 1-36. URL http://www.jstatsoft.org/v61/i06/.
}
\seealso{
  \code{\link{gama}}.
}

\examples{
# loads data about CPU execution metrics of a distributed
# version of Alternating Least Squares (ALS) algorithm
library(gama)
data(cpu.als)

# call estimation by using minimal method (Elbow graphic)
k <- gama.how.many.k (cpu.als)
print(k)

# call estimation by using broad method (NbClust)
k <- gama.how.many.k (cpu.als, method = 'broad')
print(k)

}
\keyword{file}
