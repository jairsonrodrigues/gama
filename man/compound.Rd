\name{compound}
\docType{data}
\alias{compound}
\title{Synthetic dataset of two-dimensional points.}
\description{
  This is a synthetic dataset that contains groups of different density points, varied shapes, and necks between partitions.
}
\usage{data(compound)}
\format{A data frame containing 399 observations and two dimensions, forming six partitions:

    \enumerate{
    \item x1: synthetically generated real positive values
    \item x2: synthetically generated real positive values
    }

   Originally, the dataset had contained three dimensions. We intentionally removed the third dimension that corresponds to the label which the data point belongs. All description about the data set may be found in \emph{Graph-theoretical methods for detecting and describing gestalt clusters} article, in the references.
}

\source{The dataset was collected from \href{http://cs.joensuu.fi/sipu/datasets/}{Clustering basic benchmark} site. }

\references{
C.T. Zahn, \emph{Graph-theoretical methods for detecting and describing gestalt clusters}. IEEE Transactions on Computers, 1971. 100(1): p. 68-86.

P. Franti and S. Sieranoja, \emph{K-means properties on six clustering benchmark datasets,} vol. 48, no. 12. pp. 4743-4759, 2018.

}
\keyword{datasets}
