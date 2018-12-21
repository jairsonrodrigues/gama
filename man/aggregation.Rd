\name{aggregation}
\docType{data}
\alias{aggregation}
\title{Synthetic dataset of two-dimensional points.}
\description{
  This is a synthetic dataset that contains features that are known to create difficulties for the selected algorithms such as, narrow bridges between clusters, uneven-sized clusters, etc. See references, for details.
}
\usage{data(aggregation)}
\format{A data frame containing 788 observations and two dimensions, forming seven partitions:

    \enumerate{
    \item x1: synthetically generated real positive values
    \item x2: synthetically generated real positive values
    }

   Originally, the dataset had contained three dimensions. We intentionally removed the third dimension that corresponds to the label which the data point belongs. All description about the data set may be found in \emph{Clustering Aggregation} article, in the references.
}

\source{The dataset was collected from \href{http://cs.joensuu.fi/sipu/datasets/}{Clustering basic benchmark} site. }

\references{
A. Gionis, H. Mannila, and P. Tsaparas, \emph{Clustering aggregation}. ACM Transactions on Knowledge Discovery from Data (TKDD), 2007. 1(1): p. 1-30.

P. Franti and S. Sieranoja, \emph{K-means properties on six clustering benchmark datasets,} vol. 48, no. 12. pp. 4743-4759, 2018.

}
\keyword{datasets}
