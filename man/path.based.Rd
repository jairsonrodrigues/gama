\name{path.based}
\docType{data}
\alias{path.based}
\title{Circular cluster.}
\description{
  The data set consists of a circular cluster with an opening near the bottom and two Gaussian distributed clusters inside. All description about the data set may be found in \emph{Robust path-based spectral clustering} article, and the Clustering basic benchmark site, in the references.
}
\usage{data(path.based)}
\format{A data frame containing 240 observations and two dimensions, forming two partitions:

    \enumerate{
    \item x1: real positive values
    \item x2: real positive values
    }

   Originally, the dataset had contained three dimensions. We intentionally removed the third dimension that corresponds to the label which the data point belongs.
}

\source{The dataset was collected from \href{http://cs.joensuu.fi/sipu/datasets/}{Clustering basic benchmark} site. }

\references{
H. Chang and D.Y. Yeung, \emph{Robust path-based spectral clustering}. Pattern Recognition, 2008. 41(1): p. 191-203.

P. Franti and S. Sieranoja, \emph{K-means properties on six clustering benchmark datasets,} vol. 48, no. 12. pp. 4743-4759, 2018.

}
\keyword{datasets}
