\name{cpu.pca}
\docType{data}
\alias{cpu.pca}
\title{CPU usage metrics for distributed PCA algorithm}
\description{
  Consumption metrics gathered during an execution of the Distributed Machine Learning algorithm Principal Component Analysis (PCA) in an eigth-node cluster, by using the Spark framework.
}
\usage{cpu.pca}
\format{A data frame containing 938 observations and four dimensions:

    \enumerate{
    %\item node: the identifier of each seven** nodes (from 0 to 6)
    \item user: CPU usage by the algorithm
    \item system: CPU usage spent by Operating System (O.S.)
    \item iowait: waiting time for Input/Output (I/O) operations
    \item softirq: CPU time spent by software interrupt requests
    }

    The values comprise the domain from 0 to 100, for all dimensions. The dataset contains zero-values, however there is no missing or null values.

    ** A spark cluster of N nodes has 1 (one) master node and N-1 slave nodes.
}

\source{The data was measured and collected by the author by using Intel HiBench benchmark framework in a eigth-node Spark cluster hosted in Google Cloud DataProc engine. Each node had 16-core CPU annd 106 GB RAM. The algorithm Principal Component Analytsis had consumed 11.2 min of runtime to execute over a sinteticaly generated dataset totalizing 1.68 Gigabytes. }

\references{
    J.Shlens,\emph{A Tutorial on Principal Component Analysis,} Epidemiology, vol. 2, no. c, pp. 223???228, 2005.

Jolliffe, I.T.: \emph{Principal Component Analysis}, Second Edition. Encycl. Stat. Behav. Sci. 30, 487 (2002).

S. Huang, J. Huang, J. Dai, T. Xie, and B. Huang, \emph{The HiBench benchmark suite: Characterization of the MapReduce-based data analy- sis,} in 2010 IEEE 26th International Conference on Data Engineering Workshops (ICDEW 2010), 2010, pp. 41???51.

}
\keyword{datasets}
