\name{cpu.als}
\docType{data}
\alias{cpu.als}
\title{CPU usage metrics for distributed ALS algorithm}
\description{
  Consumption metrics gathered during an execution of the Distributed Machine Learning algorithm Alternating Least Squares (ALS) in an eigth-node cluster, by using the Spark framework.
}
\usage{cpu.als}
\format{A data frame containing 308 observations and four dimensions:

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

\source{The data was measured and collected by the author by using Intel HiBench benchmark framework in a eigth-node Spark cluster hosted in Google Cloud DataProc engine. Each node had 16-core CPU annd 106 GB RAM. The algorithm Alternating Least Squares had consumed 3.7 min of runtime to execute over a sinteticaly generated dataset totalizing 1.68 Gigabytes. }

\references{

D. Goldberg, D. Nichols, B. M. Oki, and D. Terry, \emph{Using collaborative filtering to weave an information tapestry,} Commun. ACM, vol. 35, no. 12, pp. 61???70, 1992.

Y. Koren, R. Bell, and C. Volinsky, \emph{Matrix factorization techniques for recommender systems,} Computer (Long. Beach. Calif)., vol. 42, no. 8, 2009.

Y. Hu, Y. Koren, and C. Volinsky, \emph{Collaborative filtering for implicit feedback datasets,} in Data Mining, 2008. ICDM'08. Eighth IEEE International Conference on, 2008, pp. 263???272.

S. Huang, J. Huang, J. Dai, T. Xie, and B. Huang, \emph{The HiBench benchmark suite: Characterization of the MapReduce-based data analysis,} in 2010 IEEE 26th International Conference on Data Engineering Workshops (ICDEW 2010), 2010, pp. 41???51.

}
\keyword{datasets}
