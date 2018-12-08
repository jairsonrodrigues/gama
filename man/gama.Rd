\name{gama}
\alias{gama}
\title{Segment a dataset by using genetic search.}
\description{
  This function realizes a genetic search in order to segment a given dataset. The genetic algorithm evaluates individuals of a population to maximize (or minimize) a one of four user-specified criterion: Average Silhouette Width, Calinski-Harabasz index, C-index or Dunn index. The algorithm receives the optimal number of partitions (best k value) from the user or automatically suggests it.

}
\usage{
gama(data, k = NA, scale = FALSE, crossover.rate = 0.9,
                   mutation.rate = 0.01, elitism = 0.05, pop.size = 25,
                   generations = 100, seed.p = 42, fitness.criterion =
                   "ASW", penalty.function = NULL, plot.internals = TRUE, \dots)
}
\arguments{
  \item{data}{the dataset to segment into k partitions. Must be numerical, because GAMA only works for integer ou real values.}
  \item{k}{the best value for the optimal number of partitions. The user may ommit this argument, in this way, the algorithm will use internal validation criterias for automatically estimate the best k value for the data.}
  \item{scale}{if the dataset will be scaled and normalized before segmentation. Default value is FALSE.}
  \item{crossover.rate}{the probability of crossover between pairs of chromosomes. Default value is 0.9.}
  \item{mutation.rate}{the probability of mutation in a parent chromosome. Usually mutation occurs with a small probability. Default value is 0.01.}
  \item{elitism}{the number of best individuals to survive. Default value is the top 5\% individuals of population.}
  \item{pop.size}{the number of individuals in the population. Default value is 25. Observation: this argument have a great impact on the performance of the search, due to increased number of matrix calculations when the population grows.}
  \item{generations}{the number of generations to execute the search. Default value is 100.}
  \item{seed.p}{an integer value containing the random number generator.}
  \item{fitness.criterion}{the key point of the genetic search. The algorithm will search for the ideal centroids that maximizes one of the pre-specified criteria:
    \enumerate{
      \item "ASW": Average Silhouette Width;
      \item "CH": Calinski Harabasz index;
      \item "DI": Dunn index;
      \item "CI": C-index.
    }

    The default value is "ASW" and will be assumed if the user does not supply value or put an invalid entry.
  }
  \item{penalty.function}{...}
  \item{plot.internals}{...}
  \item{\dots}{other arguments that user may pass to the function.}

}
\references{
Scrucca, L. (2013) GA: \emph{A Package for Genetic Algorithms in R}. Journal of Statistical Software, 53(4), 1-37. http://www.jstatsoft.org/v53/i04/

P. J. Rousseeuw, \emph{Silhouettes: A graphical aid to the interpretation and validation of cluster analysis,} J. Comput. Appl. Math., vol. 20, no. C, pp. 53???65, 1987.

T. Calinski and J. Harabasz. A dendrite method for cluster analysis. Com- munications in Statistics, 3, no. 1:1???27, 1974.

J. Dunn. Well separated clusters and optimal fuzzy partitions. Journal of Cybernetics, 4:95???104, 1974.

Hubert, L.J., Levin, J.R. A general statistical framework for assessing categorical clustering in free recall. Psychol. Bull., 1976, 83, 1072-1080

}
\examples{
\dontrun{
# loads data about CPU execution metrics of a distributed
# version of Alternating Least Squares (ALS) algorithm
data(cpu.als)

# a specific function to calculate penalties for CPU execution metrics
# whose does not allow the sum of loads above 100\%
my.penalty <- function(m.individual,...) {

  penalty <- 0

  # counts how many centroids results in overflow (inequality > 100)
  sums <- apply(m.individual, 1, sum)
  overflow <- which(sums > 100)
  num_constraints = length(overflow)

  # if there are overflows, subtract each dimension
  # by the maximum proportion of the excess x the number of overflows
  if (num_constraints > 0) {
     penalty <- num_constraints * max(abs(sums[overflow] -100)/sums[overflow])
  }

  return (penalty)
}

# call the gama clustering to maximize ASW criterion (default)
# and delegates to GAMA choose the best k value
gama.obj <- gama(data = cpu.als, penalty.function = my.penalty)

print(gama.obj)
}
}

\keyword{file}
