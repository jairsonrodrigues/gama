\name{gama}
\alias{gama}
\title{Segments a dataset by using genetic search.}
\description{
  This function realizes a genetic search in order to segment a given dataset. The genetic algorithm evaluates individuals of a population to maximize (or minimize) one of four user-specified criteria: Average Silhouette Width, Calinski-Harabasz index, C-index or Dunn index. The algorithm receives the optimal number of partitions (best k value) from the user or automatically suggests it.

}
\usage{
gama(dataset = NULL, k = "broad", scale = FALSE, crossover.rate = 0.9,
                   mutation.rate = 0.01, elitism = 0.05, pop.size = 25,
                   generations = 100, seed.p = 42, fitness.criterion =
                   "ASW", penalty.function = NULL, plot.internals = TRUE, \dots)
}
\arguments{
  \item{dataset}{the data to segment into k partitions. Must be numerical, because GAMA only works for integer ou real values.}
  \item{k}{the best value for the optimal number of partitions. May be:
  \enumerate{
    \item an integer positive value: when known, the user may specify directly the best value for k;

    \item a string value: 'minimal', for estimation based on 'elbow' in Within-cluster Sum of Squares Error graph or 'broad', the alternative to estimate the best k by majority voting between 24 distinct internal validation criteria for automatically estimate the best k value for the data.
  }

  The user may omit this argument, in this way, the default value 'broad' will be used.
  }



  \item{scale}{if the dataset will be scaled and normalized before segmentation. The default value is FALSE.}
  \item{crossover.rate}{the probability of crossover between pairs of chromosomes. The default value is 0.9.}
  \item{mutation.rate}{the probability of mutation in a parent chromosome. Usually mutation occurs with a small probability. The default value is 0.01.}
  \item{elitism}{the number of best individuals to survive. The default value is the top 5\% individuals of population.}
  \item{pop.size}{the number of individuals in the population. The default value is 25. Observation: this argument have a great impact on the performance of the search, due to increased number of matrix calculations when the population grows.}
  \item{generations}{the number of generations to execute the search. The default value is 100.}
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
  \item{penalty.function}{an optional user-specified function to be applied over fitness results, to penalize undesirable or incongruent individuals.}
  \item{plot.internals}{if TRUE, the evolution of the fitness values and a Silhouette graph will be display after the genetic search. The default value is TRUE.}
  \item{\dots}{other arguments that user may pass to the function.}


}
\details{

A call to \code{gama} function should at least contain the argument \code{data}. In this case, the clustering will be done without penalities over individuals, aiming to maximize the \emph{silhouette} criterion, and the function will suggest the best value for \emph{k}, the optimal number of partitions. For the others arguments, the default values applied will be as defined in the call: no normalization or scaling, rates for crossover, mutation, and elitism equals 90\%, 1\%, and 5\%, respectivelly, population size of 25 individuals, 100 generations, and plot.internals = TRUE.

A lot of other combinations of calls is possible, especially choosing the criteria to guide the maximization of the search ("ASW", "CH", "CI", "DI"), the k value, if known, and adjustments in the genetic parameters, like number of generations, population size, and mutation/crossover/elitism rates.

The function \code{gama} returns an S4 object of class "gama". This object contains the following slots:

\enumerate{
	\item {original.data: } {the original dataset used for clustering.}
	\item {centers: }{the k (partitions) x d (dimensions) matrix representation of the centroids who best segment the data into k partitions.}
	\item {cluster: }{a vector of integers (from 1:k) indicating the cluster to which each point is allocated.}
	\item {silhouette: }{the value for Average Silhouette Width index.}
	\item {calinski_harabasz: }{the value for Calinski\_Harabasz index.}
	\item {c_index: }{the value for C-index.}
	\item {dunn_index: }{the value for Dunn index.}
	\item {runtime: }{the total time spent by the clustering algorithm.}
	\item {call: }{the string representation of the user call to the function \code{gama} and its parameters.}
}


}
\references{
Scrucca, L. (2013) GA: \emph{A Package for Genetic Algorithms in R}. Journal of Statistical Software, 53(4), 1-37. http://www.jstatsoft.org/v53/i04/

P. J. Rousseeuw, \emph{Silhouettes: A graphical aid to the interpretation and validation of cluster analysis,} J. Comput. Appl. Math., vol. 20, no. C, pp. 53???65, 1987.

T. Calinski and J. Harabasz. A dendrite method for cluster analysis. Com- munications in Statistics, 3, no. 1:1-27, 1974.

J. Dunn. Well separated clusters and optimal fuzzy partitions. Journal of Cybernetics, 4:95-104, 1974.

Hubert, L.J., Levin, J.R. A general statistical framework for assessing categorical clustering in free recall. Psychol. Bull., 1976, 83, 1072-1080

}
\examples{

# loading flame dataset
data(flame)

# segmentation of the flame dataset in k = 2 partitions
# the 'plot.internals' says to the program do not plot the graphs about
# genetic evolution search and silhouette index
gama.obj <- gama(flame, k = 4, plot.internals = FALSE, generations = 30)
# ** use at least 100 generations for simple datasets and 500 for complex datasets

# it draws the partitions to which each element belongs
gama.plot.partitions(gama.obj)

\dontrun{

# loads data about CPU execution metrics of a distributed
# version of Alternating Least Squares (ALS) algorithm
data(cpu.als)

# a user-defined function to calculate penalties for CPU execution metrics
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

# call the gama clustering to maximize Dunn index criterion
# by using 500 generations and delegates to GAMA to choose the best k value
gama.obj <- gama(data = cpu.als, fitness.criterion = "DI",
                generations = 500, penalty.function = my.penalty)

print(gama.obj)

}
}

\keyword{file}
