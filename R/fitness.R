# Maximizing solutions through Average Silhouette Width (ASW) criterion;
# from -1 (the worst) to 1 (the best).
fitness.asw <- function(individual, penalty.function = NULL, ...) {

  k <- gama.env$k
  data <- gama.env$dataset
  d2 <- gama.env$d2

  dims <- length(individual)/k
  fitness.value <- NA

  m.individual <- matrix(individual, nrow = k, ncol = dims)
  which.dists <- apply(Rfast::dista(data, m.individual, "euclidean", square = TRUE), 1, which.min)

  # to avoid a convergence for configurations different from user-specified k
  if (length(unique(which.dists)) < k) {
    # maximum penalty
    fitness.value = -1

  } else {

    # calculate the average silhouette width
    asw <- cluster::silhouette(which.dists, d2)

    # try summarize the silhouette, returns (0) if error
    fitness.value <- tryCatch(summary(asw)$avg.width, error = function (e) { return (0)})

    # optionally, some inequalities may be resolved by applying a penality
    if (!is.null(penalty.function)) {
      penalty <- penalty.function(m.individual)
      fitness.value <- fitness.value - penalty
    }
  }
  return (fitness.value)
}

# Maximizing solutions through Calinski Harabasz (CH) criterion.
# The higher the value, the "better" is the solution.
fitness.ch <- function(individual, penalty.function = NULL, ...) {

  k <- gama.env$k
  data <- gama.env$dataset
  d2 <- gama.env$d2

  dims <- length(individual)/k
  fitness.value <- NA

  m.individual <- matrix(individual, nrow = k, ncol = dims)
  which.dists <- apply(Rfast::dista(data, m.individual, "euclidean", square = TRUE), 1, which.min)

  # to avoid a convergence for configurations different from user-specified k
  # the Calinski_Harabasz index does not have a specific boundary, like ASW or C-Index
  # thus, the -10^5 assignment as maximum penalty.
  if (length(unique(which.dists)) < k) {
    fitness.value = 0

  } else {

    # calculate the Calinski_Harabasz index
    ch <- clusterCrit::intCriteria(as.matrix(data), which.dists, "Calinski_Harabasz")
    fitness.value <- ch$calinski_harabasz

    # optionally, some inequalities may be resolved by applying a penality
    if (!is.null(penalty.function)) {
      penalty <- penalty.function(m.individual)
      fitness.value <- fitness.value - penalty
    }
  }
  return (fitness.value)

}

# Evaluate individuals according C-Indez (CI) criterion.
# Minimize solutions through C-Index (CI) criterion; from 1 (the worst) to 0 (the best).
#
# C-Index varies between 0 and 1. The lower is the value the better is the cluster partition.
# REF: Hubert, L.J., Levin, J.R. A general statistical framework for assessing categorical clustering
# in free recall. Psychol. Bull., 1976, 83, 1072-1080
#
# https://stats.stackexchange.com/questions/343878/computation-of-c-index-for-cluster-validation
fitness.ci <- function(individual, penalty.function = NULL, ...) {

  k <- gama.env$k
  data <- gama.env$dataset
  d2 <- gama.env$d2

  dims <- length(individual)/k
  fitness.value <- NA

  m.individual <- matrix(individual, nrow = k, ncol = dims)
  which.dists <- apply(Rfast::dista(data, m.individual, "euclidean", square = TRUE), 1, which.min)

  # to avoid a convergence for configurations different from user-specified k
  if (length(unique(which.dists)) < k) {
    # maximum penalty
    fitness.value = 1

  } else {

    # calculate the C index
    ch <- clusterCrit::intCriteria(as.matrix(data), which.dists, "C_Index")
    fitness.value <- ch$c_index

    # optionally, some inequalities may be resolved by applying a penality
    if (!is.null(penalty.function)) {
      penalty <- penalty.function(m.individual)
      fitness.value <- fitness.value + penalty
    }
  }
  return (1 - fitness.value)

}

# Maximize solutions through Dunn index (CH) criterion.
# Dunn index are usually used to identify the "compact and well separated clusters".
# The main drawback of Dunn's index is computational since calculating becomes
# computationally very expansive as number of clusters and number of data points
# increase. In the case of overlapped clusters the values of Dunn Index is not really
# reliable because of re-partitioning the results with the hard partition method. (Balasko, et al, 2005)
fitness.di <- function(individual, penalty.function = NULL, ...) {

  k <- gama.env$k
  data <- gama.env$dataset
  d2 <- gama.env$d2

  dims <- length(individual)/k
  fitness.value <- NA

  m.individual <- matrix(individual, nrow = k, ncol = dims)
  which.dists <- apply(Rfast::dista(data, m.individual, "euclidean", square = TRUE), 1, which.min)

  # to avoid a convergence for configurations different from user-specified k
  if (length(unique(which.dists)) < k) {
    fitness.value = 0

  } else {

    # calculate the Calinski_Harabasz index
    ch <- clusterCrit::intCriteria(as.matrix(data), which.dists, "Dunn")
    fitness.value <- ch$dunn

    # optionally, some inequalities may be resolved by applying a penality
    if (!is.null(penalty.function)) {
      penalty <- penalty.function(m.individual)
      fitness.value <- fitness.value - penalty
    }
  }
  return (fitness.value)
}
