# libraries are used along the code by explicit refer syntax package::function()
# library(GA)           # genetic algorithms
# library(cluster)      # average silhouette width
# library(ClusterR)     # distortion f(K)
# library(Rfast)        # matrix fast calculations
# library(ggplot2)      # graphics library
# library(clusterCrit)  # validation criteria

# source("R/fitness.R")
# source("R/bestk.R")
devtools::load_all()

if(getRversion() >= "2.15.1")  utils::globalVariables(c("dims", "clusters", "observation", "pc.1", "pc.2"))

# generates an initial population of centroids
pop.f <- function(object) {

  lower <- object@lower
  upper <- object@upper
  nvars <- length(lower)

  population <- matrix(as.double(NA), nrow = object@popSize, ncol = dims * k)

  # generate the initial population based on uniform distribution
  for(j in 1:nvars) {
    population[,j] <- runif(object@popSize, lower[j], upper[j])
  }

  # checks impossible individuals which sum(loads) > 100
  # subtracts the load value until it reach the possible range
  # for(i in 1:object@popSize) {
  #   repeat {
  #     # which lines (1:k) of an individual in population (1:popSize) has load > 100
  #     m <- matrix(population[i,], ncol = dims, nrow = k)
  #     sums <- apply(m, 1, sum)
  #     w.sums <- which(sums > 100)
  #
  #     # stop when all centroids of an individual have loads < 100
  #     if (length(w.sums) == 0){
  #       break
  #     } else {
  #       # subtracts the load by an proportional amount
  #       for (y in w.sums) {
  #         prop_over_all <- m[y,]/sums[y]
  #         population[i,] <- as.vector(m[y, ] - prop_over_all)
  #       }
  #     } # else
  #   }#repeat
  # } # for

  return(population)
}

#
#gama <- function(data, ...) { UseMethod("gama") }

gama <- function(data, k = NA, scale = FALSE, crossover.rate = 0.9,
                         mutation.rate = 0.01, elitism = 0.05, pop.size = 25,
                         generations = 100, seed.p = 42,
                         fitness.criterion = "ASW",
                         penalty.function = NULL,
                         plot.internals = TRUE, ...) {

  call <- match.call()

  # choose the correct fitness function
  fitness.function <- switch (fitness.criterion,
                             "ASW" = fitness.asw,
                             "CH" = fitness.ch,
                             "CI" = fitness.ci,
                             "DI" = fitness.di,
                             fitness.asw)

  if (scale) {
    data <- scale(data, center = TRUE, scale = TRUE)
  }

  if (is.na(k)) {
    k <- best.k(dataset = data)
  }

  .GlobalEnv$data <- data
  .GlobalEnv$k = k
  .GlobalEnv$dims = ncol(data)

  s <- "gareal_lsSelection"
  m <- "gareal_rsMutation"
  c <- "gareal_blxCrossover"

  elit.rate = floor(pop.size * elitism)

  # distance matrix
  d <- dist(data, method = "euclidean", diag = FALSE, upper = FALSE)
  d2 <- d^2

  .GlobalEnv$d2 <- d2
  rm(d); gc()

  lowers <- apply(data, 2, min)
  uppers <- apply(data, 2, max)

  lower_bound <- unlist(lapply(lowers, function (x) { rep(x, k) } ))
  upper_bound <- unlist(lapply(uppers, function (x) { rep(x, k) } ))

  # call GA functions
  # cm <- citation("GA")
  start.time <- Sys.time()

  .GlobalEnv$start.time <- start.time
  .GlobalEnv$checked.time <- start.time

  genetic <- GA::ga(type = "real-valued",
                    seed = seed.p,
                    population = pop.f,
                    selection = s,
                    mutation = m,
                    crossover = c,
                    popSize = pop.size,
                    elitism = elit.rate,
                    pmutation = mutation.rate,
                    pcrossover = crossover.rate,
                    maxiter = generations,
                    #maxFitness = 1.0,
                    fitness = fitness.function, penalty.function,
                    lower = lower_bound,
                    upper = upper_bound,
                    parallel = T,
                    monitor = F)

  end.time <- Sys.time()
  rt <- end.time - start.time

  num_solutions = length(genetic@solution)/(k*dims)
  if (num_solutions == 1) {
    solution <- matrix(genetic@solution, nrow = k, ncol = dims)
  } else {
    # if there is more than a single solution (they are identical, in ASW,
    # and must be close for centroids values)
    solution <- matrix(genetic@solution[1,], nrow = k, ncol = dims)
  }

  # calculates the average silhouette width
  which.dists <- apply(Rfast::dista(data, solution, "euclidean", square = TRUE), 1, which.min)

  asw <- cluster::silhouette(which.dists, d2)
  ch <- clusterCrit::intCriteria(as.matrix(data), which.dists, "Calinski_Harabasz")
  ci <- clusterCrit::intCriteria(as.matrix(data), which.dists, "C_index")
  di <- clusterCrit::intCriteria(as.matrix(data), which.dists, "Dunn")

  # builds the solution object
  solution.df <- as.data.frame(solution)
  colnames(solution.df) <- colnames(data)
  solution.df <- solution.df[with(solution.df, order(apply(solution.df, 1, sum))), ]

  # plot the results
  if (plot.internals) {
    #par(mfrow=c(1,2))
    plot(genetic, main = "Evolution")
    plot(asw, main = "ASW")
    #garbage <- dev.off()
  }

  object <- methods::new("gama",
                original.data = data,
                centers = solution.df,
                cluster = as.vector(which.dists),
                silhouette = summary(asw)$avg.width,
                calinski_harabasz = ch$calinski_harabasz,
                c_index = ci$c_index,
                dunn_index = di$dunn,
                runtime = round(rt, 2),
                call = call)

  print(object)

  # return an object of class 'gama'
  return (object)

}

setClass(Class = "gama",
         slots = c(original.data = "data.frame",
         centers = "data.frame",
         cluster = "vector",
         silhouette = "numeric",
         calinski_harabasz = "numeric",
         c_index = "numeric",
         dunn_index = "numeric",
         runtime = "ANY",
         call = "ANY"))

print.gama <- function(x, ...) {

  cat("\nDetails for the object of class 'gama':\n")

  cat("\nOriginal data (first rows):\n")
  print(head(x@original.data))
  cat("\nCluster partitions:\n")
  print(x@cluster)
  cat("\nCluster Centers:\n")
  print(x@centers)
  cat("\nAverage Silhouette Width index (ASW):\n")
  print(x@silhouette)
  cat("\nCalinski Harabasz index (CH):\n")
  print(x@calinski_harabasz)
  cat("\nC-Index (CI):\n")
  print(x@c_index)
  cat("\nDunn index (DI):\n")
  print(x@dunn_index)

  cat("\nCall:\n")
  print(x@call)
  cat("\nRuntime:\n")
  print(x@runtime)
}
setMethod("print", "gama", print.gama )

# summary.gama <- function(x = NULL, ...) {
#
#   print("to do...")
#
# }
# setMethod("summary", "gama", summary.gama)

# view.method = c("total.sum", "pca", "both")
plot <- function(x = NULL, ...) { UseMethod("plot.gama") }


plot.gama <- function(x = NULL, view.method = "pca", ... ) {

  dat <- x@original.data
  dat$clusters <- x@cluster


  if (view.method == "total.sum") {

    total.sum = apply(dat, 1, sum)

    dat$total.sum <- total.sum
    dat$observation <- 1:nrow(dat)
    g <- ggplot2::ggplot(dat, ggplot2::aes(x = observation, y = total.sum, color = factor(clusters))) +
      ggplot2::geom_point() +
      ggplot2::labs(color = "partition") +
      ggplot2::xlab("observation") +
      ggplot2::ylab("total sum of dimensions")
  } else if (view.method == "pca") {

    pca = prcomp(dat)

    dat$pc.1 <- pca$x[,"PC1"]
    dat$pc.2 <- pca$x[,"PC2"]
    g <- ggplot2::ggplot(dat, ggplot2::aes(x = pc.1, y = pc.2, color = factor(clusters))) +
      ggplot2::geom_point() +
      ggplot2::labs(color = "partition") +
      ggplot2::xlab("principal component 1") +
      ggplot2::ylab("principal component 2")
  }

  g <- g + ggplot2::ggtitle(paste("Gama partitions,", "view method = ", view.method, sep = " ")) + ggplot2::theme_minimal()
  plot(g)
}
setMethod("plot", "gama", plot.gama)

