library(GA)           # genetic algorithms
library(cluster)      # average silhouette width
library(ClusterR)     # distortion f(K)
library(Rfast)        # matrix fast calculations
library(ggplot2)      # graphics library
library(clusterCrit)  # validation criteria

source("R/fitness.R")
source("R/bestk.R")

# generates a initial population of centroids
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
  for(i in 1:object@popSize) {
    repeat {
      # which lines (1:k) of an individual in population (1:popSize) has load > 100
      m <- matrix(population[i,], ncol = dims, nrow = k)
      sums <- apply(m, 1, sum)
      w.sums <- which(sums > 100)
      
      # stop when all centroids of an individual have loads < 100
      if (length(w.sums) == 0){
        break
      } else {
        # subtracts the load by an proportional amount
        for (y in w.sums) {
          prop_over_all <- m[y,]/sums[y]
          population[i,] <- as.vector(m[y, ] - prop_over_all)
        }
      } # else
    }#repeat
  } # for
  
  return(population)
}


# tracing function
# monitor.f <- function(obj) {
#   now <- Sys.time()
#   if(now - checked.time > 5) {
#     .GlobalEnv$checked.time <- now
#     print(paste("Time elapsed: ", round(now - start.time, 2), "seconds", sep = " "))
#   } else return(FALSE) #????
# }

gama <- function(data, ...) UseMethod("gama")

gama.default <- function(data = NULL, k = NA, scale = FALSE, crossover.rate = 0.9, mutation.rate = 0.01, 
                         elitism = 0.05, pop.size = 25, generations = 100, seed.p = 42,
                         fitness.criterion = "ASW", 
                         penalty.function = NULL, 
                         plot.internals = TRUE, ...) {
  
  obj <- gama.clustering(data, k, scale, crossover.rate, mutation.rate, elitism, pop.size, generations, seed.p, 
                         fitness.criterion, penalty.function, plot.internals)
  #list(original.data = data, centroids=solution.df, cluster=as.vector(which.dists), asw.mean=summary(asw)$avg.width))
  obj$call <- match.call()
  class(obj) <- "gama"
  obj 
}

print.gama <- function(x, ...) {
  
  cat("\nResults:\n\n")
  cat("\nOriginal data:\n")
  print(x$original.data)
  cat("\nCluster partitions:\n")
  print(x$cluster)
  cat("\nCluster Centroids:\n")
  print(x$centroids)
  cat("\nAverage Silhouette Width:\n")
  print(x$silhouette)
  cat("\nCalinski_Harabasz:\n")
  print(x$calinski_harabasz)
  cat("\nC Index:\n")
  print(x$c_index)
  cat("\nDunn Index:\n")
  print(x$dunn)
  
  cat("Call:\n")
  print(x$call)
  cat("Runtime:\n")
  print(x$runtime)
}

# execute the clustering process guided by a defined criteria
gama.clustering <- function(data = NULL, k = NA, scale = FALSE, crossover.rate = 0.9, mutation.rate = 0.01, 
                            elitism = 0.05, pop.size = 25, generations = 100, seed.p = 42,
                            fitness.criterion = "ASW", penalty.function = NULL, plot.internals = TRUE) {
  
  # choose the correct fitnesse function
  fitness.function <- switch (fitness.criterion, 
                              "ASW" = fitness.asw, 
                              "CH" = fitness.ch, 
                              "CI" = fitness.ci,
                              "DI" = fitness.di,
                              fitness.asw)
  
  if (scale) {
    data <- scale(data, center = TRUE, scale = TRUE)
  }
  
  # uses distortion f(K) to choose the best k estimative
  if (is.na(k)) {
    
    # cm <- citation("ClusterR")
    # print(paste("Estimating k by using distortion f(K) method...", format(cm, style = "text")))
    #print("")
    
    opt = Optimal_Clusters_KMeans(data, 
                                  max_clusters = 10,
                                  plot_clusters = F,
                                  verbose = F,
                                  criterion = 'distortion_fK', 
                                  fK_threshold = 0.85,
                                  initializer = 'kmeans++', 
                                  tol_optimal_init = 0.2)
    
    k <- best.k(dataset = data)
    print(paste("best k suggestion = ", k))
    
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
  rm(d)
  gc()
  
  lowers <- apply(data, 2, min)
  uppers <- apply(data, 2, max)
  
  lower_bound <- unlist(lapply(lowers, function (x) { rep(x, k) } ))
  upper_bound <- unlist(lapply(uppers, function (x) { rep(x, k) } ))
  #lower_bound <-  c(rep(min(data[, 1]), k), rep(min(data[, 2]), k), rep(min(data[, 3]), k), rep(min(data[, 4]), k))
  #upper_bound <-  c(rep(max(data[, 1]), k), rep(max(data[, 2]), k), rep(max(data[, 3]), k), rep(max(data[, 4]), k))
  
  
  # call GA functions 
  # cm <- citation("GA")
  # print(paste("Starting genetic seach for centroids, by using GA Algorithms...", format(cm, style = "text")))
  #print("")
  start.time <- Sys.time()
  
  .GlobalEnv$start.time <- start.time
  .GlobalEnv$checked.time <- start.time
  
  #if (missing(penalty.function)) penalty.function <- NA
  
  genetic <- ga(type = "real-valued", 
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
                parallel = F,
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
  which.dists <- apply(dista(data, solution, "euclidean", square = TRUE), 1, which.min)
  
  asw <- silhouette(which.dists, d2)
  ch <- intCriteria(as.matrix(data), which.dists, "Calinski_Harabasz")
  ci <- intCriteria(as.matrix(data), which.dists, "C_index")
  di <- intCriteria(as.matrix(data), which.dists, "Dunn")
  
  
  print(paste("Clustering process completed in:", rt, sep = " "))
  print(paste("Acceptable (distinct) solutions: ",  length(genetic@solution)/(k*dims), sep = " "))
  print(paste("Average Silhouette Width (ASW): ", round(summary(asw)$avg.width, 4), sep = " "))
  print(paste("Calinski Harabasz (CH): ", round(ch$calinski_harabasz, 4), sep = " "))
  print(paste("C-Index (CI): ", round(ci$c_index, 4), sep = " "))
  print(paste("Dunn Index (DI): ", round(di$dunn, 4), sep = " "))
  
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
    
    # cm <- citation("cluster")
    # print(paste("Average silhouette width package...", format(cm, style = "text")))
  }
  
  # setClass("gama", slots=list(original.data = "data.frame", centroids="data.frame", cluster="vector", asw.width="numeric"))
  # return(new ("gama", "original.data" = data, "centroids" = solution.df, "cluster" = as.vector(which.dists), "asw.width" = summary(asw)$avg.width))
  return(list(original.data = data, 
              centroids = solution.df, 
              cluster = as.vector(which.dists), 
              silhouette = summary(asw)$avg.width, 
              calinski_harabasz = ch$calinski_harabasz, 
              c_index = ci$c_index,
              dunn_index = di$dunn,
              runtime = rt))
}

# view.method = c("total.sum", "pca", "both")
plot.gama <- function(x = NULL, view.method = "pca") {
  
  dat <- x$original.data
  dat$clusters <- x$cluster
  
  
  if (view.method == "total.sum") {
    
    total.sum = apply(dat, 1, sum)
    
    dat$total.sum <- total.sum
    dat$observation <- 1:nrow(dat)
    g <- ggplot(dat, aes(x = observation, y = total.sum, color = factor(clusters))) + 
      geom_point() + 
      labs(color = "partition") +
      xlab("observation") +
      ylab("total sum of dimensions")
  } else if (view.method == "pca") {
    
    pca = prcomp(dat)
    
    dat$pc.1 <- pca$x[,"PC1"]
    dat$pc.2 <- pca$x[,"PC2"]
    g <- ggplot(dat, aes(x = pc.1, y = pc.2, color = factor(clusters))) + 
      geom_point() + 
      labs(color = "partition") +
      xlab("principal component 1") +
      ylab("principal component 2")
  }
  
  g <- g + ggtitle("Cluster") + theme_minimal()
  plot(g)
}