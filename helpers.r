library(GA)       # genetic algorithms
library(cluster)  # average silhouette width
library(ClusterR) # distortion f(K)
library(Rfast)    # matrix fast calculations
library(ggplot2)  # graphics library

# evaluate individuals according Average Silhouette Width Criterion
fitness.asw <- function(individual, penality = TRUE) {
  
  dims <- length(individual)/k
  
  m.individual <- matrix(individual, nrow = k, ncol = dims)
  which.dists <- apply(dista(data, m.individual, "euclidean", square = TRUE), 1, which.min)
  
  # calculate the average silhouette width
  asw <- silhouette(which.dists, d2)

  # try summarize the silhouette, returns (0) if error
  fitness.value <- tryCatch(summary(asw)$avg.width, error = function (e) { return (0)})
  
  # optionally, some inequalities may be resolved by applying a penality
  if (penality) {
    sums <- apply(m.individual, 1, sum)
    overflow <- which(sums > 100)
    num_constraints = length(overflow)

    if (num_constraints > 0) { 
      penalty <- num_constraints * max(abs(sums[overflow] -100)/sums[overflow]) 
      fitness.value <- fitness.value - penalty
    }
  }
  
  return (fitness.value)
}

# generates a initial population of centroids
pop.f <- function(object) {
  
  lower <- object@lower
  upper <- object@upper
  nvars <- length(lower)
  
  population <- matrix(as.double(NA), nrow = object@popSize, ncol = nvars)
  
  # generate the initial population based on uniform distribution
  for(j in 1:nvars) { 
    population[,j] <- runif(object@popSize, lower[j], upper[j]) 
  }
  
  # checks impossible individuals which sum(loads) > 100
  # subtracts the load value until it reach the possible range
  for(i in 1:object@popSize) {
    repeat {
      # which lines (1:k) of an individual in population (1:popSize) has load > 100
      m <- matrix(population[i,], ncol = 4, nrow = k)
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

# execute the clustering process guided by a defined criteria
gama <- function(data = NULL, k = NA, crossover.rate = 0.9, mutation.rate = 0.01, 
                 elitism = 0.05, pop.size = 25, generations = 100, seed.p = 42,
                 fitness.function = fitness.asw, plot.results = TRUE) {
  
  
  # uses distortion f(K) to choose the best k estimative
  if (is.na(k)) {
    
    cm <- citation("ClusterR")
    print(paste("Estimating k by using distortion f(K) method...", format(cm, style = "text")))
    #print("")
  
    opt = Optimal_Clusters_KMeans(data, 
                                  max_clusters = 10,
                                  plot_clusters = F,
                                  verbose = F,
                                  criterion = 'distortion_fK', 
                                  fK_threshold = 0.85,
                                  initializer = 'kmeans++', 
                                  tol_optimal_init = 0.2)
    
    k = as.integer(which.min(opt))
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
  
  lower_bound <-  c(rep(min(data[, 1]), k), rep(min(data[, 2]), k), rep(min(data[, 3]), k), rep(min(data[, 4]), k))
  upper_bound <-  c(rep(max(data[, 1]), k), rep(max(data[, 2]), k), rep(max(data[, 3]), k), rep(max(data[, 4]), k))
  
  
  # call GA functions 
  cm <- citation("GA")
  print(paste("Starting genetic seach for centroids, by using GA Algorithms...", format(cm, style = "text")))
  #print("")
  start.time <- Sys.time()
  
  .GlobalEnv$start.time <- start.time
  .GlobalEnv$checked.time <- start.time
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
                maxFitness = 1.0, 
                fitness = fitness.function,
                lower = lower_bound,
                upper = upper_bound,
                parallel = T,
                monitor = F)
  
  end.time <- Sys.time()
  
  num_solutions = length(genetic@solution)/(k*dims)
  if (num_solutions == 1) { 
    solution <- matrix(genetic@solution,nrow = k,ncol = dims)
  } else {
    # if there is more than a single solution (they are identical, in ASW, 
    # and must be close for centroids values)
    solution <- matrix(genetic@solution[1,], nrow = k, ncol = dims)
  }
  
  # calculates the average silhouette width
  which.dists <- apply(dista(data, solution, "euclidean", square = TRUE), 1, which.min)
  asw <- silhouette(which.dists, d2)
  
  print(paste("Clustering process completed in:", round(end.time - start.time, 2), "seconds", sep = " "))
  print(paste("Acceptable (identical) solutions: ",  length(genetic@solution)/(k*dims), sep = " "))
  print(paste("Average Silhouette Width (ASW): ", round(summary(asw)$avg.width, 2), sep = " "))
  
  # builds the solution object
  solution.df <- as.data.frame(solution)
  colnames(solution.df) <- colnames(data)
  solution.df <- solution.df[with(solution.df, order(apply(solution.df, 1, sum))), ]
  
  # plot the results
  if (plot.results) {
    plot(genetic)
    plot(asw)
  }
  
  setClass("gama", slots=list(original.data = "data.frame", centroids="data.frame", cluster="vector", asw.width="numeric"))
  return(new ("gama", "original.data" = data, "centroids" = solution.df, "cluster" = as.vector(which.dists), "asw.width" = summary(asw)$avg.width))
  
}

# view.method = c("total.sum", "pca", "both")
plot.clusters <- function(gamaResObject = NULL, view.method = "total.sum") {
  
  dat <- gamaResObject@original.data
  dat$clusters <- gamaResObject@cluster
  
  
  if (view.method == "total.sum") {
    
      total.sum = apply(dat, 1, sum)
      
      dat$total.sum <- total.sum
      dat$observation <- 1:nrow(dat)
      g <- ggplot(dat, aes(x = observation, y = total.sum, color = factor(clusters))) + 
           geom_point() + 
           labs(color = "partition") +
           theme_minimal()
  } else if (view.method == "pca") {
    
      pca = prcomp(dat)
      
      dat$pc.1 <- pca$x[,"PC1"]
      dat$pc.2 <- pca$x[,"PC2"]
      g <- ggplot(dat, aes(x = pc.1, y = pc.2, color = factor(clusters))) + 
          geom_point() + 
          labs(color = "partition") +
          theme_minimal()
  }
  
  plot(g)
}