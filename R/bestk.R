# This function estimates the best k value for the number of partitions
# the dataset should be segmented.
gama.how.many.k <- function(dataset = NULL, method = "minimal") {

    # --- arguments validation --- #
    Check <- ArgumentCheck::newArgCheck()

    if (is.null(dataset))
      ArgumentCheck::addError(
        msg = "'dataset' can not be NULL",
        argcheck = Check
      )

    if (class(dataset) != 'data.frame')
      ArgumentCheck::addError(
        msg = "'dataset' must be a data.frame object.",
        argcheck = Check
      )

    if (!(method %in% c('minimal', 'broad')))
      ArgumentCheck::addError(
        msg = "'method' must be one of the values: 'minimal' or 'broad'.",
        argcheck = Check
      )

    ArgumentCheck::finishArgCheck(Check)

    # --- final of arguments validation --- #

    best.k <- -1

    if (method == "minimal") {
      best.k <- best.k.minimal(dataset)
    } else if (method == "broad") {
      best.k <- best.k.broad(dataset)
    }

    return (best.k)
}

# compute 24 methods by using NbClust package
best.k.broad <- function(dataset = NULL) {

  cat("Looking for the best k value by 'broad' method...\n")
  indices <- c("kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin", "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "frey", "mcclain", "dunn", "sdindex", "sdbw")

  nc <- c()
  count = 1
  # nc will contais the best k value suggestion for each index
  for (i in indices) {
    objNbClust <- NbClust::NbClust(dataset, min.nc = 2, max.nc = 10, method = "kmeans", index = i)
    nc[count] <- as.integer(objNbClust$Best.nc["Number_clusters"])
    count = count + 1
  }

  nc <- nc[which(nc >=2)]
  uniques <- sort(unique(nc))
  best_k <- uniques[which.max(sapply(uniques, function(x) { length(which(nc == x)) } ))]
  cat(paste("There was found", length(uniques), "suggestions for k = [", toString(uniques), "].\n", sep = " "))
  cat(paste("Best k suggestion by majority voting = ", best_k, ".\n\n", sep =""))
  return (best_k)
}

# compute the WCSSE (elbow method)
best.k.minimal <- function (dataset = NULL) {

  cat("Looking for the best k value by 'minimal' method...\n")

  k.max <- 10
  set.seed(42)
  wss <- sapply(1:k.max, function(k){
            kmeans(dataset, centers = k, algorithm = "Hartigan-Wong",  nstart=25, iter.max = 100)$tot.withinss
          })

  best.k <- where.is.knee(wss)

  plot(1:k.max, wss, type="b",
            xlab="Number of Clusters",
            ylab="Within Cluster Sum of Squares Error",
       main = paste("Elbow graphic, best k = ", best.k, sep = ""))

  abline(v = best.k, lwd=1, lty=2, col = "red")

  cat(paste("Best k suggestion by using second derivative aproximation over \n'Within-Cluster Sum of Squares Error' (WCSSE) graphic = ", best.k, ". ", sep = ""))
  cat("See 'elbow' graphic for details.\n\n")

  return(best.k)

}

# Calculates an aproximation of the second derivative of a set of points
# the maximum second derivative will be a good choice for the inflexion point (the elbow or knee)
# https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve
# https://raghavan.usc.edu/papers/kneedle-simplex11.pdf (Finding a "Kneedle" in a Haystack:
# Detecting Knee Points in System Behavior)
where.is.knee <- function(points = NULL) {

  lower.limit <- 2
  upper.limit <- length(points) -1

  second.derivative <- sapply(lower.limit:upper.limit, function(i) { points[i+1] + points[i-1] - (2 * points[i]) })

  w.max <- which.max(second.derivative)
  best.k <- w.max +1

  return(best.k)

}
