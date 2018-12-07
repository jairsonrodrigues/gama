# to avoid no visible binding for global variable
## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("k", "dat", "data", "d2"))

# calculates an aproximation of the second derivative of a set of points
# the maximum second derivative will be a good choice for the inflexion point (the elbow or knee)
# https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve
# https://raghavan.usc.edu/papers/kneedle-simplex11.pdf (Finding a “Kneedle” in a Haystack:
# Detecting Knee Points in System Behavior)
where.is.knee <- function(dataset = NULL) {

  lower.limit <- 2
  upper.limit <- length(dataset) -1

  second.derivative <- sapply(lower.limit:upper.limit, function(i) { dataset[i+1] + dataset[i-1] - 2 * dataset[i] })

  w.max <- which.max(second.derivative)

  return(w.max +1)

}




# \name{best.k_old}
# \alias{best.k_old}
# \title{Suggests a value for best k.}
# \description{
#   Applies seven methods to discover the best value for k (the optmal number of partitions), based on criteria:
#
#     \enumerate{
#       \item Variance Explained: the sum of the within-cluster-sum-of-squares of all clusters divided by the total sum of squares
#       \item WCSSE: the sum of the within-cluster-sum-of-squares of all clusters
#       \item Dissimilarity: the average intra-cluster dissimilarity of all clusters (the distance metric defaults to euclidean)
#       \item Distortion f(K) : this criterion is based on the following paper, "Selection of K in K-means clustering" (https://www.ee.columbia.edu/~dpwe/papers/PhamDN05-kmeans.pdf)
#       \item AIC: the Akaike information criterion
#       \item BIC: the Bayesian information criterion
#       \item Adjusted R Squared : the adjusted R^2 statistic
#     }
#
#   The function finds all unique values suggested by the methods, ordered in ascending manner. After that, the function returns the best k by majority voting.
#
# }
# \usage{
#   best.k_old(dataset = NULL)
# }
# \arguments{
#   \item{dataset}{a dataframe containing numerical data to be separeted in k partitions.}
# }
# \seealso{
#   \code{\link{where.is.knee}}.
# }
# \references{
#   Lampros Mouselimis (2018). ClusterR: Gaussian Mixture Models, K-Means,
#   Mini-Batch-Kmeans, K-Medoids and Affinity Propagation Clustering. R package version
#   1.1.6. https://CRAN.R-project.org/package=ClusterR
#
# }
# \examples{
#   ## CPU metrics for Alternating Least Squares
#   data("cpu.als")
#
#   ## discovers and prints the best value for the
#   ## optimal number of partitions by majority voting
#   k_values <- best.k(cpu.als)
#   print(k_values)
# }
#
# \keyword{file}



# Finds suggestions for best k based on
# Lampros Mouselimis (2018). ClusterR: Gaussian Mixture Models, K-Means, Mini-Batch-Kmeans and K-Medoids
# Clustering. R package version 1.1.0. https://CRAN.R-project.org/package=ClusterR
# best.k.old <- function(dataset = NULL) {
#
#   distortion.fk.threshold = 0.85
#   criteria <- c("distortion_fK", "variance_explained", "WCSSE", "dissimilarity", "silhouette", "AIC", "BIC", "Adjusted_Rsquared")
#   #criteria <- c("distortion_fK", "variance_explained", "WCSSE", "dissimilarity", "AIC", "BIC", "Adjusted_Rsquared")
#
#   bests <- c()
#   count <- 1
#   for (criterion in criteria) {
#     opt = ClusterR::Optimal_Clusters_KMeans(dataset,
#                                   max_clusters = 10,
#                                   #num_init = 10,
#                                   plot_clusters = F,
#                                   verbose = F,
#                                   criterion = criterion,
#                                   fK_threshold = distortion.fk.threshold,
#                                   initializer = 'kmeans++',
#                                   tol_optimal_init = 0.2)
#
#     if (criterion == "distortion_fK") {
#
#       #best_k = min(which(opt < distortion.fk.threshold))
#       best_k = which.min(opt)
#
#     } else if (criterion == "silhouette") {
#
#       best_k <- which.max(opt)
#
#     } else {
#       best_k <- where.is.knee(opt)
#     }
#
#     #print(best_k)
#     bests[count] <- best_k
#     count <- count + 1
#   }
#
#   uniques <- unique(bests)
#   best_k <- uniques[which.max(sapply(uniques, function(x) { length(which(bests == x)) } ))]
#
#   return (best_k)
# }

best.k <- function(dataset = NULL) {

  cat("Looking for the best k value...\n")
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
