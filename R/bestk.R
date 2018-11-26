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

# Finds suggestions for best k based on 
# Lampros Mouselimis (2018). ClusterR: Gaussian Mixture Models, K-Means, Mini-Batch-Kmeans and K-Medoids
# Clustering. R package version 1.1.0. https://CRAN.R-project.org/package=ClusterR
best.k <- function(dataset = NULL) {
  
  distortion.fk.threshold = 0.85
  criteria <- c("distortion_fK", "variance_explained", "WCSSE", "dissimilarity", "silhouette", "AIC", "BIC", "Adjusted_Rsquared")
  
  bests <- c()
  count <- 1
  for (criterion in criteria) {
    #pdf(g_df)
    opt = Optimal_Clusters_KMeans(dataset, 
                                  max_clusters = 10,
                                  #num_init = 10,
                                  plot_clusters = F,
                                  verbose = F,
                                  criterion = criterion, 
                                  fK_threshold = distortion.fk.threshold,
                                  initializer = 'kmeans++', 
                                  tol_optimal_init = 0.2)
    
    #garbage <- dev.off()
    
    if (criterion == "distortion_fK") {
      
      #best_k = min(which(opt < distortion.fk.threshold))
      best_k = which.min(opt)
      
    } else if (criterion == "silhouette") {
      
      best_k <- which.max(opt)
      
    } else {
      best_k <- where.is.knee(opt)
    }
    
    #print(best_k)
    bests[count] <- best_k
    count <- count + 1
  }
  
  uniques <- unique(bests)
  best_k <- uniques[which.max(sapply(uniques, function(x) { length(which(bests == x)) } ))]
  
  return (best_k)
}