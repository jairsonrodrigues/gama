# identifies if the solution has noise
has.noise.f <- function(sol) {
  return (is.element(0, as.integer(sol)))
}

# return the indices of noise elements in solution
which.noise.f <- function(sol) {
  return (which(as.integer(sol) == 0))
}


# convnerts labels choosen by algorithm by labels in order of centroids or mean values
convert.real.f <- function(sol, centroid) {
  
  df <- as.data.frame(centroid)
  tmp = as.integer(sol)
  for (i in 1:length(df$cluster.label)) {
    tmp[which(sol == df$cluster.label[i])] <- df$cluster.order[i]
  }
  
  return (tmp)
}


pertinence.f <- function(theBest, solutions, centroids, desc, toLatex = FALSE) {
  
  cat(desc, " (ASW: ", theBest[2], ", DBCV: ", theBest[4], ")", sep = "")
  cat("\n\n")
  
  sol <- solutions[theBest[1],]
  algo <- cpu.all$algo
  
  # if has noise, remove the partition from solution
  if (has.noise.f(sol)) {
    noise_indices <- which.noise.f(sol)
    sol <- sol[-noise_indices]
    algo <- cpu.all$algo[-noise_indices]
  }
  
  theCentroid <- centroids[[theBest[1]]]
  sol_converted <- convert.real.f(sol, theCentroid)
  
  ct <- table(algo, as.integer(sol_converted))
  
  #med$cluster <- c(1:best_k)
  for (i in 1:nrow(ct)) {
    major = max(ct[i,])/sum(ct[i,])
    w.major = which.max(ct[i,])
    all = which(ct[i,] >0)
    if (toLatex) {
       cat(" &", round(major  * 100, 2) , "\\% $|" , w.major , "| ^\\{", paste(c(all[-which(all == w.major)]), collapse=", ") , "\\}$\n", sep = " ")
     } else {
       cat(row.names(ct)[i], ": ", round(major  * 100, 2) , "% [" , w.major , "], distortion set = {", paste(c(all[-which(all == w.major)]), collapse=", ") , "}\n", sep = "")
     }
  }
  
  cat("\n\n")
  
}

load(file = "RData/solutions.RData")
load(file = "RData/criteria.RData")
cpu <- read.csv("data/cpu-large.csv")
cpu.all <- read.csv("data/all_exp2_large.csv")

# k-means
pertinence.f(theBestOnes$kmeans, km_sol, km_centroids, "k-means, large", T)
pertinence.f(theBestOnes$gama, gama_sol, gama_centroids, "gama, large", T)
pertinence.f(theBestOnes$dbscan, dbscan_sol, dbscan_centroids, "dbscan, large", T)
pertinence.f(theBestOnes$nkhga, nkhga_sol, nkhga_centroids, "nkhga, large", T)
