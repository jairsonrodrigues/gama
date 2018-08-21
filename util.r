#
# Utils code to support project
# Author: [Jairson Rodrigues](jairsonrodrigues@gmail.com)
#

# returns TRUE if the solution set contains something
has.noise.f <- function(sol) {
  return (is.element(0, as.integer(sol)))
}

# return the indices of noise elements in solution
which.noise.f <- function(sol) {
  return (which(as.integer(sol) == 0))
}


# convert labels choosen by algorithm by labels in order of centroids or mean values
convert.real.f <- function(sol, centroid) {
  
  df <- as.data.frame(centroid)
  tmp = as.integer(sol)
  for (i in 1:length(df$cluster.label)) {
    tmp[which(sol == df$cluster.label[i])] <- df$cluster.order[i]
  }
  
  return (tmp)
}

# the inverse of intersect 
outersect.f <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}