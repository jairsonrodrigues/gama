#
# Utils code to support project
# Author: [Jairson Rodrigues](jairsonrodrigues@gmail.com)
#

# round up to the nearest 10 (or 100 or 1000 ...)
roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

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