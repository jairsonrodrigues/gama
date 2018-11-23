source("helpers.r")
data("als.cpu.big.1")

# a specific function to calculate penalties for CPU datasets 
# whose does not allow the sum of loads above 100%
my.cpu.penality <- function(m.individual,...) {
  
  penalty <- 0
  
  # counts the how many centroids results in overflow (inequality > 100)
  sums <- apply(m.individual, 1, sum)
  overflow <- which(sums > 100)
  num_constraints = length(overflow)
  
  # if there are overflows, subtract each dimension 
  # by the maximum proportion of the excess x the number of overflows
  if (num_constraints > 0) { 
     penalty <- num_constraints * max(abs(sums[overflow] -100)/sums[overflow]) 
  }
  
  return (penalty)
}

gama.p <- gama(als.cpu.big.1, plot.internals = T, generations = 25, seed.p = 17, penality.function = my.cpu.penality)
gama.nop <- gama(als.cpu.big.1, plot.internals = T, generations = 25, seed.p = 17)
plot(gama.p, view.method = "total.sum")
plot(gama.nop, view.method = "total.sum")


# iris tests
data(iris)
i <- iris[,1:4]

k = 2
gama.i <- gama(i, k = 2, plot.internals = F)
km.i <- kmeans(i, centers = 2)
plotcluster(i, gama.i$cluster, main = "gama")
plotcluster(i, km.i$cluster, main = "kmeans")

print(gama.i$asw.mean)
asw <- silhouette(km.i$cluster, d2)
print(summary(asw)$avg.width)

k = 3
gama.i <- gama(i, k = 3, plot.internals = F)
km.i <- kmeans(i, centers = 3)
plotcluster(i, gama.i$cluster)
plotcluster(i, km.i$cluster)

print(gama.i$asw.mean)
asw <- silhouette(km.i$cluster, d2)
print(summary(asw)$avg.width)

k = 4
gama.i <- gama(i, k = 4, plot.internals = F)
km.i <- kmeans(i, centers = 4)
plotcluster(i, gama.i$cluster)
plotcluster(i, km.i$cluster)

print(gama.i$asw.mean)
asw <- silhouette(km.i$cluster, d2)
print(summary(asw)$avg.width)

k = 5
gama.i <- gama(i, k = 5, plot.internals = F)
km.i <- kmeans(i, centers = 5)
plotcluster(i, gama.i$cluster)
plotcluster(i, km.i$cluster)

print(gama.i$asw.mean)
asw <- silhouette(km.i$cluster, d2)
print(summary(asw)$avg.width)

