library(cluster)
library(Rfast)
library(fpc)
library(factoextra)

source("util.r")

plot.f <- function(sol, d, algo = "not specified") {
  
  cmd <- cmdscale(d, eig = TRUE)
  df <- data.frame("dc1" = cmd$points[,1], "dc2" = cmd$points[,2], "cluster" = sol)
  
  noise_df <- NULL
  if (has.noise.f(sol)) {
    w.noise <- which.noise.f(sol)
    df <- df[-c(w.noise),]
    noise <- cmd$points[w.noise,]
    if (is.null(nrow(noise))) # special case: only a single noise in data set
      noise_df <- data.frame("dc1" = noise[1], "dc2" = noise[2], "cluster" = sol[w.noise])
    else 
      noise_df <- data.frame("dc1" = noise[,1], "dc2" = noise[,2], "cluster" = sol[w.noise])
  }
  
  p1=ggplot(df, aes(x=dc2, y=dc1, color=factor(cluster), shape = factor(cluster)))
  p1= p1 + geom_point(size = 3.5) + theme_bw(base_size = 18) + 
    xlab("component 2") +
    ylab("component 1") +
    #ggtitle(algo) +  
    guides(color=guide_legend(title="CPU load:"), 
           shape=guide_legend(title="CPU load:")) +
    theme(legend.position="top")  
  
  if (!is.null(noise_df)) {
    p1 = p1 + geom_point(data = noise_df, color="black", size = 12, pch = '*') +
      labs(caption = "* = noise") 
  }
  
  p1
}

cpu <- read.csv("data/cpu-large.csv")
all <- read.csv("data/all_exp2_large.csv")
d <- dist(cpu)

# read solutions
load(file = "RData/solutions.RData")
load(file = "RData/criteria.RData")

# convert solutions as vector of doubles 
# take the best ones for each clustering algorithm (ASW)
sol_km <- as.double(km_sol[theBestOnes$kmeans[1],])
sol_gama <- as.double(gama_sol[theBestOnes$gama[1],])
sol_dbscan <- as.double(dbscan_sol[theBestOnes$dbscan[1],]) #+1 # to allow printable noise points (zero values)
sol_nkhga <- as.double(nkhga_sol[theBestOnes$nkhga[1],]) #+ 1 # to allow printable noise points (zero values)


sol_km <- convert.real.f(sol_km, centroid = km_centroids[[theBestOnes$kmeans[1]]])
sol_gama <- convert.real.f(sol_gama, centroid = gama_centroids[[theBestOnes$gama[1]]])
sol_dbscan <- convert.real.f(sol_dbscan, centroid = dbscan_centroids[[theBestOnes$dbscan[1]]])
sol_nkhga <- convert.real.f(sol_nkhga, centroid = nkhga_centroids[[theBestOnes$nkhga[1]]])

pdf(file = "images/kmeans_shape.pdf")
  plot.f(sol_km, d, "k-means")
  garbage <- dev.off()

pdf(file = "images/gama_shape.pdf")
  plot.f(sol_gama, d, "GAMA")
  garbage <- dev.off()

pdf(file = "images/dbscan_shape.pdf")
  plot.f(sol_dbscan, d, "DBSCAN")
  garbage <- dev.off()

pdf(file = "images/nkhga_shape.pdf")
  plot.f(sol_nkhga, d, "NKHGA")
  garbage <- dev.off()