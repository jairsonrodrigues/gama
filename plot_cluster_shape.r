library(cluster)
library(Rfast)
library(fpc)
library(factoextra)
library(colorspace)

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
  
  pal.c = c(90,40)
  pal.l=c(33,93)
  pal.power=c(1/5, 1.5)
  cols <- rev(heat_hcl(max(df$cluster), c = pal.c, l=pal.l, power=pal.power))
  
  p1=ggplot(df, aes(x=dc2, y=dc1, color = factor(cluster), shape = factor(cluster)))
  p1= p1 + geom_point(size = 3.5, alpha = 0.1) +  
    scale_color_manual(values = cols) +
    xlab("component 2") +
    ylab("component 1") +
    #ggtitle(algo) +  
    theme_minimal(base_size = 24) +
    guides(color=guide_legend(title="CPU load:"), 
           shape=guide_legend(title="CPU load:")) +
    
    theme(legend.position="bottom", 
          axis.ticks.length=unit(0.15,"cm"),
          axis.ticks = element_line(size = 0.25),
          axis.line = element_line(colour = "black"))  
  
  if (!is.null(noise_df)) {
    p1 = p1 + geom_point(data = noise_df, color="black", size = 12, pch = '*') +
      labs(caption = "* = noise") 
  }
  
  p1
}

plot.f.alpha <- function(sol, d, algo = "not specified") {
  
  cmd <- cmdscale(d, eig = TRUE)
  df <- data.frame("dc1" = cmd$points[,1], "dc2" = cmd$points[,2], "cluster" = sol)
  
  # calculate the alpha values (transparency)
  t1<-sapply(1:5, function(x) { df$dc1[df$cluster == x] } )
  t2<-sapply(1:5, function(x) { df$dc2[df$cluster == x] } )
  dfc <- data.frame(dc1 = unlist(lapply(t1, mean)), dc2 = unlist(lapply(t2, mean)))
  
  # alphas <- sapply(1:5, function(y) { 
  #   tmp <- 1/apply(df[df$cluster == y,1:2], 1,  function(x) sqrt(sum((x - dfc[y,]) ^ 2)) )
  #   tmp <- round(tmp/max(tmp), 1)
  #   tmp
  # })
  
  alphas <- sapply(1:5, function(y) { 
    tmp <- apply(df[df$cluster == y,1:2], 1,  function(x) sqrt(sum((x - dfc[y,]) ^ 2)) )
    tmp <- round(tmp/max(tmp), 1)
    tmp < 1- tmp
    tmp
  })
  
  df$alphas <- unlist(alphas)
  
  df$alphas <- sapply(df$alphas, function(x) { 
    if (x >= 0.66) { 
      ret = 1; 
    } else if (x<=0.33) { 
      ret = 0.1 
      } else { 
        ret = 0.5 
      } 
    ret
  })
  
  
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
  
  pal.c = c(90,40)
  pal.l=c(33,93)
  pal.power=c(1/5, 1.5)
  cols <- rev(heat_hcl(max(df$cluster), c = pal.c, l=pal.l, power=pal.power))
  
  p1=ggplot(df, aes(x=dc2, y=dc1, color = factor(cluster), shape = factor(cluster)))
  #p1=ggplot(df, aes(x=dc2, y=dc1, color = factor(cluster)))
  p1= p1 + geom_point(aes(alpha = alphas), size = 3.5) +  
    scale_color_manual(values = cols) +
    xlab("component 2") +
    ylab("component 1") +
    #ggtitle(algo) +  
    theme_minimal(base_size = 24) +
    guides(color=guide_legend(title="CPU load:"), 
           shape=guide_legend(title="CPU load:")) +
    
    theme(legend.position="bottom", 
          axis.ticks.length=unit(0.15,"cm"),
          axis.ticks = element_line(size = 0.25),
          axis.line = element_line(colour = "black"))  
  
  # if (!is.null(noise_df)) {
  #   p1 = p1 + geom_point(data = noise_df, color="black", size = 12, pch = '*') +
  #     labs(caption = "* = noise") 
  # }
  
  p1 = p1 + geom_point(data = dfc, color="black", size = 3.5, pch = 17)
  #p1 = p1 + scale_alpha("alphas" )
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

# pdf(file = "images/kmeans_shape.pdf")
#   plot.f(sol_km, d, "k-means")
#   garbage <- dev.off()
# 
# pdf(file = "images/gama_shape.pdf")
#   plot.f(sol_gama, d, "GAMA")
#   garbage <- dev.off()
# 
# pdf(file = "images/dbscan_shape.pdf")
#   plot.f(sol_dbscan, d, "DBSCAN")
#   garbage <- dev.off()
# 
# pdf(file = "images/nkhga_shape.pdf")
#   plot.f(sol_nkhga, d, "NKHGA")
#   garbage <- dev.off()

plot.f.alpha(sol_km, d, "k-means")