#
# Source code to plot heatmaps and cluster partitions over the nodes graphics
# Author: [Jairson Rodrigues](jairsonrodrigues@gmail.com)
#  

# import util functions of the project
source("util.r")

# import libraries
library(ggplot2)
# a more convenient space color
library(colorspace) 

# load variables needed to cluster analysis
load(file = "RData/solutions.RData")
load(file = "RData/criteria.RData")
load(file = "RData/metrics.RData")

# plot the CPU load heatmap for a DML algorithm
plot.heatmap.load.f <- function(dataset, algo) {
  
    target <- which(dataset$algo == algo)
    
    #dataset$total[target] <- dataset$user[target] + dataset$system[target] + dataset$iowait[target] + dataset$softirq[target]
    #target <- unique(dataset$total)
    nodes <- as.integer(sapply(dataset$hostname[target], function(x) { substr(x, 15, 15) })) + 1
    load <- dataset$user[target] + dataset$system[target] + dataset$iowait[target] + dataset$softirq[target]
    #load <- dataset$user[target] 
    #load <- dataset$system[target] 
    #load <- round(dataset$iowait[target], 2)
    #load <- round(dataset$softirq[target], 2)
    
    moment <- round((dataset$timestamp[target] - min(dataset$timestamp[target])) + 1, 0)
    # limite = 600
    # nodes <- nodes[1:limite]
    # moment <- moment[1:limite]
    # load <- load[1:limite]

    # choose length(unique(target) different values for colors in heat.colors
    heat = data.frame("nodes" = nodes, "load" = load, "moment" = moment)
    
    # uses heat_cld palette from colorspace package
    #heat_chl(12, c = c(80,30), l=c(30,90), power=c(1/5, 1.5))
    #cols <- rev(heat_hcl(length(unique(target)), alpha = 1))
    pal.c = c(90,40)
    pal.l=c(33,93)
    pal.power=c(1/5, 1.5)
    cols <- rev(heat_hcl(length(unique(target)), c = pal.c, l=pal.l, power=pal.power))
    cols <- rev(heat_hcl(20, c = pal.c, l=pal.l, power=pal.power))
    
    # set the x-axis breaks only to odd indices in moments vector (to save space in plot)
    x_breaks <- seq(1,max(moment),max(10, roundUpNice(max(moment/10))))-1
    #leg_breaks <- c(min(moment), roundUpNice(max(moment)/2), max(moment))
    
    c <- ggplot(heat, aes(y = factor(nodes), x = moment))
    c + geom_raster(aes(fill = load), interpolate = FALSE) +
      #scale_fill_gradientn(colors = cols, breaks = leg_breaks) + 
      scale_fill_gradientn(colors = cols) + 
      scale_x_continuous(breaks = x_breaks) +
      labs(fill = "CPU load:") +
      xlab("time (seconds)") +
      ylab("cluster machine (node)") + 
      #ggtitle(paste("heatmap for", algo)) +
      theme_minimal(base_size = 16) +
      theme(legend.position="bottom",
            axis.ticks.length=unit(0.1,"cm"),
            axis.ticks = element_line(size = 0.25)) 
}

# plot the partitions found by each one of the clustering algorithms 
plot.heatmap.cluster.f <- function(dataset, algo, clalgo, sol, centroid, noise.value = 0) {
  
  target <- which(dataset$algo == algo)
  
  nodes <- as.integer(sapply(dataset$hostname[target], function(x) { substr(x, 15, 15) })) + 1
  moment <- round((dataset$timestamp[target] - min(dataset$timestamp[target])) + 1, 0)
  clusters <-  convert.real.f(as.integer(sol)[target], centroid)
  
  heat = data.frame("nodes" = nodes, "clusters"= as.character(clusters), "moment" = moment)
  
  #cols <- rev(heat_hcl(max(unique(clusters))))
  pal.c = c(90,40)
  pal.l=c(33,93)
  pal.power=c(1/5, 1.5)
  cols <- rev(heat_hcl(max(unique(clusters)), c = pal.c, l=pal.l, power=pal.power))
  lbls <- sort(unique(clusters))
  
  # if there is more than one cluster and the solution has noise
  if( (length(unique(clusters)) > 1) && has.noise.f(sol)) {
    
    # adds one (+1) to cols vector due to the max(unique(clusters)) having 0 (noise)
    cols <- rev(heat_hcl(max(unique(clusters)) +1, c = pal.c, l=pal.l, power=pal.power))
    noise.i = which(unique(clusters) == noise.value)
    cols[1] = "gray30"
    lbls[1] = "noise"
    
    # if there is only a single clusters
  } else if (length(unique(clusters)) == 1) {
    n.color = 100
    cols[1] = rev(heat_hcl(n.color), c = pal.c, l=pal.l, power=pal.power)[n.color/2] # close to yellow (pastel)
    
    # if there is no an exact correspondece for all partition numbers
    # (remember: k-means vs. gama case, for clusters 1,2,4,5 and 1,2,3,5)
  } else if (length(cols) > length(lbls)) {
    t <- 1:length(cols)
    #print(paste(algo, cols, lbls))
    cols <- cols [-c(outersect.f(t, lbls))]
  }
  
  
  c <- ggplot(heat, aes(y = factor(nodes), x = moment)) 
  c + geom_raster(aes(fill = clusters), interpolate = FALSE) +
      scale_fill_manual(values = cols, labels = lbls) +
      labs(fill = "cluster:") +
      xlab("time (seconds)") +
      ylab("cluster machine (node)") + 
      #ggtitle(paste("clustering of", algo, "by", clalgo)) +
    theme_minimal(base_size = 24) +
    theme(legend.position="bottom",
          axis.ticks.length=unit(0.1,"cm"),
          axis.ticks = element_line(size = 0.25)) 
}

plot.time.series.f <- function(dataset, algo) {
  
  target <- which(dataset$algo == algo)
  nodes <- as.integer(sapply(dataset$hostname[target], function(x) { substr(x, 15, 15) })) + 1
  load <- dataset$user[target] + dataset$system[target] + dataset$iowait[target] + dataset$softirq[target]
  moment <- round((dataset$timestamp[target] - min(dataset$timestamp[target])) + 1, 0)
  
  # choose length(unique(target) different values for colors in heat.colors
  heat = data.frame("nodes" = nodes, "load" = load, "moment" = moment)
  
  plot(x = heat[heat$nodes == 1,]$moment, y = heat[heat$nodes == 1,]$load, type = "l", ylim=c(0,100))
  for (i in 2:7)
    points(x = heat[heat$nodes == i,]$moment, y = heat[heat$nodes == i,]$load, type = "l")
  
}

# choose here the DML algorithm to plot heatmap and clusters
# als    bayes  gbt    kmeans lda    linear lr     pca    rf     svd    svm
dml = "bayes"

# uncomment these lines to generate graphics for all DML algorithms
# caution: may cause problems with R graphics devices, it is better to run each DML separately
#dmls <- as.character(unique(cpu.large.all$algo))
#for (dml in dmls) {
    cat("Plotting resource load heatmap for DML algorithm \n")
    pdf(file = paste("images/", dml, "-heatmap.pdf", sep=""), width=6,height=5,paper='special')
    plot.heatmap.load.f(cpu.large.all, dml)
    #plot.time.series.f(cpu.large.all, dml)
    garbage <- dev.off()

    # cat("Plotting partitions map (k-means)\n")
    # png(file = paste("images/", dml, "-cl-kmeans.png", sep=""))
    # plot.heatmap.cluster.f(cpu.large.all, dml, "k-means", km_sol[theBestOnes$kmeans[1],], km_centroids[theBestOnes$kmeans[1]])
    # garbage <- dev.off()
    # 
    # cat("Plotting partitions map (GAMA)\n")
    # png(file = paste("images/", dml, "-cl-gama.png", sep=""))
    # plot.heatmap.cluster.f(cpu.large.all, dml, "GAMA", gama_sol[theBestOnes$gama[1],], gama_centroids[theBestOnes$gama[1]])
    # garbage <- dev.off()
    # 
    # cat("Plotting partitions map (DBSCAN)\n")
    # png(file = paste("images/", dml, "-cl-dbscan.png", sep=""))
    # plot.heatmap.cluster.f(cpu.large.all, dml, "DBSCAN", dbscan_sol[theBestOnes$dbscan[1],], dbscan_centroids[theBestOnes$dbscan[1]])
    # garbage <- dev.off()
    # 
    # cat("Plotting partitions map (NKHGA)\n")
    # png(file = paste("images/", dml, "-cl-nkhga.png", sep=""))
    # plot.heatmap.cluster.f(cpu.large.all, dml, "NKHGA", nkhga_sol[theBestOnes$nkhga[1],], nkhga_centroids[theBestOnes$nkhga[1]])
    # garbage <- dev.off()
    # 
    # cat("\n\nWell done!\n\n")
    
#}
    
    

 