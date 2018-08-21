# [1] "cluster 1"
# cluster      user     system     iowait    softirq
# 5       5  3.786517  0.4483773 0.03808826 0.01324445
# 2       2 20.544683  2.2433864 0.06552283 0.12082967
# 1       1 34.502336  2.7702786 0.04011255 0.26831005
# 4       4 57.218328  3.6768818 0.01621823 0.27408665
# 3       3 27.309430 38.8017971 0.02792401 0.04602310
parse.km.f <- function(raw) {
  
  centroids = list()
  for (i in 0:24) {
    tmp <- raw[(2+(7*i)):(6+(7*i)),1]
    m_matrix = matrix(NA, nrow= 5, ncol = 5)
    for (j in 1:5){
      m_matrix[j,] <- unlist(lapply(strsplit(as.vector(tmp[j]), " ", fixed = TRUE), function(x){x[!x ==""]}))[2:6]
    }
    colnames(m_matrix) <- c("cluster.label", "user", "system", "iowait", "softirq")
    df <- as.data.frame(m_matrix)
    df$cluster.order <- 1:5
    centroids[[i+1]] <- apply(df, 2, as.numeric)
  }
  
  return (centroids)
}

# user    system    iowait   softirq
# 5  4.995023  3.103209 1.9617822 0.5702648
# 4 27.671183  6.941746 0.5749052 0.5985109
# 3 20.164766 32.996861 1.4910562 0.5615550
# 1 59.700555 14.559688 0.3301937 0.8035897
# 2 43.502298 46.378935 1.7966425 0.6967906

parse.gama.f <- function(raw) {
  
  centroids = list()
  for (i in 0:24) {
    tmp <- raw[(1+(6*i)):(5+(6*i)),1]
    m_matrix = matrix(NA, nrow= 5, ncol = 5)
    for (j in 1:5){
      m_matrix[j,] <- unlist(lapply(strsplit(as.vector(tmp[j]), " ", fixed = TRUE), function(x){x[!x ==""]}))[1:5]
    }
    colnames(m_matrix) <- c("cluster.label", "user", "system", "iowait", "softirq")
    df <- as.data.frame(m_matrix)
    df$cluster.order <- 1:5
    centroids[[i+1]] <- apply(df, 2, as.numeric)
  }
  
  return (centroids)
}

# infer the load level for DBSCAN and NKHGA (not based on centroids)
# OBS: dbscan lables the clusters in format 0..n, no salts
parse.dbscan.f <- function(raw) {
  
  centroids <- list()
  
  for (i in 1:nrow(raw)) {
    
    sols <- raw[i,]
    uniques <- unique(as.integer(sols))
    noise <- which(uniques == 0)
    
    if (length(noise) > 0)
      uniques <- uniques[-noise]
    
    # recover indices for each element by cluster
    w.indices <- sapply(uniques, function (x) { which(sols ==x) })
    
    mean.values <- c()
    cluster.labels <-c()
    for (el in w.indices) {
      
          mean.values <- c(mean.values, sum(cpu.large[el,])/nrow(cpu.large[el,]))
          cluster.labels <- c(cluster.labels, unique(as.integer(sols[el])))
    }
    
    # if (i == 5)
    #   print(i)
    
    tmp <- data.frame(cluster.labels, mean.values)
    tmp <- tmp[order(mean.values),]
    tmp$cluster.order <- 1:length(uniques)
    
    centroids[[i]] <- tmp
  }
  
  return (centroids)
}

# load metrics
load(file = "RData/metrics.RData")

# load cluster solutions
km_sol <- read.csv("solutions/kmeans/sol_data.dat", header = FALSE, sep = " ") 
gama_sol <- read.csv("solutions/gama/sol_data.dat", header = FALSE, sep = " ")
nkhga_sol <- read.csv("solutions/nkhga/sol_data.dat", header = FALSE, sep = " ")
dbscan_sol <- read.delim("solutions/dbscan/sol_data.dat", header = FALSE, sep = " ")
nkhga_sol <- nkhga_sol[- (length(nkhga_sol))] # remove a última coluna (em branco)

# load centroids for k-means and GAMA from solutions file
km_centroids_raw <- read.csv("solutions/kmeans/points_data.dat")
gama_centroids_raw <- read.csv("solutions/gama/points_data.dat")

# parse centroids for k-means and GAMA
km_centroids <- parse.km.f(km_centroids_raw)
gama_centroids <- parse.gama.f(gama_centroids_raw)

# infer clusters load levels from solutions (there is no centroids)
dbscan_centroids <- parse.dbscan.f(dbscan_sol)
nkhga_centroids <- parse.dbscan.f(nkhga_sol)

save(km_centroids, gama_centroids, dbscan_centroids, nkhga_centroids, 
     km_sol, gama_sol, dbscan_sol, nkhga_sol,
     file = "RData/solutions.RData")