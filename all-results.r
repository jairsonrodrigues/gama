library(cluster)
library(dplyr)

load(file = "RData/metrics.RData")
load(file = "RData/solutions.RData")

# distance matrix
d <- dist(cpu.large, method = "euclidean", diag = FALSE, upper = FALSE)
d2 <- d^2

# compute ASW
asw.kmeans <- c()
asw.dbscan <- c()
asw.nkhga <- c()
asw.gama <- c()

for (i in 1:25) {
  asw.kmeans[i] <- summary(silhouette(km_sol[i,], d2))$avg.width
  asw.dbscan[i] <- summary(silhouette(dbscan_sol[i,], d2))$avg.width
  asw.nkhga[i] <- summary(silhouette(nkhga_sol[i,], d2))$avg.width
  asw.gama[i] <- summary(silhouette(gama_sol[i,], d2))$avg.width
}

# load Matlab DBCV
dbcv.kmeans <- read.csv("dbcv/dbcv-kmeans.dat", header = FALSE, col.names = FALSE)
dbcv.nkhga <- read.csv("dbcv/dbcv-nkhga.dat", header = FALSE, col.names = FALSE)
dbcv.dbscan <- read.csv("dbcv/dbcv-dbscan.dat", header = FALSE, col.names = FALSE)
dbcv.gama <- read.csv("dbcv/dbcv-gama.dat", header = FALSE, col.names = FALSE)

eval.asw <- data.frame( 
  group = rep(c("kmeans", "nkhga", "dbscan", "gama"), each = 25),
  v = c(t(asw.kmeans),  t(asw.nkhga), t(asw.dbscan), t(asw.gama))
)

eval.dbcv <- data.frame( 
  group = rep(c("kmeans", "nkhga", "dbscan", "gama"), each = 25),
  v = c(t(dbcv.kmeans),  t(dbcv.nkhga), t(dbcv.dbscan), t(dbcv.gama))
)

print(group_by(eval.asw, group) %>%
  summarise(
    #count = n(),
    min = min(v, na.rm = TRUE),
    max = max(v, na.rm = TRUE),
    median = median(v, na.rm = TRUE),
    mean = mean(v, na.rm = TRUE),
    sd = sd(v, na.rm = TRUE)))

print(group_by(eval.dbcv, group) %>%
  summarise(
    #count = n(),
    min = min(v, na.rm = TRUE),
    max = max(v, na.rm = TRUE),
    median = median(v, na.rm = TRUE),
    mean = mean(v, na.rm = TRUE),
    sd = sd(v, na.rm = TRUE)))

wt_gama_kmeans_asw <- wilcox.test(asw.gama, asw.kmeans, paired = TRUE, exact = FALSE)
wt_gama_nkhga_asw <- wilcox.test(asw.gama, asw.nkhga, paired = TRUE, exact = FALSE)
wt_gama_dbscan_asw <- wilcox.test(asw.gama, asw.dbscan, paired = TRUE, exact = FALSE)

wt_gama_kmeans_dbcv <- wilcox.test(t(dbcv.gama), t(dbcv.kmeans), paired = TRUE, exact = FALSE)
wt_gama_nkhga_dbcv <- wilcox.test(t(dbcv.gama), t(dbcv.nkhga), paired = TRUE, exact = FALSE)
wt_gama_dbscan_dbcv <- wilcox.test(t(dbcv.gama), t(dbcv.dbscan), paired = TRUE, exact = FALSE)

print(wt_gama_kmeans_asw)
print(wt_gama_nkhga_asw)
print(wt_gama_dbscan_asw)
print(wt_gama_kmeans_dbcv)
print(wt_gama_nkhga_dbcv)
print(wt_gama_dbscan_dbcv)

algos <- c("kmeans", "gama", "dbscan", "nkhga")
theBestOnes <- as.data.frame(sapply(algos, function (x) {
      c(which.max(eval.asw[eval.asw$group == x,]$v), 
        max(eval.asw[eval.asw$group == x,]$v),
        which.max(eval.dbcv[eval.dbcv$group == x,]$v),
        max(eval.dbcv[eval.dbcv$group == x,]$v))
  }
))

save(eval.asw, eval.dbcv, theBestOnes, file = "RData/criteria.RData")
