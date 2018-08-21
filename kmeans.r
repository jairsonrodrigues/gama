library(cluster)    

    file.name <- "data/cpu-large.csv"
    sol.cluster = "solutions/kmeans/sol_data.dat"
    sol.points  = "solutions/kmeans/points_data.dat"
    sol.summary = "solutions/kmeans/summary_data.dat"
    
    cpu <- read.csv(file.name, header = TRUE)

    # ==== PARÂMETROS ==== 
    k = 5 # cluster
    dimen = ncol(cpu)
    
    d <- dist(cpu, method = "euclidean", diag = FALSE, upper = FALSE)
    d2 <- d^2
    
    silhouettes <- vector()
    
    for (i in 1:25) {
      
      set.seed(i)  
      fit.km <- kmeans(cpu, algorithm = "Hartigan-Wong", nstart = 1, centers = k, iter.max = 100)
      centers <- aggregate(cpu, by=list(cluster=fit.km$cluster), mean)
      cluster <- fit.km$cluster
      
      centers <- centers[with(centers, order(user + system + softirq + iowait)), ]
      
      asw   <- silhouette(cluster, d2)  
      
      silhouettes[i] <- summary(asw)$avg.width
      
      sink(file=sol.cluster, append = TRUE, split=FALSE)
      cat(as.vector(cluster))
      cat("\n")
      sink()
      
      sink(file=sol.points, append = TRUE, split=FALSE)
      print(paste("cluster", i), sep = " ")
      print(centers)
      cat("\n")
      sink()
    }
    
    sink(file=sol.summary, append = TRUE)
    cat(silhouettes)
    cat("\n")
    cat(sd(silhouettes))
    cat("\n")
    print(summary(silhouettes))
    sink()