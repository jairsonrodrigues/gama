#R CMD BATCH '--args large 1000 25' gamas.r gamas.out

library(GA)
library(cluster)
library(Rfast)

# fitness_asw_nopenalty <- function(palpite) {
#   m.palpite <- matrix(palpite,nrow = k,ncol = dimen)
#   which.dists <- apply(dista(cpu.norm, m.palpite, 
#                              "euclidean", square = TRUE), 1, which.min)
#   
#   asw   <- silhouette(which.dists, d2)
#   
#   # tenta o cálculo da silhueta 
#   # retorna neutro (0) em caso de erro
#   # observação: não descobri qual o erro
#   sm <- tryCatch(summary(asw)$avg.width, error = function (e) { return (0)})
#   
#   return (sm)
# }

fitness_asw <- function(palpite) {
  
  m.palpite <- matrix(palpite,nrow = k,ncol = dimen)
  
  sums <- apply(m.palpite, 1, sum)
  overflow <- which(sums > 100)
  num_constraints = length(overflow)
  
  penalty = 0
  if (num_constraints > 0) { 
    penalty <- num_constraints * max(abs(sums[overflow] -100)/sums[overflow]) 
  }
  
  which.dists <- apply(dista(cpu.norm, m.palpite, "euclidean", square = TRUE), 1, which.min)

  asw   <- silhouette(which.dists, d2)

  # tenta o cálculo da silhueta 
  # retorna neutro (0) em caso de erro
  sm <- tryCatch(summary(asw)$avg.width, error = function (e) { return (0)})
  
  return (sm - penalty)
}


pop.f <- function(object) {
  
  lower <- object@lower
  upper <- object@upper
  nvars <- length(lower)
  
  population <- matrix(as.double(NA), nrow = object@popSize, ncol = nvars)
  for(j in 1:nvars) 
  { population[,j] <- runif(object@popSize, lower[j], upper[j]) }
  
  # valida cada conjunto de pontos (population [i,1:4])
  for(i in 1:object@popSize) { 
    
    # diminui o valor de cada linha até que o somatório seja <= 100
    repeat {
      
      # quais linhas (1:k) do individuo na populacao (1:popSize) sao maiores que 100
      m <- matrix(population[i,], ncol = 4, nrow = 5)
      
      sums <- apply(m, 1, sum)
      w.sums <- which(sums > 100)
      
      # para o laco repeat quando a quantidade de linhas com somatorio > 100 for igual a 0
      if (length(w.sums) == 0){
        break
      } else {
        # diminui o valor de cada linha em uma proporção equivalente a cada dimensão sobre o total 
        for (y in w.sums) {
          
          # calcula a proporção sobre a o indivíduo desnormalizado
          prop_over_all <- m[y,]/sums[y]
          m[y, ] <- m[y, ] - prop_over_all
          population[i,] <- as.vector(m)
        }
      }
    }
    
  }
  
  return(population)
}

# ==== PARÂMETROS ==== 
# R CMD BATCH '--args scale = large'
args=(commandArgs(TRUE))
if (length(args)!=3) {
  stop("Give the data volume scale, the number of generations and iterations.n", call.=FALSE)
} 

#eval(parse(text=args[[1]]))
scale = args[[1]]
generations = as.integer(args[[2]])
iterations = as.integer(args[[3]])
print(args)

#for (scale in scales) {

      k = 5 # cluster
      #cpu <- read.csv(paste("data/cpu-", scale, ".csv", sep=""))
      cpu.norm <- read.csv(paste("../data/cpu-", scale, ".csv", sep=""))
      
      # dimensions
      dimen = ncol(cpu.norm)
      
      # distance matrix
      d <- dist(cpu.norm, method = "euclidean", diag = FALSE, upper = FALSE)
      d2 <- d^2
      rm(d)
      gc()
      
      # lower and upper bounds for the operators
      lower_bound <-  c(rep(min(cpu.norm[, 1]), k), rep(min(cpu.norm[, 2]), k), rep(min(cpu.norm[, 3]), k), rep(min(cpu.norm[, 4]), k))
      upper_bound <-  c(rep(max(cpu.norm[, 1]), k), rep(max(cpu.norm[, 2]), k), rep(max(cpu.norm[, 3]), k), rep(max(cpu.norm[, 4]), k))
       
      s <- "gareal_lsSelection"
      m <- "gareal_rsMutation"
      c <- "gareal_blxCrossover"
      
      sl = "linScaling"
      cl = "blend"
      ml = "randAroundSol"
      
      
      pop.r <- 50
      cross.r <- 0.9
      mut.r <- 0.01
      #generations <- 10
      
      # elitismo = 5% da populacao
      elit.r = floor(pop.r * 0.05)
      
      # numero de iteracoes
      #iterations = 1
      
      asw <- c()
      solution <- list()
      
      for (w in 1:iterations) {
        
              print(paste(scale, "Iter: ", w, sep = (" ")))
              
              # para que a seed seja diferente dos primeiros experimentos
              seed = w*100
              start.time <- Sys.time()
              
              genetic <- ga(type = "real-valued", 
                            seed = w, 
                            population = pop.f,
                            selection = s, mutation = m, crossover = c,
                            popSize = pop.r, elitism = elit.r, pmutation = mut.r, pcrossover = cross.r, 
                            maxiter = generations, maxFitness = 1.0, fitness = fitness_asw,
                            lower = lower_bound,
                            upper = upper_bound,
                            monitor = FALSE)
              
              end.time <- Sys.time()
              
              # ==== COMPUTA O ASW =====
              # se o GA retornou apenas uma solução
              individuo = NA
              num_sols = length(genetic@solution)/(k*dimen)
              if (num_sols == 1) { 
                individuo <- matrix(genetic@solution,nrow = k,ncol = dimen)
              } else {
                # se o GA retornou mais que uma solução (toma a primeira, igual fitness)
                individuo <- matrix(genetic@solution[1,],nrow = k,ncol = dimen)
              }
              
              which.dists <- apply(dista(cpu.norm, individuo, "euclidean", square = TRUE), 1, which.min)
              asw[w] <- summary(silhouette(which.dists, d2))$avg.width
              solution[[w]] <- which.dists
              
              t <- as.data.frame(individuo)
              colnames(t) <- colnames(cpu.norm)
              t <- t[with(t, order(user + system + softirq + iowait)), ]
              
              sink(file=paste(scale,"/points_data.dat", sep = ""), append = TRUE)
              print(t)
              cat("\n")
              sink()
              
              sink(file=paste(scale,"/sol_data.dat", sep=""), append = TRUE)
              cat(solution[[w]])
              cat("\n")
              sink()
              
              sink(file=paste(scale,"/summary_data.dat", sep=""), append = TRUE)
              cat(asw[[w]])
              cat("\n")
              sink()
              
              rm(genetic, which.dists)
        
      }
      
      sink(file=paste(scale,"/summary_data.dat", sep=""), append = TRUE)
      print(summary(asw))
      cat("\n")
      sink()
      
      rm(d2)
      gc()
#}