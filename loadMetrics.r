#
# Source code to load data metrics of experiments and save in memory variables.
# Author: [Jairson Rodrigues](jairsonrodrigues@gmail.com)
#  

# load data
cat("\n\nLoading data metrics ...")
cpu.large <- read.csv("data/cpu-large.csv")
cpu.large.all <- read.csv("data/all_exp2_large.csv")

cat("Saving in memory variables.")
save(cpu.large, cpu.large.all, file = "RData/metrics.RData")

cat("OK!\n\n")