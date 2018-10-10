source("helpers.r")

raw <- read.csv("data/cpu-all.csv")
raw <- raw[raw$exp == 'exp2',]
raw <- raw[raw$scale == 'bigdata',]
raw <- raw[raw$algo == 'kmeans',]
raw <- raw[,c("user", "system", "iowait", "softirq")]


gamaRes <- gama(raw)

plot.clusters(gamaRes)
