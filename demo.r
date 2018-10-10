source("helpers.r")


raw <- read.csv("data/cpu-all.csv")
raw <- raw[raw$exp == 'exp2',]
raw <- raw[raw$scale == 'bigdata',]
raw <- raw[raw$algo == 'als',]
raw <- raw[,c("user", "system", "iowait", "softirq")]

k = 3
obj <- gama(data = raw, k, generations = 100)

# === OBS melhorar o eixo X ====
d <-data
d$clusters <- obj@cluster
d$moment <- 1:nrow(d)

g <- ggplot(d, aes(y = user + system + iowait + softirq, x = moment, color = factor(clusters))) + 
  geom_point() + theme_minimal()
g