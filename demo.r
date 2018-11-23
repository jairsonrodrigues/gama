source("helpers.r")
data("pca.cpu.big.2")

gama.res <- gama(pca.cpu.big.2, plot.internals = T, generations = 250, seed.p = 17)
plot(gama.res, view.method = "pca")

