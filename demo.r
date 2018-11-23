source("helpers.r")
data("als.cpu.gig")

gama.res <- gama(als.cpu.gig, plot.internals = T, generations = 100)
plot(gama.res, view.method = "pca") 

nkn
