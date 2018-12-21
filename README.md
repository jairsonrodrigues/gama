# gama
We present GAMA, a Genetic Approach to MAximize a clustering criteria --- a R package for evolutionary hard partitional clustering, by using guided operators (those guided by some kind of information about the quality of individual clusters), for a fixed a priori known number of partitions, encoded as real-valued centroid based partitions.


To use, just call:

data(cpu.als)

# minimal call to gama, the best value for k will be estimated internally
obj <- gama(dataset = cpu.als)

# call gama to segment data into three partitions and 200 generations
obj <- gama(dataset = cpu.als, k = 3, generations = 200)

See documentation for more details.
