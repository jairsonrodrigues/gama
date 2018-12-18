gama.plot.partitions <- function(gama.obj = NULL, view = "pca", ...) {

  dat <- gama.obj@original.data
  dat$clusters <- gama.obj@cluster

  if (view == "total.sum") {

    total.sum = apply(dat, 1, sum)

    dat$total.sum <- total.sum
    dat$observation <- 1:nrow(dat)
    g <- ggplot2::ggplot(dat, ggplot2::aes(x = observation, y = total.sum, color = factor(clusters))) +
      ggplot2::geom_point() +
      ggplot2::labs(color = "partition") +
      ggplot2::xlab("observation") +
      ggplot2::ylab("total sum of dimensions")

  } else if (view == "pca") {

    dims <- ncol(gama.obj@original.data)

    if (dims > 2) {
      pca = prcomp(dat)
      x1 <- pca$x[,"PC1"]
      x2 <- pca$x[,"PC2"]
      lbl.1 <- "pc 1"
      lbl.2 <- "pc 2"
    } else {
      x1 <- dat[,1]
      x1 <- dat[,2]
      lbl.1 <- "x1"
      lbl.2 <- "x2"
    }

    g <- ggplot2::ggplot(dat, ggplot2::aes(x = x1, y = x2, color = factor(clusters))) +
      ggplot2::geom_point() +
      ggplot2::labs(color = "partition") +
      ggplot2::xlab(lbl.1) +
      ggplot2::ylab(lbl.2)
  }

  g <- g + ggplot2::ggtitle(paste("Gama partitions,", "view = ", view, sep = " ")) + ggplot2::theme_minimal()
  plot(g)
}

# gama.plot.dataset <- function(data = NULL) {
#
#   dims <- ncol(data)
#
#   if (dims > 2) {
#     pca = prcomp(data)
#     x1 <- pca$x[,"PC1"]
#     x2 <- pca$x[,"PC2"]
#     lbl.1 <- "pc 1"
#     lbl.2 <- "pc 2"
#   } else {
#     x1 <- data[,1]
#     x1 <- data[,2]
#     lbl.1 <- "x1"
#     lbl.2 <- "x2"
#   }
#
#   g <- ggplot2::ggplot(data, ggplot2::aes(x = x1, y = x2)) +
#     ggplot2::geom_point() +
#     ggplot2::labs(color = "partition") +
#     ggplot2::xlab(lbl.1) +
#     ggplot2::ylab(lbl.2) + ggplot2::theme_minimal()
#
#   g
# }
