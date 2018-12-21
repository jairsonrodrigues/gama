# This function takes a gama object and plots the partitions found by the algorithm.
# The partitions will be coloured accordingly the distribution of the clusters.
gama.plot.partitions <- function(gama.obj = NULL, view = "pca", ...) {

  # --- arguments validation --- #
  Check <- ArgumentCheck::newArgCheck()

  if (is.null(gama.obj))
    ArgumentCheck::addError(
      msg = "'gama.obj' can not be NULL",
      argcheck = Check
    )

  if (class(gama.obj) != 'gama')
    ArgumentCheck::addError(
      msg = "'gama.obj' must be a gama object.",
      argcheck = Check
    )

  if (is.null(view)) {
    ArgumentCheck::addError(
      msg = "'view' can not be NULL",
      argcheck = Check
    )
  } else if (!(view %in% c('pca', 'total.sum')))
    ArgumentCheck::addError(
      msg = "'view' must be etiher 'pca' or 'total.sum'.",
      argcheck = Check
    )

  ArgumentCheck::finishArgCheck(Check)

  # --- final of arguments validation --- #

  # to avoid R CMD check problems (the checker considers
  # these variables as global not declared)
  observation <- NA
  total.sum <- NA
  clusters <- NA

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
      # remove the last column (cluster ids) from PCA calculus
      pca = prcomp(dat[-(dims +1)])
      x1 <- pca$x[,"PC1"]
      x2 <- pca$x[,"PC2"]
      lbl.1 <- "pc 1"
      lbl.2 <- "pc 2"
    } else {
      x1 <- dat[,1]
      x2 <- dat[,2]
      lbl.1 <- "x1"
      lbl.2 <- "x2"
    }

    df <- data.frame(x1, x2, dat$clusters)
    colnames(df) <- c("x1", "x2", "clusters")

    y.min <- min(x2)
    y.max <- max(x2)

    g <-  ggplot2::ggplot(df, ggplot2::aes(x = x1, y = x2, color = factor(clusters))) +
          ggplot2::geom_point(size = 1.8) +
          ggplot2::labs(color = "partition") +
          ggplot2::xlab(lbl.1) +
          ggplot2::ylab(lbl.2) +
          ggplot2::ylim(y.min, y.max)
  }

  g <- g + ggplot2::ggtitle(paste("gama partitions", sep = "")) +
       ggplot2::theme_minimal(base_size = 16)
  plot(g)
}
