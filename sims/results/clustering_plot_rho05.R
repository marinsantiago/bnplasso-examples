# ------------------------------------------------------------------------------
# Clustering plot of the regression coefficients when rho = 0.5 and sigma2 = 1
# ------------------------------------------------------------------------------

library(bnplasso)
source("./R/gen_data.R"); source("./R/io_results.R")
source("./R/cls_plt.R"); source("./R/vi.R"); source("./R/vio_plt.R")

rho <- 0.5
reps <- 50L
n.vals <- c(100, 250, 500)
lambdas.zero.coef <- lapply(seq_along(n.vals), \(.) rep(NA, reps))
lambdas.nonzero.coef <- lapply(seq_along(n.vals), \(.) rep(NA, reps))
clusters <- lapply(seq_along(n.vals), \(.) matrix(nrow = reps, ncol = 200))
post.clusters.k <- lapply(seq_along(n.vals), \(.) rep(NA, reps))

for (i in seq_along(n.vals)) {
  n <- n.vals[i]
  for (r in seq_len(reps)) {
    set.seed(r)
    # Make sure to use the correct values of rho and n
    data <- gen_data(n = n, rho = rho)
    true.beta <- data$true.beta; true.sigma2 <- data$true.sigma2 # True params.
    X <- data$X; X.test <- data$X.test; y <- data$y; y.test <- data$y.test
    rm(data); gc()
    set.seed(1)
    out.bnp <- bnplasso::bnplasso.lm(X, y, a = 0.1, b = 0.01/n, alpha = 0.6)
    clusts <- vi(out.bnp$Post.clust_idx) # Point estimate of clustering
    clusts.idx <- unique(clusts)
    k.mean.post <- mean(out.bnp$Post.K) # Posterior mean number of clusters
    post.clusters.k[[i]][r] <- k.mean.post 
    lambda.sqrt <- sqrt(out.bnp$Post.lambda2)
    lambdas <- sapply(clusts.idx, \(c) mean(lambda.sqrt[,clusts == c]))
    lambdas.zero.coef[[i]][r] <- max(lambdas)
    lambdas.nonzero.coef[[i]][r] <- min(lambdas)
    # Label switching
    if (length(clusts.idx) == 2L) {
      if (clusts.idx[1] > clusts.idx[2]) {
        # Re-label
        clusts.new <- rep(NA, 200)
        clusts.new[clusts == 2] <- 1
        clusts.new[clusts == 1] <- 2
        clusts <- clusts.new
      }
    }
    clusters[[i]][r,] <- clusts
    rm(out.bnp, lambda.sqrt, clusts); gc()
    print(paste("n: ", n, "Replication: ", r))
  }
}

round(sapply(lambdas.nonzero.coef, mean), 2) # n = 100, 250, 500
round(sapply(lambdas.zero.coef, mean), 2) # n = 100, 250, 500

violin_plot_k(post.clusters.k, rho, 1)

p1 <- clusters_plot(clusters[[1]], n = n.vals[1])
p2 <- clusters_plot(clusters[[2]], n = n.vals[2])
p3 <- clusters_plot(clusters[[3]], n = n.vals[3])

# Arrange the plots in a 3x1 grid ----------------------------------------------
#png("./clustering_rho05.png", width = 6000, height = 7000, res = 900)
grid.arrange(p1, NULL, p2, NULL, p3, ncol = 1, heights = c(1, 0.1, 1, 0.1, 1))
#dev.off()
