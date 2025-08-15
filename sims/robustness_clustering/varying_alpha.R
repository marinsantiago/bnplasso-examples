# ------------------------------------------------------------------------------
# Clustering results for different values of alpha = 0.1, 1.0, 10.0
# ------------------------------------------------------------------------------

library(bnplasso)
source("./R/metrics.R"); source("./R/cls_plt.R"); source("./R/io_results.R")
source("./R/gen_data.R"); source("./R/vio_plt.R"); source("./R/vi.R")

# Simulation settings ----------------------------------------------------------
n <- 250
p <- 200
rho <- 0.3
reps <- 50L
alpha.vals <- c(0.1, 1.0, 10.0)

# Prepare the returns ----------------------------------------------------------
CLUSTS <- lapply(seq_along(alpha.vals), \(.) matrix(NA, nrow = reps, ncol = p))
LAMBDAS.ZERO <- LAMBDAS.NONZERO <- lapply(CLUSTS, \(.) rep(NA, reps))
post.clusters.k <- lapply(seq_along(alpha.vals), \(.) rep(NA, reps))  

# Run the simulations ----------------------------------------------------------
for (.al in seq_along(alpha.vals)) {
  alpha <- alpha.vals[.al]
  for (r in seq_len(reps)) {
    set.seed(r)
    data <- gen_data(n = n, rho = rho)
    true.beta <- data$true.beta; true.sigma2 <- data$true.sigma2 # True params.
    X <- data$X; X.test <- data$X.test; y <- data$y; y.test <- data$y.test
    rm(data); gc()
    set.seed(1)
    # Make sure to use the current alpha value.
    out.bnp <- bnplasso::bnplasso.lm(X, y, a = 0.1, b = 0.01/n, alpha = alpha)
    # Clustering of the regression coefficients
    clusts <- vi(out.bnp$Post.clust_idx) # Point estimate of clustering
    clusts.idx <- unique(clusts)
    k.mean.post <- mean(out.bnp$Post.K) # Posterior mean number of clusters
    post.clusters.k[[.al]][r] <- k.mean.post 
    lambda.sqrt <- sqrt(out.bnp$Post.lambda2)
    lambdas <- sapply(clusts.idx, \(c) mean(lambda.sqrt[,clusts == c]))
    LAMBDAS.ZERO[[.al]][r] <- max(lambdas) 
    LAMBDAS.NONZERO [[.al]][r] <- min(lambdas)
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
    CLUSTS[[.al]][r,] <- clusts
    rm(out.bnp, lambda.sqrt); gc()
    print(paste("alpha: ", alpha, "Replication: ", r))
  }
}

# Export the results -----------------------------------------------------------
path <- "./sims/robustness_clustering/"
file.type <- "_varying_alpha.Rdata"
save(CLUSTS, file = paste0(path, "CLUSTS", file.type))
save(LAMBDAS.ZERO, file = paste0(path, "LAMBDAS.ZERO", file.type))
save(LAMBDAS.NONZERO, file = paste0(path, "LAMBDAS.NONZERO", file.type))
save(post.clusters.k, file = paste0(path, "post_clusters", file.type))

# Load the results -------------------------------------------------------------
load(paste0(path, "CLUSTS", file.type))
load(paste0(path, "LAMBDAS.ZERO", file.type))
load(paste0(path, "LAMBDAS.NONZERO", file.type))
load(paste0(path, "post_clusters", file.type))

# Individual plots -------------------------------------------------------------
p1 <- clusters_plot_alpha(CLUSTS[[1]], alpha = alpha.vals[1])
p2 <- clusters_plot_alpha(CLUSTS[[2]], alpha = alpha.vals[2])
p3 <- clusters_plot_alpha(CLUSTS[[3]], alpha = alpha.vals[3])

# Arrange the plots in a 3x1 grid ----------------------------------------------
#png("./clustering_alphas.png", width = 6000, height = 7000, res = 900)
grid.arrange(p1, NULL, p2, NULL, p3, ncol = 1, heights = c(1, 0.1, 1, 0.1, 1))
#dev.off()

# Lambda values ----------------------------------------------------------------
sapply(LAMBDAS.ZERO, mean) # alpha = 0.1, 1.0, 10.0
sapply(LAMBDAS.NONZERO, mean) # alpha = 0.1, 1.0, 10.0

# Violin plots -----------------------------------------------------------------
violin_plot_k_varying_alpha(
  post.clusters.k, "0.1", "0.01/n", n = 250, rho_ = 0.3,
)
