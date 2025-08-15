# ------------------------------------------------------------------------------
# Asymptotic behavior of lambdas
# ------------------------------------------------------------------------------

library(bnplasso)
source("./R/metrics.R"); source("./R/gen_data.R"); source("./R/vi.R")

rho <- 0.3
reps <- 80L
n.vals <- seq(100, 2000, 100)
p <- 200

# Prepare the returns ----------------------------------------------------------
lambdas.zero.coef <- lapply(seq_along(n.vals), \(.) rep(NA, reps))
lambdas.nonzero.coef <- lapply(seq_along(n.vals), \(.) rep(NA, reps))

# Run the simulations ----------------------------------------------------------
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
    lambda.sqrt <- sqrt(out.bnp$Post.lambda2)
    lambdas <- sapply(clusts.idx, \(c) mean(lambda.sqrt[,clusts == c]))
    if (length(clusts.idx) == 2L) {
      lambdas.zero.coef[[i]][r] <- max(lambdas)
      lambdas.nonzero.coef[[i]][r] <- min(lambdas)
    }
    rm(out.bnp, lambda.sqrt); gc()
    print(paste("n: ", n, "Replication: ", r))
  }
}

# Store the results ------------------------------------------------------------
save(lambdas.zero.coef, file = "./sims/asymptotic_lambdas/lambda1.Rdata")
save(lambdas.nonzero.coef, file = "./sims/asymptotic_lambdas/lambda2.Rdata")

# Load the results -------------------------------------------------------------
load("./sims/asymptotic_lambdas/lambda1.Rdata")
load("./sims/asymptotic_lambdas/lambda2.Rdata")

# Plots ------------------------------------------------------------------------
lambda1.med <- sapply(lambdas.zero.coef, median)
lambda1.upr <- sapply(lambdas.zero.coef, \(x) quantile(x, 0.75))
lambda1.lwr <- sapply(lambdas.zero.coef, \(x) quantile(x, 0.25))
lambda1.exp <- expression(lambda[1])

lambda2.med <- sapply(lambdas.nonzero.coef, median)
lambda2.upr <- sapply(lambdas.nonzero.coef, \(x) quantile(x, 0.75))
lambda2.lwr <- sapply(lambdas.nonzero.coef, \(x) quantile(x, 0.25))
lambda2.exp <- expression(lambda[2])

# Lambda1
plot(
  x = n.vals, y = lambda1.upr, type = "n", xlab = "Sample size",
  ylab = lambda1.exp, main = lambda1.exp, ylim = c(50, 315), cex = 1.2,
  cex.main = 1.55, cex.lab = 1.2, cex.axis = 1.2, font.main = 2
)
polygon(
  x = c(n.vals, rev(n.vals)), y = c(lambda1.lwr, rev(lambda1.upr)),
  col = "#ADD8E6", border = NA
)
lines(x = n.vals, lambda1.med, type = "b", pch = 16, lwd = 2, col = "navy")


# Lambda 2
plot(
  x = n.vals, y = lambda2.lwr, type = "n", xlab = "Sample size",
  ylab = lambda2.exp, main = lambda2.exp, ylim = c(0.278, 0.329), cex = 1.2,
  cex.main = 1.55, cex.lab = 1.2, cex.axis = 1.2, font.main = 2
)
polygon(
  x = c(n.vals, rev(n.vals)), y = c(lambda2.lwr, rev(lambda2.upr)),
  col = "#FF00004D", border = NA
)
lines(x = n.vals, lambda2.med, type = "b", pch = 16, lwd = 2, col = "darkred")

