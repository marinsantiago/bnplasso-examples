# ------------------------------------------------------------------------------
# Data generating mechanism
# ------------------------------------------------------------------------------

library(mvtnorm) |> suppressPackageStartupMessages() |> suppressWarnings()

gen_data <- function(n, rho, p = 200, n.test = 1000, sigma2 = 1) {
  
  # Generate vector of regression coefficients ---------------------------------
  p.clust1 <- 0.025 * p
  p.clust2 <- 0.075 * p
  p.clust3 <- p - (p.clust1 + p.clust2)
  true.beta <- c(rep(10, p.clust1), rep(2, p.clust2), rep(0, p.clust3))
  
  # Generate matrix of predictors ----------------------------------------------
  Sigma.x <- rho^abs(matrix(1:p - 1, p, p, byrow = TRUE) - (1:p - 1))
  X <- mvtnorm::rmvnorm(n, rep(0, p), Sigma.x)
  X.test <- mvtnorm::rmvnorm(n.test, rep(0, p), Sigma.x)
  
  # Generate vector of responses -----------------------------------------------
  y <- rnorm(n, X %*% true.beta, sqrt(sigma2))
  y.test <- rnorm(n.test, X.test %*% true.beta, sqrt(sigma2))
  
  list(
    true.beta = true.beta, true.sigma2 = sigma2,
    X = X, X.test = X.test, y = y, y.test = y.test
  )
}
