# ------------------------------------------------------------------------------
# Co-clustering analysis for different values of n with rho = 0.3
# ------------------------------------------------------------------------------

library(bnplasso)
source("./R/co_clustering.R"); source("./R/gen_data.R")

# n = 100 ----------------------------------------------------------------------
set.seed(1)
p <- 200
rho <- 0.3
n <- 100
data <- gen_data(n = n, rho = rho)
true.beta <- data$true.beta; true.sigma2 <- data$true.sigma2 # True params.
X <- data$X; X.test <- data$X.test; y <- data$y; y.test <- data$y.test
out.bnp.100 <- bnplasso::bnplasso.lm(X, y, a = 0.1, b = 0.01/n, alpha = 0.6)
posterior.co_clustering.probs(out.bnp.100, "n = 100")
# Plot size: 10 x 12. Landscape

# n = 250 ----------------------------------------------------------------------
set.seed(1)
n <- 250
data <- gen_data(n = n, rho = rho)
true.beta <- data$true.beta; true.sigma2 <- data$true.sigma2 # True params.
X <- data$X; X.test <- data$X.test; y <- data$y; y.test <- data$y.test
out.bnp.250 <- bnplasso::bnplasso.lm(X, y, a = 0.1, b = 0.01/n, alpha = 0.6)
posterior.co_clustering.probs(out.bnp.250, "n = 250")
# Plot size: 10 x 12. Landscape

# n = 500 ----------------------------------------------------------------------
set.seed(1)
n <- 500
data <- gen_data(n = n, rho = rho)
true.beta <- data$true.beta; true.sigma2 <- data$true.sigma2 # True params.
X <- data$X; X.test <- data$X.test; y <- data$y; y.test <- data$y.test
out.bnp.500 <- bnplasso::bnplasso.lm(X, y, a = 0.1, b = 0.01/n, alpha = 0.6)
posterior.co_clustering.probs(out.bnp.500, "n = 500")
# Plot size: 10 x 12. Landscape
