# ------------------------------------------------------------------------------
# Co-clustering analysis for different values of alpha = 0.1, 1.0, 10.0
# ------------------------------------------------------------------------------

library(bnplasso)
source("./R/co_clustering.R"); source("./R/gen_data.R")

# Simulation settings ----------------------------------------------------------
set.seed(1)
n <- 250
p <- 200
rho <- 0.3
data <- gen_data(n = n, rho = rho)
true.beta <- data$true.beta; true.sigma2 <- data$true.sigma2 # True params.
X <- data$X; X.test <- data$X.test; y <- data$y; y.test <- data$y.test

# alpha = 0.1 ------------------------------------------------------------------
set.seed(1)
out.bnp.0.1 <- bnplasso::bnplasso.lm(X, y, a = 0.1, b = 0.01/n, alpha = 0.1)
posterior.co_clustering.probs(out.bnp.0.1, bquote(alpha == .("0.1")))
# Plot size: 10 x 12. Landscape

# alpha = 1.0 ------------------------------------------------------------------
set.seed(1)
out.bnp.1 <- bnplasso::bnplasso.lm(X, y, a = 0.1, b = 0.01/n, alpha = 1.0)
posterior.co_clustering.probs(out.bnp.1, bquote(alpha == .("1.0")))
# Plot size: 10 x 12. Landscape

# alpha = 10.0 -----------------------------------------------------------------
set.seed(1)
out.bnp.10 <- bnplasso::bnplasso.lm(X, y, a = 0.1, b = 0.01/n, alpha = 10.0)
posterior.co_clustering.probs(out.bnp.10, bquote(alpha == .("10.0")))
# Plot size: 10 x 12. Landscape
