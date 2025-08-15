# ------------------------------------------------------------------------------
# Co-clustering analysis for different tails in G0
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

# (a = 0.01, b = 0.01/n) -------------------------------------------------------
set.seed(1)
out.bnp.0.1 <- bnplasso::bnplasso.lm(X, y, a = 0.01, b = 0.01/n, alpha = 0.6)
posterior.co_clustering.probs(out.bnp.0.1, "(a = 0.01, b = 0.01/n)")
# Plot size: 10 x 12. Landscape

# (a = 0.1, b = 0.1/n) ---------------------------------------------------------
set.seed(1)
out.bnp.1 <- bnplasso::bnplasso.lm(X, y, a = 0.1, b = 0.1/n, alpha = 0.6)
posterior.co_clustering.probs(out.bnp.1, "(a = 0.1, b = 0.1/n)")
# Plot size: 10 x 12. Landscape

# (a = 1, b = 0.1/n) -----------------------------------------------------------
set.seed(1)
out.bnp.10 <- bnplasso::bnplasso.lm(X, y, a = 1, b = 0.1/n, alpha = 0.6)
posterior.co_clustering.probs(out.bnp.10, "(a = 1, b = 0.1/n)")
# Plot size: 10 x 12. Landscape
