# ------------------------------------------------------------------------------
# Simulation setting with n = 500, rho = 0.3 and sigma2 = 5, BB-SSL (RUN 1/10)
# ------------------------------------------------------------------------------

setwd("~/Desktop/bnplasso-examples")

library(bnplasso)
source("./R/metrics.R"); source("./R/bbssl.R")
source("./R/gen_data.R"); source("./R/io_results.R")

# Simulation settings ----------------------------------------------------------
n <- 500
rho <- 0.3
sigma2 <- 5
reps <- 200

# Prepare the returns ----------------------------------------------------------
MSE.bbssl <- rep(NA, reps)
ACC.bbssl <- MSPE.bbssl <- TIME.bbssl <- ELPPD.bbssl <- MSE.bbssl

print("Starting batch 1")

# Run the simulations ----------------------------------------------------------
for (r in 1:20) {
  
  # Generate the data ----------------------------------------------------------
  set.seed(r)
  data <- gen_data(n = n, rho = rho, sigma2 = sigma2)
  true.beta <- data$true.beta; true.sigma2 <- data$true.sigma2 # True parameters
  X <- data$X; X.test <- data$X.test; y <- data$y; y.test <- data$y.test
  rm(data); gc()
  
  # Fit the models -------------------------------------------------------------
  set.seed(1)
  # BBSSL
  out.bbssl <- bbssl(X, y)
  
  # Point estimates ------------------------------------------------------------
  # beta.hat
  beta.hat.bbssl <- bnplasso::point.estimates(out.bbssl, retain = "mean")
  # y.hat
  y.hat.bbssl <- X.test %*% beta.hat.bbssl
  
  # Assessment metrics ---------------------------------------------------------
  MSE.bbssl[r] <- mse(beta.hat.bbssl, true.beta)
  ACC.bbssl[r] <- sel.accuracy(beta.hat.bbssl, true.beta)
  MSPE.bbssl[r] <- mspe(y.hat.bbssl, y.test)
  ELPPD.bbssl[r] <- bnplasso::elppd(out.bbssl, X.test, y.test)
  TIME.bbssl[r] <- out.bbssl$elapsed
  
  rm(out.bbssl); gc()
  print(paste("Replication: ", r))
}

gc()

print("Ending batch 1")

# Export the results -----------------------------------------------------------

export.results.bbssl(n, rho, sigma2 = sigma2, reps, run = 1)
