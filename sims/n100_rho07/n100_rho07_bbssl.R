# ------------------------------------------------------------------------------
# Simulation setting with n = 100 and rho = 0.7, on BB-SSL
# ------------------------------------------------------------------------------

library(bnplasso)
source("./R/metrics.R"); source("./R/bbssl.R")
source("./R/gen_data.R"); source("./R/io_results.R")

# Simulation settings ----------------------------------------------------------
n <- 100
rho <- 0.7
reps <- 200
run <- 1

# Prepare the returns ----------------------------------------------------------
MSE.bbssl <- rep(NA, reps)
ACC.bbssl <- MSPE.bbssl <- TIME.bbssl <- ELPPD.bbssl <- MSE.bbssl

# Run the simulations ----------------------------------------------------------
for (r in seq_len(reps)) {
  
  # Generate the data ----------------------------------------------------------
  set.seed(r)
  data <- gen_data(n = n, rho = rho)
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

# Export the results -----------------------------------------------------------

export.results.bbssl(n, rho, sigma2 = 1, reps, run)

# Load the results -------------------------------------------------------------

path <- paste0(subfolder.path(n, rho, sigma2 = 1, reps, bbssl = T), "/")
paste0(path, current.set(n, rho, sigma2 = 1, reps), "_run", run, ".Rdata") |> load()

# Print the results ------------------------------------------------------------

# BB-SSL n = 100. rho = 0.7, sigma2 = 1, run = 1/1.

out.bbssl.1$MSE$MSE.bbssl |> mean(x = _, na.rm = T) |> round(x = _, 3)
out.bbssl.1$ACC$ACC.bbssl |> mean(x = _, na.rm = T) |> round(x = _, 3)
out.bbssl.1$MSPE$MSPE.bbssl |> mean(x = _, na.rm = T) |> round(x = _, 3)
out.bbssl.1$ELPPD$ELPPD.bbssl |> mean(x = _, na.rm = T) |> round(x = _, 3)
out.bbssl.1$TIME$TIME.bbssl |> mean(x = _, na.rm = T) |> round(x = _, 3)
