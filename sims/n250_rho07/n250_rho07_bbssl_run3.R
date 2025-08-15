# ------------------------------------------------------------------------------
# Simulation setting with n = 250 and rho = 0.7, on BB-SSL (RUN 3/3)
# ------------------------------------------------------------------------------

setwd("~/Desktop/bnplasso-examples")

library(bnplasso)
source("./R/metrics.R"); source("./R/bbssl.R")
source("./R/gen_data.R"); source("./R/io_results.R")

# Simulation settings ----------------------------------------------------------
n <- 250
rho <- 0.7
reps <- 200

# Load the results from the previous run ---------------------------------------
path <- paste0(subfolder.path(n, rho, sigma2 = 1, reps, bbssl = T), "/")
current.setting <- current.set(n, rho, sigma2 = 1, reps)
paste0(path, current.setting, "_run", 1, ".Rdata") |> load()

# Extract and prepare the returns ----------------------------------------------
MSE.bbssl <- out.bbssl.1$MSE$MSE.bbssl
ACC.bbssl <- out.bbssl.1$ACC$ACC.bbssl
MSPE.bbssl <- out.bbssl.1$MSPE$MSPE.bbssl
TIME.bbssl <- out.bbssl.1$TIME$TIME.bbssl
ELPPD.bbssl <- out.bbssl.1$ELPPD$ELPPD.bbssl

print("Starting batch 3")

# Run the simulations ----------------------------------------------------------
for (r in 141:200) { 
  
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

gc()

print("Ending batch 3")

# Export the results -----------------------------------------------------------

export.results.bbssl(n, rho, sigma2 = 1, reps, run = 1)

# Load the results -------------------------------------------------------------

path <- paste0(subfolder.path(n, rho, sigma2 = 1, reps, bbssl = T), "/")
paste0(path, current.set(n, rho, sigma2 = 1, reps), "_run", 1, ".Rdata") |> load()

# Print the results ------------------------------------------------------------

# BB-SSL n = 250, rho = 0.7, sigma2 = 1, run = 3/3, reps = 200

out.bbssl.1$MSE$MSE.bbssl |> mean(x = _, na.rm = T) |> round(x = _, 3)
out.bbssl.1$ACC$ACC.bbssl |> mean(x = _, na.rm = T) |> round(x = _, 3)
out.bbssl.1$MSPE$MSPE.bbssl |> mean(x = _, na.rm = T) |> round(x = _, 3)
out.bbssl.1$ELPPD$ELPPD.bbssl |> mean(x = _, na.rm = T) |> round(x = _, 3)
out.bbssl.1$TIME$TIME.bbssl |> mean(x = _, na.rm = T) |> round(x = _, 3)
