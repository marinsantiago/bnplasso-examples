# ------------------------------------------------------------------------------
# Simulation setting with n = 500 and rho = 0.5, on MCMC-based models
# ------------------------------------------------------------------------------

library(bnplasso)
source("./R/metrics.R"); source("./R/cls_plt.R"); source("./R/io_results.R")
source("./R/gen_data.R"); source("./R/horseshoe.R"); source("./R/vi.R")

# Simulation settings ----------------------------------------------------------
n <- 500
rho <- 0.5
reps <- 200

# Prepare the returns ----------------------------------------------------------
MSE.bnp <- MSE.bayes <- MSE.adapt <- MSE.hs <- rep(NA, reps)
ACC.bnp <- ACC.bayes <- ACC.adapt <- ACC.hs <- rep(NA, reps)
MSPE.bnp <- MSPE.bayes <- MSPE.adapt <- MSPE.hs <- rep(NA, reps)
TIME.bnp <- TIME.bayes <- TIME.adapt <- TIME.hs <- rep(NA, reps)
ELPPD.bnp <- ELPPD.bayes <- ELPPD.adapt <- ELPPD.hs <- rep(NA, reps)
CLUSTERS.bnp <- matrix(NA, nrow = min(reps, 25), ncol = 200)

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
  # BNP-Lasso
  out.bnp <- bnplasso::bnplasso.lm(X, y, a = 0.1, b = 0.01/n, alpha = 0.6)
  # B-Lasso
  out.bayes <- bnplasso::blasso.lm(X, y, a = 0.1, b = 0.01)
  # BA-Lasso
  out.bayes.adapt <- bnplasso::balasso.lm(X, y, a = 0.1, b = 0.01)
  # Horseshoe
  out.hs <- hs_prior(X, y)
  # BBSSL
  #out.bbssl <- bbssl(X, y)
  
  # Point estimates ------------------------------------------------------------
  # beta.hat
  beta.hat.bnp <- bnplasso::point.estimates(out.bnp, retain = "mean")
  beta.hat.bayes <- bnplasso::point.estimates(out.bayes, retain = "mean")
  beta.hat.adapt <- bnplasso::point.estimates(out.bayes.adapt, retain = "mean")
  beta.hat.hs <- bnplasso::point.estimates(out.hs, retain = "mean")
  #beta.hat.bbssl <- bnplasso::point.estimates(out.bbssl, retain = "mean")
  # y.hat
  y.hat.bnp <- X.test %*% beta.hat.bnp
  y.hat.bayes <- X.test %*% beta.hat.bayes
  y.hat.adapt <- X.test %*% beta.hat.adapt
  y.hat.hs <- X.test %*% beta.hat.hs
  #y.hat.bbssl <- X.test %*% beta.hat.bbssl
  # Clustering of the regression coefficients
  if (r <= 25) CLUSTERS.bnp[r,] <- vi(out.bnp$Post.clust_idx)
  
  # Assessment metrics ---------------------------------------------------------
  # MSE
  MSE.bnp[r] <- mse(beta.hat.bnp, true.beta)
  MSE.bayes[r] <- mse(beta.hat.bayes, true.beta)
  MSE.adapt[r] <- mse(beta.hat.adapt, true.beta)
  MSE.hs[r] <- mse(beta.hat.hs, true.beta)
  #MSE.bbssl[r] <- mse(beta.hat.bbssl, true.beta)
  # ACCURACY
  ACC.bnp[r] <- sel.accuracy(beta.hat.bnp, true.beta)
  ACC.bayes[r] <- sel.accuracy(beta.hat.bayes, true.beta)
  ACC.adapt[r] <- sel.accuracy(beta.hat.adapt, true.beta)
  ACC.hs[r] <- sel.accuracy(beta.hat.hs, true.beta)
  #ACC.bbssl[r] <- sel.accuracy(beta.hat.bbssl, true.beta)
  # MSPE
  MSPE.bnp[r] <- mspe(y.hat.bnp, y.test)
  MSPE.bayes[r] <- mspe(y.hat.bayes, y.test)
  MSPE.adapt[r] <- mspe(y.hat.adapt, y.test)
  MSPE.hs[r] <- mspe(y.hat.hs, y.test)
  #MSPE.bbssl[r] <- mspe(y.hat.bbssl, y.test)
  # ELPPD
  ELPPD.bnp[r] <- bnplasso::elppd(out.bnp, X.test, y.test)
  ELPPD.bayes[r] <- bnplasso::elppd(out.bayes, X.test, y.test)
  ELPPD.adapt[r] <- bnplasso::elppd(out.bayes.adapt, X.test, y.test)
  ELPPD.hs[r] <- bnplasso::elppd(out.hs, X.test, y.test)
  #ELPPD.bbssl[r] <- bnplasso::elppd(out.bbssl, X.test, y.test)
  # ELAPSED TIME
  TIME.bnp[r] <- out.bnp$elapsed
  TIME.bayes[r] <- out.bayes$elapsed
  TIME.adapt[r] <- out.bayes.adapt$elapsed
  TIME.hs[r] <- out.hs$elapsed
  #TIME.bbssl[r] <- out.bbssl$elapsed
  
  rm(out.bnp, out.bayes, out.bayes.adapt, out.hs); gc()
  print(paste("Replication: ", r))
}

# Export the results -----------------------------------------------------------

export.results.mcmc(n, rho, sigma2 = 1, reps)

# Load the results -------------------------------------------------------------

path <- paste0(subfolder.path(n, rho, sigma2 = 1, reps, bbssl = F), "/")
paste0(path, current.set(n, rho, sigma2 = 1, reps), ".Rdata") |> load()

# Print the results ------------------------------------------------------------

cat("n =", n, "rho =", rho, "var = ", 1)

out.mcmc$MSE$MSE.bnp |> mean(x = _, na.rm = T) |> round(x = _, 4)
out.mcmc$MSE$MSE.hs |> mean(x = _, na.rm = T) |> round(x = _, 3)
out.mcmc$MSE$MSE.bayes |> mean(x = _, na.rm = T) |> round(x = _, 3)
out.mcmc$MSE$MSE.adapt |> mean(x = _, na.rm = T) |> round(x = _, 3)

out.mcmc$ACC$ACC.bnp |> mean(x = _, na.rm = T) |> round(x = _, 3)
out.mcmc$ACC$ACC.hs |> mean(x = _, na.rm = T) |> round(x = _, 3)
out.mcmc$ACC$ACC.bayes |> mean(x = _, na.rm = T) |> round(x = _, 3)
out.mcmc$ACC$ACC.adapt |> mean(x = _, na.rm = T) |> round(x = _, 3)

out.mcmc$MSPE$MSPE.bnp |> mean(x = _, na.rm = T) |> round(x = _, 3)
out.mcmc$MSPE$MSPE.hs |> mean(x = _, na.rm = T) |> round(x = _, 3)
out.mcmc$MSPE$MSPE.bayes |> mean(x = _, na.rm = T) |> round(x = _, 3)
out.mcmc$MSPE$MSPE.adapt |> mean(x = _, na.rm = T) |> round(x = _, 3)

out.mcmc$ELPPD$ELPPD.bnp |> mean(x = _, na.rm = T) |> round(x = _, 3)
out.mcmc$ELPPD$ELPPD.hs |> mean(x = _, na.rm = T) |> round(x = _, 3)
out.mcmc$ELPPD$ELPPD.bayes |> mean(x = _, na.rm = T) |> round(x = _, 3)
out.mcmc$ELPPD$ELPPD.adapt |> mean(x = _, na.rm = T) |> round(x = _, 3)

out.mcmc$TIME$TIME.bnp |> mean(x = _, na.rm = T) |> round(x = _, 3)
out.mcmc$TIME$TIME.hs |> mean(x = _, na.rm = T) |> round(x = _, 3)
out.mcmc$TIME$TIME.bayes |> mean(x = _, na.rm = T) |> round(x = _, 3)
out.mcmc$TIME$TIME.adapt |> mean(x = _, na.rm = T) |> round(x = _, 3)
