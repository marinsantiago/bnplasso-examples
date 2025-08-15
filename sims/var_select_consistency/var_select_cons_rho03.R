# ------------------------------------------------------------------------------
# Posterior concentration
# ------------------------------------------------------------------------------

library(bnplasso)
source("./R/metrics.R"); source("./R/gen_data.R"); source("./R/vi.R")

rho <- 0.3
reps <- 80L
n.vals <- seq(100, 2000, 100)
tot.n <- length(n.vals)
p <- 200

# Prepare the returns ----------------------------------------------------------
ACC.Tang <- matrix(NA, nrow = reps, ncol = tot.n)
ACC.Li_Lin <- matrix(NA, nrow = reps, ncol = tot.n)
ACC.bnp_clust <- matrix(NA, nrow = reps, ncol = tot.n)
FP.Tang <- FP.Li_Lin <- FP.bnp_clust <- matrix(NA, nrow = reps, ncol = tot.n)
FN.Tang <- FN.Li_Lin <- FN.bnp_clust <- matrix(NA, nrow = reps, ncol = tot.n)
CONCENTRATION <- matrix(NA, nrow = reps, ncol = tot.n)

# Variable selection procedure from Tang et al. (2018) -------------------------
beta.hat.Tang <- function(X, y, n, p, object) {
  if (p > n) return(rep(NA, p))
  beta.ols <- c(lm(y ~ X - 1)$coef)
  beta.pm <- bnplasso::point.estimates(object, "post.mean", "mean")
  rule <- ifelse(abs(beta.pm / beta.ols) > 0.5, 1.0, 0.0)
  beta.pm * rule
}

# Variable selection based on the clustering results from the bnplasso ---------
beta.hat.bnp_clust <- function(object) {
  clusts <- vi(object$Post.clust_idx)
  clusts.idx <- unique(clusts)
  lambdas <- sapply(clusts.idx, \(c) mean(object$Post.lambda2[,clusts == c]))
  beta.pm <- bnplasso::point.estimates(object, "post.mean", "mean")
  rule <- ifelse(clusts == which.max(lambdas), 0.0, 1.0)
  beta.pm * rule
}

# Run the simulations ----------------------------------------------------------
for (i in seq_len(tot.n)) {
  n <- n.vals[i]
  for (r in seq_len(reps)) {
    # Generate the data --------------------------------------------------------
    set.seed(r)
    data <- gen_data(n = n, rho = rho)
    true.beta <- data$true.beta; true.sigma2 <- data$true.sigma2 # True params
    X <- data$X; X.test <- data$X.test; y <- data$y; y.test <- data$y.test
    rm(data); gc()
    # Fit the model ------------------------------------------------------------
    set.seed(1)
    out.bnp <- bnplasso::bnplasso.lm(X, y, a = 0.1, b = 0.01/n, alpha = 0.6)
    # Compare various variable selection procedures ----------------------------
    estimates.bnp.clsuts <- beta.hat.bnp_clust(out.bnp)
    estimates.Li_Lin <- bnplasso::point.estimates(out.bnp, retain = "mean")
    estimates.Tang <- beta.hat.Tang(X, y, n, p, out.bnp)
    ACC.Tang[r, i] <- sel.accuracy(estimates.Tang, true.beta)
    ACC.Li_Lin[r, i] <- sel.accuracy(estimates.Li_Lin, true.beta)
    ACC.bnp_clust[r, i] <- sel.accuracy(estimates.bnp.clsuts, true.beta)
    # False positives ----------------------------------------------------------
    FP.Tang[r, i] <- sum(estimates.Tang[21:200] != 0) / 180
    FP.Li_Lin[r, i] <- sum(estimates.Li_Lin[21:200] != 0) / 180
    FP.bnp_clust[r, i] <- sum(estimates.bnp.clsuts[21:200] != 0) / 180
    # False negatives ----------------------------------------------------------
    FN.Tang[r, i] <- sum(estimates.Tang[1:20] == 0) / 20
    FN.Li_Lin[r, i] <- sum(estimates.Li_Lin[1:20] == 0) / 20
    FN.bnp_clust[r, i] <- sum(estimates.bnp.clsuts[1:20] == 0) / 20
    # Posterior concentration --------------------------------------------------
    p <- out.bnp$n.preds
    dd <- out.bnp$n.draws
    concentration <- rep(NA, p)
    for (j in seq_len(p)) {
      beta.post <- out.bnp$Post.beta[,j]
      true.b <- true.beta[j]
      ss <- sum(((true.b-0.1) < beta.post) & ((true.b+0.1) > beta.post))
      concentration[j] <- ss/dd
    }
    CONCENTRATION[r, i] <- mean(concentration)
    paste("n =", n, "rep =", r) |> print()
  }
}

out_var_select_cons_rho03 <- list(
  ACC.Tang = ACC.Tang, ACC.Li_Lin = ACC.Li_Lin, ACC.bnp_clust = ACC.bnp_clust,
  FP.Tang = FP.Tang, FP.Li_Lin = FP.Li_Lin, FP.bnp_clust = FP.bnp_clust,
  FN.Tang = FN.Tang, FN.Li_Lin = FN.Li_Lin, FN.bnp_clust = FN.bnp_clust,
  CONCENTRATION = CONCENTRATION
)

# Export the results -----------------------------------------------------------
save(
  out_var_select_cons_rho03, 
  file = "./sims/var_select_consistency/var_select_cons_rho03.Rdata"
)
