# ------------------------------------------------------------------------------
# Bayesian bootstrap spike-and-slab Lasso
# ------------------------------------------------------------------------------

library(BBSSL) |> suppressPackageStartupMessages() |> suppressWarnings()

bbssl <- function(X, y, lambda1 = length(y), lambda2 = 1, # Default lambda vals
                  pi_a = 1, pi_b = ncol(X),
                  annealing = TRUE, # Whether to use annealing strategy
                  max.iters = 3000L) {
  
  # Pre-compute constants
  p <- ncol(X)
  lambda_vec <- c(lambda1, lambda2) # Make sure that lambda1 >> lambda2
  
  # Fit the model --------------------------------------------------------------
  start <- Sys.time()
  bbssl_fit <- BBSSL::BB_SSL(
    y = y, X = X, lambda = lambda_vec, NSample = max.iters,
    a = pi_a, b = pi_b, initial.beta = rep(0, p), burn.in	= annealing
  ) |> suppressWarnings()
  end <- Sys.time()
  elapsed <- difftime(end, start, units = "secs")
  
  # Prepare the returns --------------------------------------------------------
  bbssl_out <- list(
    # By default, BBSSL fixes sigma2 == 1
    Post.beta = bbssl_fit$beta, Post.sigma2 = rep(1, nrow(bbssl_fit$beta)),
    Post.gamma = bbssl_fit$gamma, elapsed = elapsed, intercept = FALSE, 
    max.iters = max.iters, n.obs = length(y), n.preds = p, n.draws = max.iters
  )
  
  # Posterior predictive fitted values and residuals ---------------------------
  # Pre-compute all linear predictors (without the intercept)
  linPreds <- tcrossprod(X, bbssl_out$Post.beta)
  # Note: In "post_pred_fits" and "post_pred_res", each row corresponds to an
  # BBSSL draw and each column to an observation.
  n.draws <- bbssl_out$n.draws
  n.obs <- bbssl_out$n.obs
  post_pred_fits <- matrix(data = NA, nrow = n.draws, ncol = n.obs)
  post_pred_res <- post_pred_fits
  post_sigmas <- sqrt(abs(bbssl_out$Post.sigma2) + 1e-16)
  for (s in seq_len(n.draws)) {
    fit_val_s <- linPreds[,s] + rnorm(n.obs, 0, post_sigmas[s])
    post_pred_fits[s,] <- fit_val_s
    post_pred_res[s,] <- y - fit_val_s
  }
  rm(linPreds, post_sigmas, fit_val_s); gc()
  bbssl_out[["X"]] <- X
  bbssl_out[["y"]] <- y
  bbssl_out[["post.pred.fitted.values"]] <- post_pred_fits
  bbssl_out[["post.pred.residuals"]] <- post_pred_res
  # Make the class consistent with the package "bnplasso"
  class(bbssl_out) <- "lmBayes"
  bbssl_out
}
