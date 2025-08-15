# ------------------------------------------------------------------------------
# Fit linear models using the horseshoe prior
# ------------------------------------------------------------------------------

library(horseshoe) |> suppressPackageStartupMessages() |> suppressWarnings()

hs_prior <- function(X, y, max.iters = 6000, burn.in = 1000) {
  
  # Fit the model via MCMC -----------------------------------------------------
  start <- Sys.time()
  hs_fit <- horseshoe::horseshoe(
    y = y, X = X, burn = burn.in, nmc = max.iters - burn.in, 
  )
  end <- Sys.time()
  elapsed <- end - start
  
  # Prepare the returns --------------------------------------------------------
  hs_out <- list(
    Post.beta = t(hs_fit$BetaSamples), 
    # By default, horseshoe fixes sigma2 == 1
    Post.sigma2 = rep(1, nrow(t(hs_fit$BetaSamples))),
    elapsed = elapsed, intercept = FALSE, max.iters = max.iters, 
    burn.in = burn.in, thin = 1L, n.obs = length(y),
    n.preds = ncol(X), n.draws = nrow(t(hs_fit$BetaSamples))
  )
  
  # Posterior predictive fitted values and residuals ---------------------------
  # Pre-compute all linear predictors (without the intercept)
  linPreds <- tcrossprod(X, hs_out$Post.beta)
  # Note: In "post_pred_fits" and "post_pred_res", each row corresponds to an
  # MCMC draw and each column to an observation.
  n.draws <- hs_out$n.draws
  n.obs <- hs_out$n.obs
  post_pred_fits <- matrix(data = NA, nrow = n.draws, ncol = n.obs)
  post_pred_res <- post_pred_fits
  post_sigmas <- sqrt(abs(hs_out$Post.sigma2) + 1e-16)
  for (s in seq_len(n.draws)) {
    fit_val_s <- linPreds[,s] + rnorm(n.obs, 0, post_sigmas[s])
    post_pred_fits[s,] <- fit_val_s
    post_pred_res[s,] <- y - fit_val_s
  }
  rm(linPreds, post_sigmas, fit_val_s); gc()
  hs_out[["X"]] <- X
  hs_out[["y"]] <- y
  hs_out[["post.pred.fitted.values"]] <- post_pred_fits
  hs_out[["post.pred.residuals"]] <- post_pred_res
  # Make the class consistent with the package "bnplasso"
  class(hs_out) <- "lmBayes"
  hs_out
}
