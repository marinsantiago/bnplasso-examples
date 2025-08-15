# ------------------------------------------------------------------------------
# Routines for co-cluster analysis
# ------------------------------------------------------------------------------

library(ggplot2) |> suppressPackageStartupMessages() |> suppressWarnings()

posterior.co_clustering.probs <- function(object, main, return.mat = FALSE,
                                          viridis.pal = "C", axis.idx = NULL) {
  
  # Pre-compute and extract constants ------------------------------------------
  S.mcmc <- object$Post.clust_idx
  dd <- object$n.draws
  p <- object$n.preds
  
  # Compute co-clustering matrix -----------------------------------------------
  coclustering <- matrix(0, nrow = p, ncol = p)
  for (iter in seq_len(dd)) {
    cluster <- S.mcmc[iter, ]
    coclustering <- coclustering + outer(cluster, cluster, FUN = "==")
  }
  coclustering.out <- coclustering / dd 
  if (return.mat) return(coclustering.out)
  
  # Gen. plot ------------------------------------------------------------------
  if (is.null(axis.idx)) {
    idxs <- seq(0, 200, 20); idxs[1] <- 1
  } else {
    idxs <- axis.idx
  }
  
  df <- reshape2::melt(coclustering.out)
  df$Prob <- df$value
  ggplot(df, aes(Var1, Var2, fill = Prob)) +
    geom_tile() +
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank()) +
    theme(panel.background = element_rect(fill = "white")) +
    ggtitle(main) + 
    labs(x = "Coefficient index", y = "Coefficient index") +
    scale_x_continuous(breaks = idxs, expand = c(0, 0)) +
    scale_y_continuous(breaks = idxs, expand = c(0, 0)) +
    viridis::scale_fill_viridis(name = "Prob.", option = viridis.pal) +
    theme(axis.title.x = element_text(margin = margin(t = 8))) +
    theme(axis.title.y = element_text(margin = margin(r = 8))) +
    theme(axis.title = element_text(size = 27)) +
    theme(axis.text = element_text(size = 27)) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 32)) +
    theme(legend.title = element_text(margin = margin(b = 27), size = 28)) +
    theme(legend.text = element_text(size = 26)) + 
    theme(legend.key.size = unit(2.3, "lines"))
}


co_clustering <- function(S, main, return.mat = FALSE, 
                          viridis.pal = "C", axis.idx = NULL) {
  
  # Pre-compute constants ------------------------------------------------------
  p <- length(S)
  
  # Compute co-clustering matrix -----------------------------------------------
  coclusts <- outer(S, S, FUN = "==") * 1
  if (return.mat) return(coclusts)
  
  # Gen. plot ------------------------------------------------------------------
  lower.col <- viridis::viridis(9, option = viridis.pal)[1]
  upper.col <- viridis::viridis(9, option = viridis.pal)[9]
  if (is.null(axis.idx)) {
    idx <- rep("", p)
    idx.seq <- seq(0, p, 20); idx.seq[1] <- 1
    idx[idx.seq] <- as.character(idx.seq)
  } else {
    idx <- rep("", p)
    idx[axis.idx] <- as.character(axis.idx)
  }
  df <- reshape2::melt(coclusts)
  df$`Co-cluster` <- factor(
    ifelse(df$value == 1, "yes", "no"), 
    levels = c("yes", "no")
  )
  ggplot(df, aes(Var1, Var2, fill = `Co-cluster`)) + 
    geom_tile() +
    scale_fill_manual(values = c("yes" = upper.col, "no" = lower.col)) + 
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank()) +
    theme(panel.background = element_rect(fill = "white")) +
    ggtitle(main) + 
    labs(x = "Coefficient index", y = "Coefficient index") +
    scale_x_discrete(limits = factor(1:p), labels = idx, breaks = idx) + 
    scale_y_discrete(limits = factor(1:p), labels = idx, breaks = idx) + 
    theme(axis.title.x = element_text(margin = margin(t = 8))) +
    theme(axis.title.y = element_text(margin = margin(r = 8))) +
    theme(axis.text = element_text(size = 10)) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
    theme(legend.title = element_text(margin = margin(b = 15)))
}
