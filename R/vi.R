# ------------------------------------------------------------------------------
# Clustering via the posterior expected variation of information loss function
# ------------------------------------------------------------------------------

library(BNPmix) |> suppressPackageStartupMessages() |> suppressWarnings()

vi <- function(clusters.draws) {
  cls.draws <- apply(clusters.draws, 1, \(x) as.numeric(as.factor(x))) |> t()
  #object <- list(clust = cls.draws - 1)
  object <- list(clust = cls.draws)
  class(object) <- "BNPdens"
  BNPmix.out <- BNPmix::partition(object, dist = "VI")  
  min.lower.bound.expected.loss <- which.min(BNPmix.out$scores)
  BNPmix.out$partitions[min.lower.bound.expected.loss,]
}
