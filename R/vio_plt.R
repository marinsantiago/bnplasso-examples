# ------------------------------------------------------------------------------
# Violin plots of posterior mean number of clusters
# ------------------------------------------------------------------------------

library(vioplot) |> suppressPackageStartupMessages() |> suppressWarnings()

violin_plot_k <- function(means.k, rho_val, sigma2_val, n = c(100, 250, 500)) {
  m <- bquote("Setting:" ~ rho == .(rho_val) ~ "," ~ sigma^2 == .(sigma2_val))
  vioplot(
    means.k, names = paste0("n = ", n), border = NA, col = "lightblue",
    rectCol = "navy", lineCol = "navy", colMed = "navy", colMed2 = "navy",
    pchMed = 3, ylab = "", main = m, cex.main = 1.55, cex.lab = 1.35, 
    cex.axis = 1.2, font.main = 2
  )
  mtext("Posterior mean number of clusters", side = 2, line = 2.7, cex = 1.35)
}


violin_plot_k_varying_alpha <- function(means.k, a, b, n, rho_,
                                        alphas = c(0.1, 1.0, 10.0)) {
  m <- bquote(
    "Setting:" ~ rho == .(rho_) ~ "," 
    ~ "n" == .(n) ~ "," ~ "a" == .(a) ~ "," ~ "b" == .(b)
  )
  vioplot(
    means.k, names = paste0("alpha = ", alphas), 
    border = NA, col = "lightgreen", rectCol = "darkgreen", 
    lineCol = "darkgreen", colMed = "darkgreen", colMed2 = "darkgreen", 
    pchMed = 3, ylab = "", main = m, cex.main = 1.55, cex.lab = 1.35, 
    cex.axis = 1.2, font.main = 2 
  )
  mtext(
    "Posterior mean number of clusters", side = 2, line = 2.7, cex = 1.35
  )
}


violin_plot_k_varying_G0 <- function(means.k, alpha_, n, rho_,
                                     a. = c(0.01, 0.1, 1.0),
                                     b. = c(0.01, 0.1, 0.1)) {
  m <- bquote(
    "Setting:" ~ rho == .(rho_) ~ "," ~ "n" == .(n) ~ "," ~ alpha == .(alpha_) 
  )
  n <- c(
    paste("(",paste("a =",a.[1]), ", ",paste("b =", b.[1]), "/n)",sep=""),
    paste("(",paste("a =",a.[2]), ", ",paste("b =", b.[2]), "/n)",sep=""),
    paste("(",paste("a =",a.[3]), ", ",paste("b =", b.[3]), "/n)",sep="")
  )
  vioplot(
    means.k, names = n, 
    border = NA, col = "lightgreen", rectCol = "darkgreen", 
    lineCol = "darkgreen", colMed = "darkgreen", colMed2 = "darkgreen", 
    pchMed = 3, ylab = "", main = m, cex.main = 1.55, cex.lab = 1.35, 
    cex.axis = 1.2, font.main = 2 
  )
  mtext(
    "Posterior mean number of clusters", side = 2, line = 2.7, cex = 1.35
  )
}
