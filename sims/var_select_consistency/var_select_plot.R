# ------------------------------------------------------------------------------
# Posterior concentration results and plots
# ------------------------------------------------------------------------------

load("./sims/var_select_consistency/var_select_cons_rho03.Rdata")

# Posterior concentration ------------------------------------------------------
n.vals <- seq(100, 2000, 100)
conc.out <- out_var_select_cons_rho03$CONCENTRATION
conc.med <- apply(conc.out, 2, median)
conc.upq <- apply(conc.out, 2, \(x) quantile(x, 0.75))
conc.lwq <- apply(conc.out, 2, \(x) quantile(x, 0.25))

plot(
  n.vals, conc.upq, type = "n", xlab = "Sample size", ylab = "Concentration",
  main = "Posterior concentration", cex = 1.2, cex.main = 1.55,
  cex.lab = 1.35, cex.axis = 1.2, font.main = 2, ylim = c(0.90, 1)
)
polygon(
  x = c(n.vals, rev(n.vals)), y = c(conc.lwq, rev(conc.upq)),
  col = "#FF00004D", border = NA
)
lines(x = n.vals, conc.med, type = "b", pch = 16, lwd = 2, col = "darkred")


# Selection Accuracy, False Positives and False negatives ----------------------
par(mfrow = c(3, 1))

out <- out_var_select_cons_rho03

mean.acc.Tang <- colMeans(out$ACC.Tang)
mean.acc.Li_lin <- colMeans(out$ACC.Li_Lin)
mean.acc.bnp <- colMeans(out$ACC.bnp_clust)

plot(
  x = n.vals, y = mean.acc.Tang, type = "b", ylim = c(0.75, 1), pch = 15,
  xlab = "Sample size", ylab = "Selection accuracy", lwd = 2, 
  main = "Variable selection accuracy", cex = 1.2, lty = 2, cex.main = 1.55, 
  cex.lab = 1.35, cex.axis = 1.2, font.main = 2, col = "#08306B"
)
lines(
  x = n.vals, y = mean.acc.Li_lin, type = "b", ylim = c(0.75, 1),
  col = "darkred", pch = 17, lwd = 2, cex = 1.2
)
#lines(
#  x = n.vals, y = mean.acc.bnp, type = "b", ylim = c(0.75, 1),
#  pch = 16, lwd = 2, cex = 1.2, lty = 3
#)
legend(
  x = "bottomleft", c("Tang et al. (2018)", "Li and Lin (2010)"), 
  pch = c(15, 17), cex = 1.2, col = c("#08306B", "darkred"), bty = "n",
  lty = c(2, 1, 3), lwd = 1.4, y.intersp = 0.6
)

mean.fp.Tang <- colMeans(out$FP.Tang)
mean.fp.Li_lin <- colMeans(out$FP.Li_Lin)
mean.fp.bnp <- colMeans(out$FP.bnp_clust)

plot(
  x = n.vals, y = mean.fp.Tang, type = "b", pch = 15, ylim = c(-0.003, 0.1),
  xlab = "Sample size", ylab = "False positive", lwd = 2, 
  main = "False positive rate", cex = 1.2, lty = 2, cex.main = 1.55, 
  cex.lab = 1.35, cex.axis = 1.2, font.main = 2, col = "#08306B"
)
lines(
  x = n.vals, y = mean.fp.Li_lin, type = "b", ylim = c(0.75, 1),
  col = "darkred", pch = 17, lwd = 2, cex = 1.2
)
#lines(
#  x = n.vals, y = mean.fp.bnp, type = "b", ylim = c(0.75, 1),
#  pch = 16, lwd = 2, cex = 1.2, lty = 3
#)
legend(
  x = "topleft", c("Tang et al. (2018)", "Li and Lin (2010)"), 
  pch = c(15, 17), cex = 1.2, col = c("#08306B", "darkred"), bty = "n",
  lty = c(2, 1, 3), lwd = 1.4, y.intersp = 0.6
)


mean.fn.Tang <- colMeans(out$FN.Tang)
mean.fn.Li_lin <- colMeans(out$FN.Li_Lin)
mean.fn.bnp <- colMeans(out$FN.bnp_clust)


plot(
  x = n.vals, y = mean.fn.Tang, type = "b", pch = 15,
  xlab = "Sample size", ylab = "False negative", lwd = 2, 
  main = "False negative rate", cex = 1.2, lty = 2, cex.main = 1.55, 
  cex.lab = 1.35, cex.axis = 1.2, font.main = 2, col = "#08306B"
)
lines(
  x = n.vals, y = mean.fn.Li_lin, type = "b", ylim = c(0.75, 1),
  col = "darkred", pch = 17, lwd = 2, cex = 1.2
)
#lines(
#  x = n.vals, y = mean.fn.bnp, type = "b", ylim = c(0.75, 1),
#  pch = 16, lwd = 2, cex = 1.2, lty = 3
#)
legend(
  #x = "topright", 
  x = 1600, y = 0.14,
  c("Tang et al. (2018)", "Li and Lin (2010)"), 
  pch = c(15, 17), cex = 1.2, col = c("#08306B", "darkred"), bty = "n",
  lty = c(2, 1, 3), lwd = 1.4, y.intersp = 0.6
)

# Size of plot: 10 x 8.5 portrait mode
