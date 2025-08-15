# ------------------------------------------------------------------------------
# Analysis of the protein data
# ------------------------------------------------------------------------------

library(bnplasso)
library(BAS) |> suppressPackageStartupMessages() |> suppressWarnings()
library(FKSUM) |> suppressPackageStartupMessages() |> suppressWarnings()
source("./R/horseshoe.R"); source("./R/bbssl.R"); source("./R/vi.R")

data("protein")

X.full <- model.matrix(
  prot.act4 ~ 
    buf + buf * pH + buf*NaCl + buf*con + buf*ra + buf*det + buf*MgCl2 + 
    buf*temp + pH + pH*NaCl + pH*con + pH*ra + pH*det + pH*MgCl2 + pH*temp +
    NaCl + NaCl*con + NaCl*ra + NaCl*det + NaCl*MgCl2 + NaCl*temp + con +
    con*ra + con*det + con*MgCl2 +con*temp + ra + ra*det + ra*MgCl2 +
    ra*temp + det + det*MgCl2 + det*temp + MgCl2 + MgCl2*temp + I(NaCl^2) +
    I(pH^2) + I(con^2) + I(temp^2),
  data = protein
)

X.full <- X.full[,-1]
X.full <- scale(X.full, scale = FALSE)
max(abs(colMeans(X.full))) < 1e-10
sum(cor(X.full)[upper.tri(cor(X.full))] > 0.95) # 17 Highly correlated features
max(cor(X.full)[upper.tri(cor(X.full))])

y.full <- protein$prot.act4 - mean(protein$prot.act4)
mean(y.full)

dim.bf <- dim(X.full)
(n.full <- dim.bf[1]) # Number of observations
(p <- dim.bf[2]) # Number of potential features

# ------------------------------------------------------------------------------
# Run different shrinkage priors
# ------------------------------------------------------------------------------

# Benchmark - BAS --------------------------------------------------------------
set.seed(1)
bas_fit <- bas.lm(y.full ~ X.full, prior = "JZS", method = "BAS")
# Look at Posterior Inclusion Probabilities (PIP) from BAS 
pip <- summary(bas_fit)[c(2:(p+1)),1]
bas_index <- which(pip > 0.5)
plot(pip, type = "h")
abline(h = 0.5, col = 2)
# Important features recovered by BAS
(variables.bas <- colnames(X.full)[bas_index])
length(bas_index)

# BNP-Lasso --------------------------------------------------------------------
set.seed(1)
# a = 0.1, b = 1e-5
out.bnp <- bnplasso::bnplasso.lm(X.full, y.full, a = 0.1, b = 1e-5, alpha = 0.8)
bnp_index <- which(bnplasso::point.estimates(out.bnp, retain = "mean") != 0)
# Important features recovered by BNP-Lasso
(variables.bnp <- colnames(X.full)[bnp_index]) 
length(bnp_index)

# BBSSL1 -----------------------------------------------------------------------
set.seed(1)
out.bbssl1 <- bbssl(X.full, y.full, lambda1 = n.full, lambda2 = 1)
bbssl1_index <- which(bnplasso::point.estimates(out.bbssl1, retain="mean") != 0)
colnames(X.full)[bbssl1_index] # Important features recovered by BBSSL1

# BBSSL2 -----------------------------------------------------------------------
set.seed(1)
out.bbssl2 <- bbssl(X.full, y.full, lambda1 = 50, lambda2 = 1)
bbssl2_index <- which(bnplasso::point.estimates(out.bbssl2, retain="mean") != 0)
colnames(X.full)[bbssl2_index] # Important features recovered by BBSSL2

# BBSSL3 -----------------------------------------------------------------------
set.seed(1)
out.bbssl3 <- bbssl(X.full, y.full, lambda1 = 20, lambda2 = 1)
bbssl3_index <- which(bnplasso::point.estimates(out.bbssl3, retain="mean") != 0)
colnames(X.full)[bbssl3_index] # Important features recovered by BBSSL3

# BBSSL4 -----------------------------------------------------------------------
set.seed(1)
out.bbssl4 <- bbssl(X.full, y.full, lambda1 = 10, lambda2 = 1)
bbssl4_index <- which(bnplasso::point.estimates(out.bbssl4, retain="mean") != 0)
colnames(X.full)[bbssl4_index] # Important features recovered by BBSSL4

# BBSSL5 -----------------------------------------------------------------------
set.seed(1)
out.bbssl5 <- bbssl(X.full, y.full, lambda1 = 10, lambda2 = 0.15)
bbssl5_index <- which(bnplasso::point.estimates(out.bbssl5, retain="mean") != 0)
colnames(X.full)[bbssl5_index] # Important features recovered by BBSSL5

# Horseshoe --------------------------------------------------------------------
set.seed(1)
out.hs <- hs_prior(X.full, y.full)
hs_index <- which(bnplasso::point.estimates(out.hs, retain = "mean") != 0)
colnames(X.full)[hs_index] # Important features recovered by Horseshoe

# Tail behavior so that it mimics the spike from Moran et al. (2019).
a.param.lasso <- n.full^4 / 1e-5
b.param.lasso <- n.full^2 / 1e-5

# B-Lasso ----------------------------------------------------------------------
set.seed(1)
out.bayes <- bnplasso::blasso.lm(
  X = X.full, y = y.full, a = a.param.lasso, b = b.param.lasso, 
  variance.prior.type = "independent"
)
bayes_index <- which(bnplasso::point.estimates(out.bayes, retain = "mean") != 0)
colnames(X.full)[bayes_index] # Important features recovered by B-Lasso
length(bayes_index)

# BA-Lasso ---------------------------------------------------------------------
set.seed(1)
out.adapt <- bnplasso::balasso.lm(
  X = X.full, y = y.full, a = a.param.lasso, b = b.param.lasso, 
  variance.prior.type = "independent"
)
adapt_index <- which(bnplasso::point.estimates(out.adapt, retain = "mean") != 0)
colnames(X.full)[adapt_index] # Important features recovered by BA-Lasso
length(adapt_index)

# Correlation between the important features identified by BAS and BNP-Lasso
cor(X.full[,colnames(X.full) %in% union(variables.bnp, variables.bas)]) 
cor(X.full[,"detT"], X.full[,"pH:detT"]) |> round(x = _, digits = 3)  # 0.988
cor(X.full[,"detN"], X.full[,"con:detN"]) |> round(x = _, digits = 3) # 0.735

# ------------------------------------------------------------------------------
# Posterior predictive fitted values
# ------------------------------------------------------------------------------

BAS_fitted <- fitted.values(bas_fit)

par(mfrow = c(3, 3))
for (i in seq_len(9)) {
  dens.bnp <- FKSUM::fk_density(out.bnp$post.pred.fitted.values[,i])
  dens.bnp.y <- dens.bnp$y
  dens.bnp.x <- dens.bnp$x
  plot(
    x = dens.bnp.x, y = dens.bnp.y, type = "l", col = "#08306B", lwd = 2,
    main = paste("Obs: ", i), xlab = expression(hat(y)), ylab = "Density"
  )
  abline(v = BAS_fitted[i], lty = 3, lwd = 2, col = "darkred")
  legend(
    x = "topleft", c("BNP-L", "BAS"), lty = c(1, 3), lwd = 2,
    col = c("#08306B", "darkred"), bty = "n", y.intersp = 0.5
  )
} # Size of plot: 10 x 8.5 portrait mode

# Fitted densities -------------------------------------------------------------
d.eval <- \(x, x_eval) FKSUM::fk_density(x = x, x_eval = x_eval)$y

yhat.bnp <- out.bnp$post.pred.fitted.values
yhat.bbssl1 <- out.bbssl1$post.pred.fitted.values
yhat.bbssl2 <- out.bbssl2$post.pred.fitted.values
yhat.bbssl3 <- out.bbssl3$post.pred.fitted.values
yhat.bbssl4 <- out.bbssl4$post.pred.fitted.values
yhat.bbssl5 <- out.bbssl5$post.pred.fitted.values
yhat.hs <- out.hs$post.pred.fitted.values
yhat.bayes <- out.bayes$post.pred.fitted.values
yhat.adapt <- out.adapt$post.pred.fitted.values

n <- n.full
fit.bnp <- sapply(1:n, \(i) d.eval(yhat.bnp[,i], BAS_fitted[i])) |> mean()
fit.bbssl1 <- sapply(1:n, \(i) d.eval(yhat.bbssl1[,i], BAS_fitted[i])) |> mean()
fit.bbssl2 <- sapply(1:n, \(i) d.eval(yhat.bbssl2[,i], BAS_fitted[i])) |> mean()
fit.bbssl3 <- sapply(1:n, \(i) d.eval(yhat.bbssl3[,i], BAS_fitted[i])) |> mean()
fit.bbssl4 <- sapply(1:n, \(i) d.eval(yhat.bbssl4[,i], BAS_fitted[i])) |> mean()
fit.bbssl5 <- sapply(1:n, \(i) d.eval(yhat.bbssl5[,i], BAS_fitted[i])) |> mean()
fit.hs <- sapply(1:n, \(i) d.eval(yhat.hs[,i], BAS_fitted[i])) |> mean()
fit.bayes <- sapply(1:n, \(i) d.eval(yhat.bayes[,i], BAS_fitted[i])) |> mean()
fit.adapt <- sapply(1:n, \(i) d.eval(yhat.adapt[,i], BAS_fitted[i])) |> mean()

# Print the results
out <- c(
  fit.bnp, fit.bbssl1, fit.bbssl2, fit.bbssl3, fit.bbssl4, 
  fit.bbssl5, fit.hs, fit.bayes, fit.adapt
)
names(out) <- c(
  "fit.bnp", "fit.bbssl1", "fit.bbssl2", "fit.bbssl3", "fit.bbssl4", 
  "fit.bbssl5", "fit.hs", "fit.bayes", "fit.adapt"
)
round(out, 3) |> print()
