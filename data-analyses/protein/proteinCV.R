# ------------------------------------------------------------------------------
# Cross-validation analysis of the protein data
# ------------------------------------------------------------------------------

library(bnplasso)
library(BAS) |> suppressPackageStartupMessages() |> suppressWarnings()
source("./R/horseshoe.R"); source("./R/bbssl.R"); source("./R/metrics.R")

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
# 10-fold CV
# ------------------------------------------------------------------------------

CV.folds <- 10

# Prepare the returns ----------------------------------------------------------
shrinkage.priors <- c("bnp", paste0("bbssl", 1:5), "hs", "bayes", "adapt")
MSPE.CV <- lapply(seq_len(length(shrinkage.priors)), \(m) rep(NA, CV.folds))
ELPPD.CV <- lapply(seq_len(length(shrinkage.priors)), \(m) rep(NA, CV.folds))
names(MSPE.CV) <- paste0("MSPE.CV.", shrinkage.priors)
names(ELPPD.CV) <- paste0("ELPPD.CV.", shrinkage.priors)

# Run 10-fold CV ---------------------------------------------------------------
set.seed(1)
folds <- sample(rep(seq_len(CV.folds), each = n.full/CV.folds))

for (k in seq_len(CV.folds)) {
  # Train and held-out sets ----------------------------------------------------
  idx.train <- folds != k
  idx.test <- folds == k
  X <- X.full[idx.train,]
  y <- y.full[idx.train]
  X.test <- X.full[idx.test,] 
  y.test <- y.full[idx.test]
  n <- length(y)
  # Fit the models -------------------------------------------------------------
  set.seed(1)
  # BNP Lasso
  out.bnp <- bnplasso::bnplasso.lm(X, y, a = 0.1, b = 1e-5, alpha = 0.8)
  # BBSSL1
  out.bbssl1 <- bbssl(X, y, lambda1 = n, lambda2 = 1)
  # BBSSL2
  out.bbssl2 <- bbssl(X, y, lambda1 = 50, lambda2 = 1)
  # BBSSL3
  out.bbssl3 <- bbssl(X, y, lambda1 = 20, lambda2 = 1)
  # BBSSL4
  out.bbssl4 <- bbssl(X, y, lambda1 = 10, lambda2 = 1)
  # BBSSL5
  out.bbssl5 <- bbssl(X, y, lambda1 = 10, lambda2 = 0.15)
  # Horseshoe
  out.hs <- hs_prior(X, y)
  # B-Lasso
  out.bayes <- bnplasso::blasso.lm(
    X = X, y = y, a = n^4 / 1e-5, b = n^2 / 1e-5, 
    variance.prior.type = "independent"
  )
  # BA-Lasso
  out.adapt <- bnplasso::balasso.lm(
    X = X, y = y, a = n^4 / 1e-5, b = n^2 / 1e-5, 
    variance.prior.type = "independent"
  )
  # Point estimates ------------------------------------------------------------
  # beta.hat
  beta.hat.bnp <- bnplasso::point.estimates(out.bnp, retain = "mean")
  beta.hat.bbssl1 <- bnplasso::point.estimates(out.bbssl1, retain = "mean")
  beta.hat.bbssl2 <- bnplasso::point.estimates(out.bbssl2, retain = "mean")
  beta.hat.bbssl3 <- bnplasso::point.estimates(out.bbssl3, retain = "mean")
  beta.hat.bbssl4 <- bnplasso::point.estimates(out.bbssl4, retain = "mean")
  beta.hat.bbssl5 <- bnplasso::point.estimates(out.bbssl5, retain = "mean")
  beta.hat.hs <- bnplasso::point.estimates(out.hs, retain = "mean")
  beta.hat.bayes <- bnplasso::point.estimates(out.bayes, retain = "mean")
  beta.hat.adapt <- bnplasso::point.estimates(out.adapt, retain = "mean")
  # y.hat
  y.hat.bnp <- X.test %*% beta.hat.bnp
  y.hat.bbssl1 <- X.test %*% beta.hat.bbssl1
  y.hat.bbssl2 <- X.test %*% beta.hat.bbssl2
  y.hat.bbssl3 <- X.test %*% beta.hat.bbssl3
  y.hat.bbssl4 <- X.test %*% beta.hat.bbssl4
  y.hat.bbssl5 <- X.test %*% beta.hat.bbssl5
  y.hat.hs <- X.test %*% beta.hat.hs
  y.hat.bayes <- X.test %*% beta.hat.bayes
  y.hat.adapt <- X.test %*% beta.hat.adapt
  # Assessment  ----------------------------------------------------------------
  # MSPE
  MSPE.CV$MSPE.CV.bnp[k] <- mspe(y.hat.bnp, y.test)
  MSPE.CV$MSPE.CV.bbssl1[k] <- mspe(y.hat.bbssl1, y.test)
  MSPE.CV$MSPE.CV.bbssl2[k] <- mspe(y.hat.bbssl2, y.test)
  MSPE.CV$MSPE.CV.bbssl3[k] <- mspe(y.hat.bbssl3, y.test)
  MSPE.CV$MSPE.CV.bbssl4[k] <- mspe(y.hat.bbssl4, y.test)
  MSPE.CV$MSPE.CV.bbssl5[k] <- mspe(y.hat.bbssl5, y.test)
  MSPE.CV$MSPE.CV.hs[k] <- mspe(y.hat.hs, y.test)
  MSPE.CV$MSPE.CV.bayes[k] <- mspe(y.hat.bayes, y.test)
  MSPE.CV$MSPE.CV.adapt[k] <- mspe(y.hat.adapt, y.test)
  # ELPPD
  ELPPD.CV$ELPPD.CV.bnp[k] <- bnplasso::elppd(out.bnp, X.test, y.test) 
  ELPPD.CV$ELPPD.CV.bbssl1[k] <- bnplasso::elppd(out.bbssl1, X.test, y.test)
  ELPPD.CV$ELPPD.CV.bbssl2[k] <- bnplasso::elppd(out.bbssl2, X.test, y.test)
  ELPPD.CV$ELPPD.CV.bbssl3[k] <- bnplasso::elppd(out.bbssl3, X.test, y.test)
  ELPPD.CV$ELPPD.CV.bbssl4[k] <- bnplasso::elppd(out.bbssl4, X.test, y.test)
  ELPPD.CV$ELPPD.CV.bbssl5[k] <- bnplasso::elppd(out.bbssl5, X.test, y.test)
  ELPPD.CV$ELPPD.CV.hs[k] <- bnplasso::elppd(out.hs, X.test, y.test)
  ELPPD.CV$ELPPD.CV.bayes[k] <- bnplasso::elppd(out.bayes, X.test, y.test)
  ELPPD.CV$ELPPD.CV.adapt[k] <- bnplasso::elppd(out.adapt, X.test, y.test)
  print(paste("Fold: ", k))
}

# Export the results -----------------------------------------------------------
save(MSPE.CV, file = "./data-analyses/protein/MSPE_CV.Rdata")
save(ELPPD.CV, file = "./data-analyses/protein/ELPPD_CV.Rdata")

# Load the results -------------------------------------------------------------
load("./data-analyses/protein/MSPE_CV.Rdata")
load("./data-analyses/protein/ELPPD_CV.Rdata")

# Box-plots --------------------------------------------------------------------
prior.names <- c("BNP-L", paste0("BB-SSL", 1:5), "H-S", "B-L", "B-AL")
names(MSPE.CV) <- names(ELPPD.CV) <- prior.names

par(mfrow = c(2, 1))
boxplot(
  MSPE.CV, col = "#08306B", pch = 16, ylab = "CV  MSPE", 
  main = "CV  MSPE", medcol = "lightgrey"
)
boxplot(
  ELPPD.CV, col = "darkred", pch = 16, ylab = "CV  ELPPD", 
  main = "CV  ELPPD", medcol = "lightgrey"
)
# Size of plot: 10.17 x 7 - Landscape mode

# Table of results -------------------------------------------------------------
out <- rbind(
  sapply(MSPE.CV, mean),
  sapply(ELPPD.CV, mean)
) |> round(x = _, digits = 3L)

print(out)
    