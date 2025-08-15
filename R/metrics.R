# ------------------------------------------------------------------------------
# Assessment metrics
# ------------------------------------------------------------------------------

mse <- \(beta.hat, beta.true) mean((beta.hat - beta.true)^2)
mspe <- \(y.hat, y.true) mean((y.hat - y.true)^2)
sel.accuracy <- \(beta.hat, beta.true) mean(sign(beta.hat) == sign(beta.true))
