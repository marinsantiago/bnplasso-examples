# ------------------------------------------------------------------------------
# Import and export simulation results
# ------------------------------------------------------------------------------

current.set <- \(n, rho, sigma2, reps) {
  current.setting <- paste("n", n, "_rho", rho, "_var", sigma2, "_reps", reps)
  current.setting |> gsub("\\.", "", x = _) |> gsub(" ", "", x = _)
}

chck.folder <- \(path) file.exists(path) && file.info(path)$isdir

create.results.folder <- \(n, rho, sigma2, reps) {
  main.folder <- "./sims/results/"
  folder.path <- paste(main.folder, current.set(n, rho, sigma2, reps), sep = "")
  if (chck.folder(folder.path)) return()
  dir.create(folder.path)
}

subfolder.path <- \(n, rho, sigma2, reps, bbssl) {
  main.folder <- "./sims/results/"
  alg <- if (bbssl) "/bbssl" else "/mcmc"
  paste(main.folder, current.set(n, rho, sigma2, reps), alg, sep = "")
}

create.results.subfolder <- \(n, rho, sigma2, reps, bbssl = FALSE) {
  folder.path <- subfolder.path(n, rho, sigma2, reps, bbssl)
  if (chck.folder(folder.path)) return()
  dir.create(folder.path)
}

export.results.mcmc <- \(n, rho, sigma2 = 1, reps = 200) {
  create.results.folder(n, rho, sigma2, reps)
  create.results.subfolder(n, rho, sigma2, reps, bbssl = FALSE)
  out <- vector("list", 0L)
  out[["MSE"]] <- list(MSE.bnp, MSE.bayes, MSE.adapt, MSE.hs)
  out[["ACC"]] <- list(ACC.bnp, ACC.bayes, ACC.adapt, ACC.hs)
  out[["MSPE"]] <- list(MSPE.bnp, MSPE.bayes, MSPE.adapt, MSPE.hs)
  out[["TIME"]] <- list(TIME.bnp, TIME.bayes, TIME.adapt, TIME.hs)
  out[["ELPPD"]] <- list(ELPPD.bnp, ELPPD.bayes, ELPPD.adapt, ELPPD.hs)
  out[["CLUSTERS.bnp"]] <- CLUSTERS.bnp
  models <- c("bnp", "bayes", "adapt", "hs") |> paste0(".", ... = _)
  n.metrics <- length(out) - 1 # Do not count "CLUSTERS.bnp"
  for (m in seq_len(n.metrics)) names(out[[m]]) <- paste0(names(out)[m], models)
  path <- paste(subfolder.path(n, rho, sigma2, reps, bbssl = F), "/", sep = "")
  file.name <- paste(path, current.set(n, rho, sigma2, reps), sep = "")
  out.mcmc <- out
  save(out.mcmc, file = paste(file.name, ".Rdata", sep =""))
}

export.results.bbssl <- \(n, rho, sigma2 = 1, reps = 200, run = 1) {
  create.results.folder(n, rho, sigma2, reps)
  create.results.subfolder(n, rho, sigma2, reps, bbssl = TRUE)
  out <- vector("list", 0L)
  out[["MSE"]] <- list(MSE.bbssl)
  out[["ACC"]] <- list(ACC.bbssl)
  out[["MSPE"]] <- list(MSPE.bbssl)
  out[["TIME"]] <- list(TIME.bbssl)
  out[["ELPPD"]] <- list(ELPPD.bbssl)
  metrics <- names(out) |> paste0(... = _, ".bbssl")
  n.metrics <- length(metrics)
  for (m in seq_len(n.metrics)) names(out[[m]]) <- metrics[m]
  path <- paste(subfolder.path(n, rho, sigma2, reps, bbssl = T), "/", sep = "")
  file.name <- paste0(path, current.set(n, rho, sigma2, reps), "_run", run)
  out.name <- paste0("out.bbssl.", run)
  assign(out.name, out)
  save(list = out.name, file = paste(file.name, ".Rdata", sep =""))
}
