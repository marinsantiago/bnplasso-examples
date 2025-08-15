# ------------------------------------------------------------------------------
# Plot of the clustering of the regression coefficients
# ------------------------------------------------------------------------------

library(ggplot2) |> suppressPackageStartupMessages() |> suppressWarnings()
library(gridExtra) |> suppressPackageStartupMessages() |> suppressWarnings()
 
clusters_plot <- function(cls.bnp, n) {
  cls.bnp <- cls.bnp[1:25,]
  cls.dims <- dim(cls.bnp)
  xx <- seq_len(cls.dims[2])
  yy <- paste0("rep", seq_len(cls.dims[1]))
  data <- expand.grid(X = xx, Y = yy)
  data$Clusters <- as.factor(t(cls.bnp))
  xlabs <- rep("", cls.dims[2])
  xlabs[c(1, 50, 100, 150, 200)] <- as.character(c(1, 50, 100, 150, 200))
  ylabs <- rep("", 25)
  ylabs[c(1, 5, 10, 15, 20, 25)] <- as.character(c(1, 5, 10, 15, 20, 25))
  ggplot(data, aes(X, Y, fill = Clusters)) +
    geom_tile(color = "black", linewidth = 0.03) +
    ggtitle(paste0("n = ", n)) +
    #scale_fill_distiller(palette = "Blues") +
    scale_fill_manual(values = c("1" = "#08306B", "2" = "#DEEBF7")) +
    labs(x = "Coefficient index", y = "Simulation experiment") +
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank()) +
    theme(panel.background = element_rect(fill = "white")) +
    scale_y_discrete(labels = ylabs) +
    scale_x_discrete(limits = factor(1:200), labels = xlabs, breaks = xlabs) +
    theme(axis.title.x = element_text(margin = margin(t = 8))) +
    theme(axis.title.y = element_text(margin = margin(r = 8))) +
    theme(axis.text = element_text(size = 10)) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

# Clustering plot alphas -------------------------------------------------------
clusters_plot_alpha <- function(cls.bnp, alpha) {
  cls.bnp <- cls.bnp[1:25,]
  cls.dims <- dim(cls.bnp)
  xx <- seq_len(cls.dims[2])
  yy <- paste0("rep", seq_len(cls.dims[1]))
  data <- expand.grid(X = xx, Y = yy)
  data$Clusters <- as.factor(t(cls.bnp))
  xlabs <- rep("", cls.dims[2])
  xlabs[c(1, 50, 100, 150, 200)] <- as.character(c(1, 50, 100, 150, 200))
  ylabs <- rep("", 25)
  ylabs[c(1, 5, 10, 15, 20, 25)] <- as.character(c(1, 5, 10, 15, 20, 25))
  ggplot(data, aes(X, Y, fill = Clusters)) +
    geom_tile(color = "black", linewidth = 0.03) +
    ggtitle(bquote(alpha == .(alpha))) +
    #scale_fill_distiller(palette = "Blues") +
    scale_fill_manual(values = c("1" = "#08306B", "2" = "#DEEBF7")) +
    labs(x = "Coefficient index", y = "Simulation experiment") +
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank()) +
    theme(panel.background = element_rect(fill = "white")) +
    scale_y_discrete(labels = ylabs) +
    scale_x_discrete(limits = factor(1:200), labels = xlabs, breaks = xlabs) +
    theme(axis.title.x = element_text(margin = margin(t = 8))) +
    theme(axis.title.y = element_text(margin = margin(r = 8))) +
    theme(axis.text = element_text(size = 10)) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

# Clustering plot G0 -----------------------------------------------------------
clusters_plot_G0 <- function(cls.bnp, a, b) {
  cls.bnp <- cls.bnp[1:25,]
  cls.dims <- dim(cls.bnp)
  xx <- seq_len(cls.dims[2])
  yy <- paste0("rep", seq_len(cls.dims[1]))
  data <- expand.grid(X = xx, Y = yy)
  data$Clusters <- as.factor(t(cls.bnp))
  xlabs <- rep("", cls.dims[2])
  xlabs[c(1, 50, 100, 150, 200)] <- as.character(c(1, 50, 100, 150, 200))
  ylabs <- rep("", 25)
  ylabs[c(1, 5, 10, 15, 20, 25)] <- as.character(c(1, 5, 10, 15, 20, 25))
  ggplot(data, aes(X, Y, fill = Clusters)) +
    geom_tile(color = "black", linewidth = 0.03) +
    ggtitle(paste0("a = ", a, ", b = ", b, "/n")) +
    #scale_fill_distiller(palette = "Blues") +
    scale_fill_manual(values = c("1" = "#08306B", "2" = "#DEEBF7")) +
    labs(x = "Coefficient index", y = "Simulation experiment") +
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank()) +
    theme(panel.background = element_rect(fill = "white")) +
    scale_y_discrete(labels = ylabs) +
    scale_x_discrete(limits = factor(1:200), labels = xlabs, breaks = xlabs) +
    theme(axis.title.x = element_text(margin = margin(t = 8))) +
    theme(axis.title.y = element_text(margin = margin(r = 8))) +
    theme(axis.text = element_text(size = 10)) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}
