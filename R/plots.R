# Plot AIC_K values across npc candidates
plot_AIC <- function(AIC_tab, nbasis_opt, npc_opt, lambda_opt, verbose = FALSE) {
  if (is.null(AIC_tab) || nrow(AIC_tab) == 0) {
    if (verbose) cat("No AIC table to plot\n")
    return(invisible(NULL))
  }

  varied_params <- c()
  if (length(unique(AIC_tab$nbasis)) > 1) varied_params <- c(varied_params, "nbasis")
  if (length(unique(AIC_tab$npc)) > 1) varied_params <- c(varied_params, "npc")
  if (length(unique(AIC_tab$lambda)) > 1) varied_params <- c(varied_params, "lambda")

  if (length(varied_params) == 0) {
    if (verbose) cat("No parameter variation to plot\n")
    return(invisible(NULL))
  }

  if (length(varied_params) == 1) {
    param_name <- varied_params[1]
    x_vals <- AIC_tab[[param_name]]
    y_vals <- AIC_tab$AIC_K
    order_idx <- order(x_vals)
    x_vals <- x_vals[order_idx]
    y_vals <- y_vals[order_idx]

    if (param_name == "lambda") {
      x_vals_log <- round(log10(x_vals), 3)
      x_lab <- "log10(lambda)"
    } else {
      x_vals_log <- x_vals
      x_lab <- param_name
    }

    plot(x_vals_log, y_vals, type = "o", pch = 19, lwd = 2,
         xlab = x_lab, ylab = "AIC_K",
         main = sprintf("AIC_K for %s Selection", param_name),
         xaxt = "n", cex = 1.2)
    axis(1, at = x_vals_log, labels = x_vals_log)

    opt_row <- AIC_tab[which.min(AIC_tab$AIC_K), ]
    opt_x <- if (param_name == "lambda") round(log10(opt_row[[param_name]]), 3) else opt_row[[param_name]]
    points(opt_x, opt_row$AIC_K, col = "red", pch = 19, cex = 2.5)
    grid(col = "gray80", lty = "dotted")

  } else if (length(varied_params) == 2) {
    param1 <- varied_params[1]
    param2 <- varied_params[2]
    x_vals <- sort(unique(AIC_tab[[param1]]))
    y_vals <- sort(unique(AIC_tab[[param2]]))

    aic_matrix <- matrix(NA, nrow = length(y_vals), ncol = length(x_vals))
    rownames(aic_matrix) <- y_vals
    colnames(aic_matrix) <- x_vals

    for (i in seq_along(x_vals)) {
      for (j in seq_along(y_vals)) {
        idx <- which(AIC_tab[[param1]] == x_vals[i] & AIC_tab[[param2]] == y_vals[j])
        if (length(idx) > 0) aic_matrix[j, i] <- AIC_tab$AIC_K[idx[1]]
      }
    }

    n_colors <- 10
    color_palette <- terrain.colors(n_colors)
    aic_range <- range(aic_matrix, na.rm = TRUE)
    aic_breaks <- seq(aic_range[1], aic_range[2], length.out = n_colors + 1)

    par(mar = c(5, 5, 4, 6))
    plot(NA, xlim = c(0.5, length(x_vals) + 0.5), ylim = c(0.5, length(y_vals) + 0.5),
         xlab = param1, ylab = param2,
         main = sprintf("AIC_K: %s vs %s", param1, param2),
         xaxt = "n", yaxt = "n", bty = "n")
    abline(h = 1:length(y_vals), col = "gray80", lty = "dotted")
    abline(v = 1:length(x_vals), col = "gray80", lty = "dotted")

    for (i in seq_along(x_vals)) {
      for (j in seq_along(y_vals)) {
        if (!is.na(aic_matrix[j, i])) {
          color_idx <- findInterval(aic_matrix[j, i], aic_breaks)
          color_idx <- max(1, min(color_idx, n_colors))
          rect(i - 0.4, j - 0.4, i + 0.4, j + 0.4,
               col = color_palette[color_idx], border = "black", lwd = 1.5)
          text(i, j, sprintf("%.1f", aic_matrix[j, i]), cex = 0.9, font = 2, col = "black")
        } else if (param1 == "nbasis" && param2 == "npc" && y_vals[j] > x_vals[i]) {
          rect(i - 0.4, j - 0.4, i + 0.4, j + 0.4,
               col = "gray90", border = "gray70", lwd = 1)
          text(i, j, "X", col = "red", cex = 1.5, font = 2)
        }
      }
    }

    axis(1, at = seq_along(x_vals), labels = x_vals)
    axis(2, at = seq_along(y_vals), labels = y_vals, las = 2)

    legend_image <- as.raster(matrix(rev(color_palette), ncol = 1))
    x_left <- length(x_vals) + 1.2
    x_right <- x_left + 0.3
    y_bottom <- 1
    y_top <- length(y_vals)
    rasterImage(legend_image, x_left, y_bottom, x_right, y_top)
    text(x_right + 0.2, y_bottom, sprintf("%.0f", aic_range[1]), cex = 0.8, adj = c(0, 0.5))
    text(x_right + 0.2, y_top, sprintf("%.0f", aic_range[2]), cex = 0.8, adj = c(0, 0.5))
    text(x_left + 0.15, y_top + 0.5, "AIC_K", cex = 0.9, font = 2)

    opt_row <- AIC_tab[which.min(AIC_tab$AIC_K), ]
    opt_i <- which(x_vals == opt_row[[param1]])
    opt_j <- which(y_vals == opt_row[[param2]])
    if (length(opt_i) > 0 && length(opt_j) > 0) {
      rect(opt_i - 0.4, opt_j - 0.4, opt_i + 0.4, opt_j + 0.4,
           border = "red", lwd = 3, lty = "solid")
      text(opt_i, opt_j + 0.5, "OPTIMAL", col = "red", cex = 0.8, font = 2)
    }
  } else {
    if (verbose) cat("More than 2 parameters varied. Complex plot not implemented.\n")
  }

  invisible(NULL)
}

# Plot 1-D cross-validation results (lambda or nbasis)
plot_cv_results <- function(cv_results, opt_value, log_scale = TRUE) {
  if ("lambda" %in% names(cv_results)) {
    x_raw <- cv_results$lambda
    x_vals <- if (log_scale) log10(x_raw) else x_raw
    x_opt  <- if (log_scale) log10(opt_value) else opt_value
    x_lab  <- if (log_scale) "log10(lambda)" else "lambda"
    title  <- "Cross-Validation: Lambda Selection"
  } else if ("nbasis" %in% names(cv_results)) {
    x_raw  <- cv_results$nbasis
    x_vals <- x_raw
    x_opt  <- opt_value
    x_lab  <- "nbasis"
    title  <- "Cross-Validation: Nbasis Selection"
    log_scale <- FALSE
  } else {
    stop("cv_results must contain either 'lambda' or 'nbasis' column")
  }

  y_vals <- cv_results$mean_cv_error
  se     <- cv_results$se_cv_error

  plot(x_vals, y_vals, type = "o", pch = 19, lwd = 2,
       xlab = x_lab, ylab = "Mean CV Error (MSE)", main = title)
  arrows(x_vals, y_vals - se, x_vals, y_vals + se,
         angle = 90, code = 3, length = 0.05, col = "gray50")

  opt_idx <- which(x_raw == opt_value)
  if (length(opt_idx) > 0) {
    points(x_opt, y_vals[opt_idx], col = "red", pch = 19, cex = 2)
    abline(v = x_opt, col = "red", lty = 2)
  }

  legend("topright",
         legend = sprintf("Optimal: %s", format(opt_value, scientific = TRUE, digits = 2)),
         col = "red", pch = 19, bty = "n")
  grid(col = "gray80", lty = "dotted")

  invisible(NULL)
}

# Plot 2-D grid LRPsOCV results
# Handles both "nbasis_npc" (cols: nbasis, npc) and "npc_lambda" (cols: npc, lambda).
# opt_vals: named list, e.g. list(nbasis=8, npc=3) or list(npc=2, lambda=1e-4)
plot_grid_cv_results <- function(grid_cv_results, opt_vals) {
  if (is.null(grid_cv_results) || nrow(grid_cv_results) == 0) return(invisible(NULL))

  stat_cols  <- c("mean_cv_error", "se_cv_error", "n_converged")
  param_cols <- setdiff(names(grid_cv_results), stat_cols)
  if (length(param_cols) != 2) return(invisible(NULL))

  # x-axis = param_cols[2], y-axis = param_cols[1]  (e.g. npc rows, lambda cols)
  p_row_name <- param_cols[1]
  p_col_name <- param_cols[2]

  p_row_vals <- sort(unique(grid_cv_results[[p_row_name]]))
  p_col_vals <- sort(unique(grid_cv_results[[p_col_name]]))

  cv_matrix <- matrix(NA, nrow = length(p_row_vals), ncol = length(p_col_vals))
  for (i in seq_along(p_row_vals)) {
    for (j in seq_along(p_col_vals)) {
      idx <- which(grid_cv_results[[p_row_name]] == p_row_vals[i] &
                   grid_cv_results[[p_col_name]] == p_col_vals[j])
      if (length(idx) > 0) cv_matrix[i, j] <- grid_cv_results$mean_cv_error[idx[1]]
    }
  }

  n_colors      <- 10
  color_palette <- terrain.colors(n_colors)
  cv_range      <- range(cv_matrix, na.rm = TRUE)
  if (diff(cv_range) == 0) cv_range <- cv_range + c(-0.001, 0.001)
  cv_breaks     <- seq(cv_range[1], cv_range[2], length.out = n_colors + 1)

  # Format column labels (log scale for lambda)
  col_labels <- if (p_col_name == "lambda") sprintf("%.1e", p_col_vals) else
                  as.character(p_col_vals)

  par(mar = c(5, 5, 4, 7))
  plot(NA,
       xlim = c(0.5, length(p_col_vals) + 0.5),
       ylim = c(0.5, length(p_row_vals) + 0.5),
       xlab = p_col_name, ylab = p_row_name,
       main = sprintf("Grid LRPsOCV: %s x %s", p_row_name, p_col_name),
       xaxt = "n", yaxt = "n", bty = "n")
  abline(h = seq_along(p_row_vals), col = "gray80", lty = "dotted")
  abline(v = seq_along(p_col_vals), col = "gray80", lty = "dotted")

  for (i in seq_along(p_row_vals)) {
    for (j in seq_along(p_col_vals)) {
      val <- cv_matrix[i, j]
      if (!is.na(val)) {
        color_idx <- max(1, min(findInterval(val, cv_breaks), n_colors))
        rect(j - 0.45, i - 0.45, j + 0.45, i + 0.45,
             col = color_palette[color_idx], border = "black", lwd = 1)
        text(j, i, sprintf("%.3f", val), cex = 0.7, font = 2)
      }
    }
  }

  axis(1, at = seq_along(p_col_vals), labels = col_labels,
       las = if (p_col_name == "lambda") 2 else 1, cex.axis = 0.8)
  axis(2, at = seq_along(p_row_vals), labels = p_row_vals, las = 2)

  legend_image <- as.raster(matrix(rev(color_palette), ncol = 1))
  x_left <- length(p_col_vals) + 1.1
  rasterImage(legend_image, x_left, 1, x_left + 0.3, length(p_row_vals))
  text(x_left + 0.5, 1,                  sprintf("%.3f", cv_range[1]), cex = 0.75, adj = c(0, 0.5))
  text(x_left + 0.5, length(p_row_vals), sprintf("%.3f", cv_range[2]), cex = 0.75, adj = c(0, 0.5))
  text(x_left + 0.15, length(p_row_vals) + 0.5, "CV", cex = 0.85, font = 2)

  # Highlight optimal cell
  opt_i <- which(p_row_vals == opt_vals[[p_row_name]])
  opt_j <- which(p_col_vals == opt_vals[[p_col_name]])
  if (length(opt_i) > 0 && length(opt_j) > 0) {
    rect(opt_j - 0.45, opt_i - 0.45, opt_j + 0.45, opt_i + 0.45,
         border = "red", lwd = 3)
    text(opt_j, opt_i + 0.55, "OPT", col = "red", cex = 0.8, font = 2)
  }

  invisible(NULL)
}
