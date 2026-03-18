# Preprocess raw data for sparse FPCA
preprocess_data <- function(data, basis_type = "bspline", G = 101, obs_range = c(0, 1), mu_nbasis = 15) {
  if (!is.matrix(data) && !is.data.frame(data)) stop("data must be a matrix or data frame")
  if (ncol(data) != 3) stop("data must have exactly 3 columns: ID, time, value")

  aT <- obs_range[1]
  bT <- obs_range[2]
  data <- as.matrix(data)
  colnames(data) <- c("ID", "time", "value")
  data[, "time"] <- round(data[, "time"], 4)
  data <- data[order(data[, "ID"], data[, "time"]), ]
  data[, "time"] <- (data[, "time"] - aT) / (bT - aT)

  subject_ids <- unique(data[, "ID"])
  N <- length(subject_ids)
  grid <- seq(0, 1, length.out = G)

  mean_result <- get_mean(data = data, estGrid = grid, type = basis_type, knots = mu_nbasis, lambdaVec = NULL)
  mu_fits <- mean_result$fits
  mu_grid <- mean_result$mu
  mu_fxn <- mean_result$fdobj$fd

  data_centered <- data
  data_centered[, "value"] <- data[, "value"] - mu_fits

  data_list <- vector("list", N)
  for (i in seq_len(N)) {
    idx <- data_centered[, "ID"] == subject_ids[i]
    subject_times <- data_centered[idx, "time"]
    subject_values <- data_centered[idx, "value"]
    binned_values <- bin_to_grid_na(times = subject_times, values = subject_values, grid = grid, method = "nearest")
    non_na <- which(!is.na(binned_values))
    if (length(non_na) > 0) {
      data_list[[i]] <- cbind(time = grid[non_na], centered_value = binned_values[non_na])
    } else {
      data_list[[i]] <- matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("time", "centered_value")))
    }
  }

  binned_centered_data <- do.call(rbind, lapply(seq_len(N), function(i) {
    if (nrow(data_list[[i]]) == 0) return(NULL)
    cbind(ID = subject_ids[i], data_list[[i]])
  }))
  if (is.null(binned_centered_data)) {
    binned_centered_data <- matrix(0, nrow = 0, ncol = 3)
    colnames(binned_centered_data) <- c("ID", "time", "centered_value")
  }

  list(
    centered_data = binned_centered_data, data_list = data_list,
    mu_fxn = mu_fxn, grid = grid, subject_ids = subject_ids,
    n = N, G = G, basis_type = basis_type,
    mean_estimation = list(mu_grid = mu_grid, basis = mean_result$basis,
                           fdobj = mean_result$fdobj, lambda = mean_result$lambda)
  )
}

# Construct orthogonalized basis terms for all subjects
get_basis_terms <- function(nbasis, basis_type, G, data_list) {
  grid <- seq(0, 1, length.out = G)
  N <- length(data_list)

  if (basis_type == "bspline") {
    basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis, norder = 4)
  } else {
    basis <- create.fourier.basis(rangeval = c(0, 1), nbasis = nbasis)
  }

  B_grid <- eval.basis(grid, basis)
  int_BTB <- inprod(basis, basis)
  Gram_mat <- solve(chol(int_BTB))

  orthB_list <- vector("list", length = N)
  y_vec_list <- vector("list", length = N)

  for (i in 1:N) {
    sub_i <- data_list[[i]]
    if (nrow(sub_i) > 0) {
      Bi <- eval.basis(sub_i[, "time"], basis)
      orthB_list[[i]] <- Bi %*% Gram_mat
      y_vec_list[[i]] <- sub_i[, "centered_value"]
    } else {
      orthB_list[[i]] <- matrix(0, nrow = 0, ncol = nbasis)
      y_vec_list[[i]] <- numeric(0)
    }
  }

  if (basis_type != "bspline") {
    P <- fourierpen(basis, Lfdobj = 2)
  } else {
    P <- bsplinepen(basis, Lfdobj = 2)
  }

  orthP <- t(Gram_mat) %*% P %*% Gram_mat
  orthB_grid <- B_grid %*% Gram_mat

  list(basis = basis, B_grid = B_grid, int_BTB = int_BTB, Gram_mat = Gram_mat,
       orthB_list = orthB_list, y_vec_list = y_vec_list, P = P,
       orthP = orthP, orthB_grid = orthB_grid, grid = grid,
       nbasis = nbasis, basis_type = basis_type)
}

# Bin irregularly-spaced observations to a regular grid
bin_to_grid_na <- function(times, values, grid, method = "nearest", bandwidth = NULL) {
  G <- length(grid)
  binned <- rep(NA, G)
  if (method == "nearest") {
    for (j in seq_along(times)) {
      grid_idx <- which.min(abs(grid - times[j]))
      if (is.na(binned[grid_idx])) {
        binned[grid_idx] <- values[j]
      } else {
        binned[grid_idx] <- mean(c(binned[grid_idx], values[j]), na.rm = TRUE)
      }
    }
  } else if (method == "average") {
    if (is.null(bandwidth)) bandwidth <- 2 * (grid[2] - grid[1])
    for (g in 1:G) {
      distances <- abs(times - grid[g])
      weights <- exp(-distances^2 / (2 * bandwidth^2))
      if (any(weights > 0.01)) {
        weighted_sum <- sum(values * weights, na.rm = TRUE)
        total_weight <- sum(weights, na.rm = TRUE)
        if (total_weight > 0) binned[g] <- weighted_sum / total_weight
      }
    }
  }
  binned
}

# Estimate the mean function via penalized smoothing with GCV
get_mean <- function(data, estGrid, type = "bspline", knots = NULL, lambdaVec = NULL) {
  if (is.null(knots)) knots <- length(estGrid)
  if (is.null(lambdaVec)) lambdaVec <- seq(1e-4, 1, length.out = 20)

  grid <- estGrid
  raw_data <- data[, 3]
  wk_ind <- data[, 2]

  if (type == "bspline") {
    basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = knots)
  } else {
    basis <- create.fourier.basis(rangeval = c(0, 1), nbasis = knots)
  }

  gcv_vec <- numeric(length(lambdaVec))
  fdPar_list <- vector("list", length(lambdaVec))

  for (i in seq_along(lambdaVec)) {
    fdParobj <- fdPar(fdobj = basis, Lfdobj = 2, lambda = lambdaVec[i])
    smooth <- smooth.basis(argvals = wk_ind, y = raw_data, fdParobj = fdParobj)
    gcv_vec[i] <- smooth$gcv
    fdPar_list[[i]] <- fdParobj
  }

  best_idx <- which.min(gcv_vec)
  lam <- lambdaVec[best_idx]
  fdParobj <- fdPar_list[[best_idx]]
  fdobj <- smooth.basis(argvals = wk_ind, y = raw_data, fdParobj = fdParobj)
  mu <- eval.basis(grid, basis) %*% fdobj$fd$coefs
  fits <- eval.basis(wk_ind, basis) %*% fdobj$fd$coefs

  list(mu = mu, fits = fits, basis = basis, fdobj = fdobj, lambda = lam)
}

# Order eigenvectors by decreasing eigenvalue
eig_order <- function(U, Xi) {
  eigenvalues <- apply(Xi, 2, var)
  ord <- order(eigenvalues, decreasing = TRUE)
  eigenvalues <- eigenvalues[ord]
  U <- U[, ord, drop = FALSE]
  Xi <- Xi[, ord, drop = FALSE]
  list(U = U, Xi = Xi, lambda = eigenvalues, ord = ord)
}

# Reset eigenfunctions and mean to original observation range
obs_range_reset <- function(result, basis_type, obs_range, G, mu_nbasis) {
  K <- ncol(result$Xi)
  D <- nrow(result$U)
  grid <- seq(obs_range[1], obs_range[2], length = G)

  if (basis_type == "bspline") {
    mu_rng_basis <- create.bspline.basis(obs_range, mu_nbasis, norder = 4)
    eig_rng_basis <- create.bspline.basis(obs_range, D, norder = 4)
  } else {
    mu_rng_basis <- create.fourier.basis(obs_range, mu_nbasis)
    eig_rng_basis <- create.fourier.basis(obs_range, D)
  }

  eigenfunctions <- lapply(1:K, function(k) fd(result$U_original[, k], eig_rng_basis))
  mu_fd <- fd(result$data$prep_data$mean_estimation$fdobj$fd$coefs, mu_rng_basis)

  list(grid = grid, mu_fd = mu_fd, eigenfunctions = eigenfunctions)
}
