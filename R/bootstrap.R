# Parametric or residual bootstrap confidence intervals for reconstructed trajectories
#' @keywords internal
bootstrap_xHat_CI <- function(
    optimal_fit,
    prep_data,
    nbasis_val,
    npc_val,
    lambda_val,
    tau,
    basis_type,
    grid_size,
    solver,
    solver_control,
    mu_nbasis,
    B,
    alpha,
    method,
    parallel,
    n_cores,
    verbose
) {

  data_list_orig <- prep_data$data_list
  N              <- prep_data$n
  G              <- prep_data$G

  boot_solver_control <- list(maxit = 5000, trace = 0, Tolerance = 1e-04)
  boot_init <- list(U = optimal_fit$U, Xi = optimal_fit$Xi)

  fitted_list <- optimal_fit$fitted_values

  if (method == "parametric") {
    sig2 <- optimal_fit$sig2
  } else {
    residuals_list <- optimal_fit$residuals
  }

  run_one_bootstrap <- function(b) {
    boot_prep_data <- prep_data
    boot_data_list <- vector("list", N)

    for (i in seq_len(N)) {
      mi <- nrow(data_list_orig[[i]])
      if (is.null(mi) || mi == 0) {
        boot_data_list[[i]] <- data_list_orig[[i]]
        next
      }

      if (method == "parametric") {
        eps_i <- rnorm(mi, mean = 0, sd = sqrt(sig2))
        boot_centered <- fitted_list[[i]] + eps_i
      } else {
        boot_idx       <- sample.int(mi, size = mi, replace = TRUE)
        boot_residuals <- residuals_list[[i]][boot_idx]
        boot_centered  <- fitted_list[[i]] + boot_residuals
      }

      boot_data_list[[i]] <- data_list_orig[[i]]
      boot_data_list[[i]][, "centered_value"] <- boot_centered
    }

    boot_prep_data$data_list <- boot_data_list

    boot_fit <- tryCatch(
      fit_model(
        prep_data      = boot_prep_data,
        nbasis_val     = nbasis_val,
        npc_val        = npc_val,
        lambda_val     = lambda_val,
        tau            = tau,
        basis_type     = basis_type,
        grid_size      = grid_size,
        init_coef      = boot_init,
        solver         = solver,
        solver_control = boot_solver_control,
        mu_nbasis      = mu_nbasis
      ),
      error = function(e) NULL
    )

    if (is.null(boot_fit)) return(NULL)
    boot_fit$xHat
  }

  if (parallel) {
    if (is.null(n_cores)) n_cores <- max(1, detectCores() - 1)
    cl <- tryCatch(init_cluster(n_cores), error = function(e) {
      warning(sprintf("Bootstrap cluster init failed: %s. Falling back to sequential.",
                      e$message))
      NULL
    })

    if (!is.null(cl)) {
      # if (verbose) {
      #   cat(sprintf("  Running %d %s bootstrap replicates on %d workers...\n",
      #               B, method, length(cl)))
      # }

      common_vars <- c(
        "prep_data", "data_list_orig", "fitted_list", "N", "method",
        "nbasis_val", "npc_val", "lambda_val", "tau",
        "basis_type", "grid_size", "solver", "boot_solver_control", "mu_nbasis"
      )
      if (method == "parametric") {
        export_vars <- c(common_vars, "sig2")
      } else {
        export_vars <- c(common_vars, "residuals_list")
      }

      xHat_boot_list <- tryCatch({
        clusterExport(cl, varlist = export_vars, envir = environment())
        parLapply(cl, seq_len(B), run_one_bootstrap)
      }, finally = {
        safe_stop_cluster(cl)
      })
    } else {
      xHat_boot_list <- lapply(seq_len(B), function(b) {
        # if (verbose && b %% 50 == 0) {
        #   cat(sprintf("  Bootstrap replicate %d/%d\n", b, B))
        # }
        run_one_bootstrap(b)
      })
    }
  } else {
    # if (verbose) {
    #   cat(sprintf("  Running %d %s bootstrap replicates sequentially...\n", B, method))
    # }
    xHat_boot_list <- lapply(seq_len(B), function(b) {
      # if (verbose && b %% 50 == 0) {
      #   cat(sprintf("  Bootstrap replicate %d/%d\n", b, B))
      # }
      run_one_bootstrap(b)
    })
  }

  valid_boots <- Filter(Negate(is.null), xHat_boot_list)
  B_actual    <- length(valid_boots)

  if (B_actual < 2) {
    warning("Fewer than 2 bootstrap replicates succeeded; cannot compute CIs")
    return(NULL)
  }

  if (B_actual < B) {
    warning(sprintf("Only %d/%d bootstrap replicates succeeded", B_actual, B))
  }

  # if (verbose) {
  #   cat(sprintf("  Computing quantiles from %d successful replicates\n", B_actual))
  # }

  xHat_array <- array(NA, dim = c(N, G, B_actual))
  for (b in seq_len(B_actual)) {
    xHat_array[, , b] <- valid_boots[[b]]
  }

  rm(valid_boots, xHat_boot_list)

  lower_prob <- alpha / 2
  upper_prob <- 1 - alpha / 2

  xHat_CI <- vector("list", N)
  for (i in seq_len(N)) {
    subject_matrix <- xHat_array[i, , ]
    ci_matrix <- t(apply(subject_matrix, 1, quantile,
                          probs = c(lower_prob, upper_prob), na.rm = TRUE))
    colnames(ci_matrix) <- c("lower", "upper")
    xHat_CI[[i]] <- ci_matrix
  }

  xHat_CI
}
