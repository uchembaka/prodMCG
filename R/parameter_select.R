# Initialize a parallel cluster with the sparseFPCACG package loaded
init_cluster <- function(n_cores) {
  cl <- makeCluster(n_cores)
  tryCatch({
    clusterEvalQ(cl, { library(sparseFPCACG) })
  }, error = function(e) {
    stopCluster(cl)
    stop(sprintf("Failed to initialize cluster: %s", e$message))
  })
  cl
}

# Safely stop a parallel cluster
safe_stop_cluster <- function(cl) {
  if (!is.null(cl)) {
    tryCatch(stopCluster(cl), error = function(e) NULL)
  }
}

# Select optimal tau via eigenvalue ordering
select_tau_eigen_order <- function(
    data, nbasis, npc, lambda, tau_seq, grid_size,
    basis_type, solver, solver_control, init_coef, mu_nbasis,
    verbose = FALSE, parallel = FALSE, n_cores = NULL, cl = NULL
) {

  if (parallel) {
    if (is.null(n_cores)) n_cores <- max(1, detectCores() - 1)

    # if (verbose) {
    #   cat(sprintf("\n=== Eigen-Order Tau Selection (Parallel: %d workers) ===\n",
    #               if (!is.null(cl)) length(cl) else n_cores))
    #   cat(sprintf("nbasis=%d, npc=%d, lambda=%.2e, n=%d\n",
    #               nbasis, npc, lambda, data$n))
    #   cat(sprintf("Testing %d tau values: %s\n",
    #               length(tau_seq), paste(format(tau_seq, digits = 4), collapse = ", ")))
    #   cat("\n")
    # }

    fit_one_tau <- function(tau_val) {
      fit <- tryCatch(
        fit_model(data, nbasis, npc, lambda, tau_val, basis_type,
                  grid_size, init_coef, solver, solver_control, mu_nbasis),
        error = function(e) NULL
      )
      if (is.null(fit)) {
        return(list(tau = tau_val, AIC = NA, AIC_K = NA, fval = NA, aDF = NA,
                    convergence = NA, iter = NA, correct_order = NA, ord = NULL, fit = NULL))
      }
      list(tau = tau_val, AIC = fit$AIC, AIC_K = fit$AIC_K %||% NA,
           fval = fit$fval, aDF = fit$aDF, convergence = fit$convergence,
           iter = fit$iter,
           correct_order = identical(as.integer(fit$ord), 1:npc),
           ord = fit$ord, fit = fit)
    }

    own_cluster <- is.null(cl)
    if (own_cluster) {
      cl <- tryCatch(init_cluster(n_cores), error = function(e) {
        stop(sprintf("Failed to initialize cluster: %s", e$message))
      })
    }

    results_list <- tryCatch({
      clusterExport(cl, varlist = c(
        "data", "nbasis", "npc", "lambda", "basis_type", "grid_size",
        "init_coef", "solver", "solver_control", "mu_nbasis"
      ), envir = environment())
      parLapply(cl, tau_seq, fit_one_tau)
    }, finally = {
      if (own_cluster) safe_stop_cluster(cl)
    })

    tau_tab <- data.frame(
      tau = sapply(results_list, `[[`, "tau"),
      AIC = sapply(results_list, `[[`, "AIC"),
      AIC_K = sapply(results_list, `[[`, "AIC_K"),
      fval = sapply(results_list, `[[`, "fval"),
      aDF = sapply(results_list, `[[`, "aDF"),
      convergence = sapply(results_list, `[[`, "convergence"),
      iter = sapply(results_list, `[[`, "iter"),
      correct_order = sapply(results_list, `[[`, "correct_order")
    )

    if (verbose) {
      for (i in seq_along(tau_seq)) {
        cat(sprintf("  tau=%.4f: AIC_K=%.2f, order=%s\n",
                    tau_tab$tau[i],
                    if (!is.na(tau_tab$AIC_K[i])) tau_tab$AIC_K[i] else NA,
                    if (!is.na(tau_tab$correct_order[i]) && tau_tab$correct_order[i])
                      "correct" else "incorrect"))
      }
    }

    flag <- FALSE
    for (i in seq_along(tau_seq)) {
      co <- tau_tab$correct_order[i]
      if (!is.na(co) && co) {
        if (flag) {
          if (verbose) {
            cat(sprintf("\nSelected tau=%.4f (two consecutive correct orderings)\n",
                        tau_tab$tau[i]))
          }
          return(list(tau = tau_tab$tau[i], fit = results_list[[i]]$fit,
                      tau_results = tau_tab))
        } else {
          flag <- TRUE
        }
      } else {
        flag <- FALSE
      }
    }

    warning("No tau gave two consecutive correct orderings. Setting tau to 0.1")
    last_fit <- results_list[[length(results_list)]]$fit
    return(list(tau = 0.1, fit = last_fit, tau_results = tau_tab))

  } else {
    lambda_order <- list()
    l <- 1
    flag <- FALSE
    tau_tab <- data.frame()

    # if (verbose) {
    #   cat(sprintf("\n=== Eigen-Order Tau Selection (Sequential) ===\n"))
    #   cat(sprintf("nbasis=%d, npc=%d, lambda=%.2e, n=%d\n",
    #               nbasis, npc, lambda, data$n))
    #   cat(sprintf("Testing %d tau values\n\n", length(tau_seq)))
    # }

    for (tau in tau_seq) {
      if (verbose) cat(sprintf("  tau=%.4f ... ", tau))

      fit <- fit_model(data, nbasis, npc, lambda, tau, basis_type,
                       grid_size, init_coef, solver, solver_control, mu_nbasis)
      lambda_order[[l]] <- fit$ord
      correct_order <- identical(as.integer(lambda_order[[l]]), 1:npc)

      tau_tab <- rbind(tau_tab, data.frame(
        tau = tau,
        AIC = if (!is.null(fit$AIC)) fit$AIC else NA,
        AIC_K = if (!is.null(fit$AIC_K)) fit$AIC_K else NA,
        fval = if (!is.null(fit$fval)) fit$fval else NA,
        aDF = if (!is.null(fit$aDF)) fit$aDF else NA,
        convergence = if (!is.null(fit$convergence)) fit$convergence else NA,
        iter = if (!is.null(fit$iter)) fit$iter else NA,
        correct_order = correct_order
      ))

      if (verbose) {
        cat(sprintf("AIC_K=%.2f, order=%s\n",
                    if (!is.null(fit$AIC_K)) fit$AIC_K else NA,
                    if (correct_order) "correct" else "incorrect"))
      }

      l <- l + 1

      if (correct_order) {
        if (flag) {
          return(list(tau = tau, fit = fit, tau_results = tau_tab))
        } else {
          flag <- TRUE
        }
      } else {
        flag <- FALSE
      }

      if (tau == tau_seq[length(tau_seq)] && !correct_order) {
        warning("No tau gave correct order. Setting tau to 0.1")
        return(list(tau = 0.1, fit = fit, tau_results = tau_tab))
      }
    }

    if (!flag && tau == tau_seq[length(tau_seq)]) {
      if (identical(as.integer(lambda_order[[l - 1]]), 1:npc)) {
        warning("tau at end of vector")
        return(list(tau = tau, fit = fit, tau_results = tau_tab))
      }
    }

    warning("Unexpected end of select_tau_eigen_order")
    list(tau = tau_seq[length(tau_seq)], fit = fit, tau_results = tau_tab)
  }
}

# Select optimal tau via AIC_K (always parallel)
select_tau_aic <- function(
    data, nbasis, npc, lambda, tau_seq, grid_size,
    basis_type, solver, solver_control, init_coef, mu_nbasis,
    verbose = FALSE, n_cores = NULL, cl = NULL
) {

  # if (verbose) {
  #   cat(sprintf("\n=== AIC_K-Based Tau Selection (Parallel: %d workers) ===\n",
  #               if (!is.null(cl)) length(cl) else (n_cores %||% (detectCores() - 1))))
  #   cat(sprintf("nbasis=%d, npc=%d, lambda=%.2e, n=%d\n",
  #               nbasis, npc, lambda, data$n))
  #   cat(sprintf("Testing %d tau values: %s\n",
  #               length(tau_seq), paste(format(tau_seq, digits = 4), collapse = ", ")))
  #   cat("\n")
  # }

  fit_one_tau <- function(tau_val) {
    fit <- tryCatch(
      fit_model(data, nbasis, npc, lambda, tau_val, basis_type,
                grid_size, init_coef, solver, solver_control, mu_nbasis),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(list(tau = tau_val, AIC = NA, AIC_K = NA, fval = NA, aDF = NA,
                  convergence = NA, iter = NA, correct_order = NA, fit = NULL))
    }
    list(tau = tau_val, AIC = fit$AIC, AIC_K = fit$AIC_K %||% NA,
         fval = fit$fval, aDF = fit$aDF, convergence = fit$convergence,
         iter = fit$iter,
         correct_order = identical(as.integer(fit$ord), 1:npc), fit = fit)
  }

  own_cluster <- is.null(cl)
  if (own_cluster) {
    if (is.null(n_cores)) n_cores <- max(1, detectCores() - 1)
    cl <- tryCatch(init_cluster(n_cores), error = function(e) {
      stop(sprintf("Failed to initialize cluster: %s", e$message))
    })
  }

  results_list <- tryCatch({
    clusterExport(cl, varlist = c(
      "data", "nbasis", "npc", "lambda", "basis_type", "grid_size",
      "init_coef", "solver", "solver_control", "mu_nbasis"
    ), envir = environment())
    parLapply(cl, tau_seq, fit_one_tau)
  }, finally = {
    if (own_cluster) safe_stop_cluster(cl)
  })

  tau_tab <- data.frame(
    tau = sapply(results_list, `[[`, "tau"),
    AIC = sapply(results_list, `[[`, "AIC"),
    AIC_K = sapply(results_list, `[[`, "AIC_K"),
    fval = sapply(results_list, `[[`, "fval"),
    aDF = sapply(results_list, `[[`, "aDF"),
    convergence = sapply(results_list, `[[`, "convergence"),
    iter = sapply(results_list, `[[`, "iter"),
    correct_order = sapply(results_list, `[[`, "correct_order")
  )

  ordered_idx <- which(!is.na(tau_tab$correct_order) & tau_tab$correct_order &
                         is.finite(tau_tab$AIC_K))

  if (length(ordered_idx) > 0) {
    opt_idx  <- ordered_idx[which.min(tau_tab$AIC_K[ordered_idx])]
    tau_opt  <- tau_tab$tau[opt_idx]
    best_fit <- results_list[[opt_idx]]$fit

    # if (verbose) {
    #   cat(sprintf("\n--- Tau AIC_K Results ---\n"))
    #   cat(sprintf("Optimal tau: %.4f (AIC_K: %.2f, correct ordering)\n",
    #               tau_opt, tau_tab$AIC_K[opt_idx]))
    # }
  } else {
    warning("No tau with correct eigenvalue ordering. Falling back to eigen-order approach.")

    flag    <- FALSE
    opt_idx <- NULL
    for (i in seq_along(tau_seq)) {
      co <- tau_tab$correct_order[i]
      if (!is.na(co) && co) {
        if (flag) {
          opt_idx <- i
          break
        } else {
          flag <- TRUE
        }
      } else {
        flag <- FALSE
      }
    }

    if (!is.null(opt_idx)) {
      tau_opt  <- tau_tab$tau[opt_idx]
      best_fit <- results_list[[opt_idx]]$fit

      # if (verbose) {
      #   cat(sprintf("\n--- Tau AIC_K Results (eigen-order fallback) ---\n"))
      #   cat(sprintf("Selected tau: %.4f (two consecutive correct orderings)\n", tau_opt))
      # }
    } else {
      warning("Eigen-order fallback also failed. Setting tau to 0.1.")
      tau_opt  <- 0.1
      best_fit <- results_list[[length(results_list)]]$fit

      # if (verbose) {
      #   cat(sprintf("\n--- Tau AIC_K Results (default fallback) ---\n"))
      #   cat(sprintf("Using default tau: %.4f\n", tau_opt))
      # }
    }
  }

  tau_tab$delta_AIC_K <- tau_tab$AIC_K - min(tau_tab$AIC_K, na.rm = TRUE)

  list(tau = tau_opt, fit = best_fit, tau_results = tau_tab)
}

# Select optimal tau via k-fold cross-validation (always parallel)
select_tau_cv <- function(
    data, nbasis, npc, lambda, tau_seq, grid_size,
    basis_type, solver, solver_control, init_coef, mu_nbasis,
    nfolds = 5, seed = 123, verbose = FALSE,
    n_cores = NULL, cl = NULL
) {

  n_subjects <- data$n
  if (n_subjects < nfolds) {
    nfolds <- n_subjects
    warning(sprintf("Using %d folds (n < %d)", nfolds, nfolds))
  }

  set.seed(seed)
  fold_ids <- sample(rep(1:nfolds, length.out = n_subjects))

  # if (verbose) {
  #   cat(sprintf("\n=== %d-Fold CV for Tau Selection (Parallel: %d workers) ===\n",
  #               nfolds, if (!is.null(cl)) length(cl) else (n_cores %||% (detectCores() - 1))))
  #   cat(sprintf("nbasis=%d, npc=%d, lambda=%.2e, n=%d\n",
  #               nbasis, npc, lambda, n_subjects))
  #   cat(sprintf("Testing %d tau values: %s\n",
  #               length(tau_seq), paste(format(tau_seq, digits = 4), collapse = ", ")))
  #   cat("\n")
  # }

  sel_solver_control <- list(maxit = 300000, trace = 0, Tolerance = 1e-04)

  cv_one_tau <- function(tau_val) {
    fold_errors <- numeric(nfolds)
    for (fold in 1:nfolds) {
      test_ids <- which(fold_ids == fold)
      train_ids <- which(fold_ids != fold)
      train_data <- data
      train_data$data_list <- data$data_list[train_ids]
      train_data$n <- length(train_ids)
      train_fit <- tryCatch(
        fit_model(train_data, nbasis, npc, lambda, tau_val,
                  basis_type, grid_size, init_coef, solver,
                  sel_solver_control, mu_nbasis),
        error = function(e) NULL
      )
      if (is.null(train_fit)) { fold_errors[fold] <- NA; next }
      U_train <- train_fit$U
      basis_data <- get_basis_terms(nbasis, basis_type, grid_size, data$data_list)
      squared_errors <- c()
      for (i in test_ids) {
        y_i <- basis_data$y_vec_list[[i]]
        orthB_i <- basis_data$orthB_list[[i]]
        BU <- orthB_i %*% U_train
        xi_hat <- as.vector(solve(crossprod(BU) + tau_val * diag(npc), crossprod(BU, y_i)))
        y_pred <- orthB_i %*% U_train %*% xi_hat
        squared_errors <- c(squared_errors, (y_i - y_pred)^2)
      }
      fold_errors[fold] <- mean(squared_errors)
    }
    fold_errors
  }

  own_cluster <- is.null(cl)
  if (own_cluster) {
    if (is.null(n_cores)) n_cores <- max(1, detectCores() - 1)
    cl <- tryCatch(init_cluster(n_cores), error = function(e) {
      stop(sprintf("Failed to initialize cluster: %s", e$message))
    })
  }

  cv_errors_list <- tryCatch({
    clusterExport(cl, varlist = c(
      "data", "nbasis", "npc", "lambda", "basis_type", "grid_size",
      "init_coef", "solver", "sel_solver_control", "mu_nbasis",
      "fold_ids", "nfolds"
    ), envir = environment())
    parLapply(cl, tau_seq, cv_one_tau)
  }, finally = {
    if (own_cluster) safe_stop_cluster(cl)
  })

  cv_errors <- do.call(cbind, cv_errors_list)
  colnames(cv_errors) <- sprintf("%.4f", tau_seq)
  rownames(cv_errors) <- paste0("fold_", 1:nfolds)

  mean_cv_errors <- colMeans(cv_errors, na.rm = TRUE)
  se_cv_errors <- apply(cv_errors, 2, function(x) {
    sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  })

  tau_tab <- data.frame(
    tau = tau_seq, mean_cv_error = mean_cv_errors,
    se_cv_error = se_cv_errors,
    n_converged = apply(cv_errors, 2, function(x) sum(!is.na(x)))
  )

  valid_idx <- which(!is.na(tau_tab$mean_cv_error))
  if (length(valid_idx) == 0) stop("All tau candidates failed CV")

  opt_idx <- valid_idx[which.min(tau_tab$mean_cv_error[valid_idx])]
  tau_opt <- tau_tab$tau[opt_idx]
  min_error <- tau_tab$mean_cv_error[opt_idx]
  min_se <- tau_tab$se_cv_error[opt_idx]
  threshold <- min_error + min_se
  candidates_1se <- which(tau_tab$mean_cv_error <= threshold & tau_tab$tau <= tau_opt)
  tau_1se <- min(tau_tab$tau[candidates_1se])

  # if (verbose) {
  #   cat(sprintf("\n--- Tau CV Results ---\n"))
  #   cat(sprintf("Optimal tau:   %.4f (CV error: %.4f +/- %.4f)\n", tau_opt, min_error, min_se))
  #   cat(sprintf("1-SE rule tau: %.4f\n\n", tau_1se))
  # }

  best_fit <- fit_model(data, nbasis, npc, lambda, tau_opt, basis_type,
                        grid_size, init_coef, solver, solver_control, mu_nbasis)

  list(tau = tau_opt, tau_1se = tau_1se, fit = best_fit,
       tau_results = tau_tab, cv_errors = cv_errors)
}

# Select optimal npc via AIC_K (convergence-filtered)
select_npc_aic <- function(
    data, nbasis_val, npc_seq, lambda_val, tau, grid_size,
    basis_type, solver, solver_control, init_coef, mu_nbasis,
    verbose = FALSE, parallel = FALSE, n_cores = NULL, cl = NULL
) {

  if (length(npc_seq) < 2) stop("npc_seq must contain at least 2 values for AIC selection")

  valid_npc <- npc_seq[npc_seq <= nbasis_val]
  if (length(valid_npc) == 0) stop("All npc candidates exceed nbasis")

  # if (verbose) {
  #   cat(sprintf("\n=== AIC_K-Based NPC Selection ===\n"))
  #   cat(sprintf("nbasis=%d, lambda=%.2e, tau=%.3f, n=%d\n",
  #               nbasis_val, lambda_val, tau, data$n))
  #   cat(sprintf("Testing %d npc values: %s\n",
  #               length(valid_npc), paste(valid_npc, collapse = ", ")))
  #   if (parallel) cat(sprintf("Parallel: %d workers\n", if (!is.null(cl)) length(cl) else (n_cores %||% (detectCores() - 1))))
  #   cat("\n")
  # }
  sel_solver_control <- list(maxit = 50000, trace = 0, Tolerance = 1e-03)
  fit_one_npc <- function(np) {
    fit <- tryCatch(
      fit_model(prep_data = data, nbasis_val = nbasis_val, npc_val = np,
                lambda_val = lambda_val, tau = tau, basis_type = basis_type,
                grid_size = grid_size, init_coef = init_coef, solver = solver,
                solver_control = sel_solver_control, mu_nbasis = mu_nbasis),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(list(npc = np, AIC = NA, AIC_K = NA, fval = NA, aDF = NA,
                  convergence = NA, iter = NA, fit = NULL))
    }
    list(npc = np, AIC = fit$AIC, AIC_K = fit$AIC_K %||% NA,
         fval = fit$fval, aDF = fit$aDF, convergence = fit$convergence,
         iter = fit$iter, fit = fit)
  }

  if (parallel) {
    own_cluster <- is.null(cl)
    if (own_cluster) {
      if (is.null(n_cores)) n_cores <- max(1, detectCores() - 1)
      cl <- tryCatch(init_cluster(n_cores), error = function(e) {
        stop(sprintf("Failed to initialize cluster: %s", e$message))
      })
    }

    results_list <- tryCatch({
      clusterExport(cl, varlist = c(
        "data", "nbasis_val", "lambda_val", "tau", "basis_type", "grid_size",
        "init_coef", "solver", "solver_control", "mu_nbasis"
      ), envir = environment())
      parLapply(cl, valid_npc, fit_one_npc)
    }, finally = {
      if (own_cluster) safe_stop_cluster(cl)
    })
  } else {
    results_list <- lapply(seq_along(valid_npc), function(idx) {
      if (verbose) cat(sprintf("  Fitting npc=%d (%d/%d) ... ", valid_npc[idx], idx, length(valid_npc)))
      res <- fit_one_npc(valid_npc[idx])
      if (verbose) {
        cat(sprintf("AIC_K=%.2f, converged=%s\n",
                    if (!is.na(res$AIC_K)) res$AIC_K else NA,
                    if (!is.na(res$convergence) && res$convergence == TRUE) "yes" else "no"))
      }
      res
    })
  }

  AIC_tab <- data.frame(
    nbasis = nbasis_val,
    npc = sapply(results_list, `[[`, "npc"),
    lambda = lambda_val,
    fn = sapply(results_list, `[[`, "fval"),
    AIC = sapply(results_list, `[[`, "AIC"),
    AIC_K = sapply(results_list, `[[`, "AIC_K"),
    convergence = sapply(results_list, `[[`, "convergence"),
    aDF = sapply(results_list, `[[`, "aDF"),
    iter = sapply(results_list, `[[`, "iter")
  )

  if (nrow(AIC_tab) == 0) stop("No valid models found for npc selection")

  converged_idx <- which(!is.na(AIC_tab$convergence) & AIC_tab$convergence == TRUE)
  if (length(converged_idx) == 0) {
    warning("No converged models found for npc selection; using all finite AIC_K values")
    converged_idx <- seq_len(nrow(AIC_tab))
  }

  valid_idx <- converged_idx[is.finite(AIC_tab$AIC_K[converged_idx])]
  if (length(valid_idx) == 0) {
    valid_idx <- converged_idx[is.finite(AIC_tab$AIC[converged_idx])]
    if (length(valid_idx) == 0) stop("All AIC values are non-finite among converged models")
    opt_idx <- valid_idx[which.min(AIC_tab$AIC[valid_idx])]
  } else {
    opt_idx <- valid_idx[which.min(AIC_tab$AIC_K[valid_idx])]
  }

  npc_opt <- AIC_tab$npc[opt_idx]
  if (verbose) cat(sprintf("  Optimal npc (AIC_K, converged): %d\n", npc_opt))

  list(npc_opt = npc_opt, AIC_tab = AIC_tab[order(AIC_tab$npc), ])
}

# Select optimal lambda via AIC
select_lambda_aic <- function(
    data, nbasis_val, npc_val, lambda_seq, tau, grid_size,
    basis_type, solver, solver_control, init_coef, mu_nbasis,
    verbose = FALSE, parallel = FALSE, n_cores = NULL, cl = NULL
) {

  if (length(lambda_seq) < 2) stop("lambda_seq must contain at least 2 values for AIC selection")

  # if (verbose) {
  #   cat(sprintf("\n=== AIC-Based Lambda Selection ===\n"))
  #   cat(sprintf("nbasis=%d, npc=%d, tau=%.3f, n=%d\n", nbasis_val, npc_val, tau, data$n))
  #   cat(sprintf("Testing %d lambda values: %s\n",
  #               length(lambda_seq),
  #               paste(format(lambda_seq, scientific = TRUE, digits = 2), collapse = ", ")))
  #   if (parallel) cat(sprintf("Parallel: %d workers\n", if (!is.null(cl)) length(cl) else (n_cores %||% (detectCores() - 1))))
  #   cat("\n")
  # }

  fit_one_lambda <- function(lam) {
    fit <- tryCatch(
      fit_model(data, nbasis_val, npc_val, lam, tau,
                basis_type, grid_size, init_coef, solver, solver_control, mu_nbasis),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(list(lambda = lam, AIC = NA, AIC_K = NA, fval = NA, aDF = NA,
                  convergence = NA, iter = NA, fit = NULL))
    }
    list(lambda = lam, AIC = fit$AIC, AIC_K = fit$AIC_K %||% NA,
         fval = fit$fval, aDF = fit$aDF, convergence = fit$convergence,
         iter = fit$iter, fit = fit)
  }

  if (parallel) {
    own_cluster <- is.null(cl)
    if (own_cluster) {
      if (is.null(n_cores)) n_cores <- max(1, detectCores() - 1)
      cl <- tryCatch(init_cluster(n_cores), error = function(e) {
        stop(sprintf("Failed to initialize cluster: %s", e$message))
      })
    }

    results_list <- tryCatch({
      clusterExport(cl, varlist = c(
        "data", "nbasis_val", "npc_val", "tau", "basis_type", "grid_size",
        "init_coef", "solver", "solver_control", "mu_nbasis"
      ), envir = environment())
      parLapply(cl, lambda_seq, fit_one_lambda)
    }, finally = {
      if (own_cluster) safe_stop_cluster(cl)
    })
  } else {
    results_list <- lapply(seq_along(lambda_seq), function(idx) {
      if (verbose) cat(sprintf("  lambda=%.2e (%d/%d): ", lambda_seq[idx], idx, length(lambda_seq)))
      res <- fit_one_lambda(lambda_seq[idx])
      if (verbose) cat(sprintf("AIC=%.2f\n", if (!is.na(res$AIC)) res$AIC else NA))
      res
    })
  }

  AIC_tab <- data.frame(
    nbasis = nbasis_val, npc = npc_val,
    lambda = sapply(results_list, `[[`, "lambda"),
    fn = sapply(results_list, `[[`, "fval"),
    AIC = sapply(results_list, `[[`, "AIC"),
    AIC_K = sapply(results_list, `[[`, "AIC_K"),
    convergence = sapply(results_list, `[[`, "convergence"),
    aDF = sapply(results_list, `[[`, "aDF"),
    iter = sapply(results_list, `[[`, "iter")
  )

  valid_idx <- which(is.finite(AIC_tab$AIC))
  if (length(valid_idx) == 0) stop("All lambda candidates failed AIC-based selection")

  opt_idx <- valid_idx[which.min(AIC_tab$AIC[valid_idx])]
  lambda_opt <- AIC_tab$lambda[opt_idx]
  AIC_tab$delta_AIC <- AIC_tab$AIC - min(AIC_tab$AIC, na.rm = TRUE)

  # if (verbose) {
  #   cat(sprintf("\n--- Lambda AIC Results ---\n"))
  #   cat(sprintf("Optimal lambda: %.2e (AIC: %.2f)\n", lambda_opt, AIC_tab$AIC[opt_idx]))
  # }

  list(lambda_opt = lambda_opt, AIC_tab = AIC_tab[order(AIC_tab$lambda), ])
}

# Select optimal nbasis via AIC
select_nbasis_aic <- function(
    data, nbasis_seq, npc_val, lambda_val, tau, grid_size,
    basis_type, solver, solver_control, init_coef, mu_nbasis,
    verbose = FALSE, parallel = FALSE, n_cores = NULL, cl = NULL
) {

  if (length(nbasis_seq) < 2) stop("nbasis_seq must contain at least 2 values for AIC selection")

  valid_nbasis <- nbasis_seq[nbasis_seq >= npc_val]
  if (length(valid_nbasis) == 0) stop("All nbasis candidates are smaller than npc")

  # if (verbose) {
  #   cat(sprintf("\n=== AIC-Based Nbasis Selection ===\n"))
  #   cat(sprintf("npc=%d, lambda=%.2e, tau=%.3f, n=%d\n", npc_val, lambda_val, tau, data$n))
  #   cat(sprintf("Testing %d nbasis values: %s\n",
  #               length(valid_nbasis), paste(valid_nbasis, collapse = ", ")))
  #   if (parallel) cat(sprintf("Parallel: %d workers\n", if (!is.null(cl)) length(cl) else (n_cores %||% (detectCores() - 1))))
  #   cat("\n")
  # }

  fit_one_nbasis <- function(nb) {
    fit <- tryCatch(
      fit_model(data, nb, npc_val, lambda_val, tau,
                basis_type, grid_size, init_coef, solver, solver_control, mu_nbasis),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(list(nbasis = nb, AIC = NA, AIC_K = NA, fval = NA, aDF = NA,
                  convergence = NA, iter = NA, fit = NULL))
    }
    list(nbasis = nb, AIC = fit$AIC, AIC_K = fit$AIC_K %||% NA,
         fval = fit$fval, aDF = fit$aDF, convergence = fit$convergence,
         iter = fit$iter, fit = fit)
  }

  if (parallel) {
    own_cluster <- is.null(cl)
    if (own_cluster) {
      if (is.null(n_cores)) n_cores <- max(1, detectCores() - 1)
      cl <- tryCatch(init_cluster(n_cores), error = function(e) {
        stop(sprintf("Failed to initialize cluster: %s", e$message))
      })
    }

    results_list <- tryCatch({
      clusterExport(cl, varlist = c(
        "data", "npc_val", "lambda_val", "tau", "basis_type", "grid_size",
        "init_coef", "solver", "solver_control", "mu_nbasis"
      ), envir = environment())
      parLapply(cl, valid_nbasis, fit_one_nbasis)
    }, finally = {
      if (own_cluster) safe_stop_cluster(cl)
    })
  } else {
    results_list <- lapply(seq_along(valid_nbasis), function(idx) {
      if (verbose) cat(sprintf("  nbasis=%d (%d/%d): ", valid_nbasis[idx], idx, length(valid_nbasis)))
      res <- fit_one_nbasis(valid_nbasis[idx])
      if (verbose) cat(sprintf("AIC=%.2f\n", if (!is.na(res$AIC)) res$AIC else NA))
      res
    })
  }

  AIC_tab <- data.frame(
    nbasis = sapply(results_list, `[[`, "nbasis"),
    npc = npc_val, lambda = lambda_val,
    fn = sapply(results_list, `[[`, "fval"),
    AIC = sapply(results_list, `[[`, "AIC"),
    AIC_K = sapply(results_list, `[[`, "AIC_K"),
    convergence = sapply(results_list, `[[`, "convergence"),
    aDF = sapply(results_list, `[[`, "aDF"),
    iter = sapply(results_list, `[[`, "iter")
  )

  valid_idx <- which(is.finite(AIC_tab$AIC))
  if (length(valid_idx) == 0) stop("All nbasis candidates failed AIC-based selection")

  opt_idx <- valid_idx[which.min(AIC_tab$AIC[valid_idx])]
  nbasis_opt <- AIC_tab$nbasis[opt_idx]
  AIC_tab$delta_AIC <- AIC_tab$AIC - min(AIC_tab$AIC, na.rm = TRUE)

  # if (verbose) {
  #   cat(sprintf("\n--- Nbasis AIC Results ---\n"))
  #   cat(sprintf("Optimal nbasis: %d (AIC: %.2f)\n", nbasis_opt, AIC_tab$AIC[opt_idx]))
  # }

  list(nbasis_opt = nbasis_opt, AIC_tab = AIC_tab[order(AIC_tab$nbasis), ])
}

# Select optimal lambda via k-fold cross-validation
select_lambda_cv <- function(
    data, nbasis_val, npc_val, lambda_seq, tau, grid_size,
    basis_type, solver, solver_control, init_coef, mu_nbasis,
    nfolds = 5, seed = 123, verbose = TRUE,
    parallel = FALSE, n_cores = NULL, cl = NULL
) {

  n_subjects <- data$n
  if (n_subjects < nfolds) {
    nfolds <- n_subjects
    warning(sprintf("Using %d folds (n < %d)", nfolds, nfolds))
  }

  set.seed(seed)
  fold_ids <- sample(rep(1:nfolds, length.out = n_subjects))

  # if (verbose) {
  #   cat(sprintf("\n=== %d-Fold CV for Lambda Selection ===\n", nfolds))
  #   cat(sprintf("nbasis=%d, npc=%d, tau=%.3f, n=%d\n", nbasis_val, npc_val, tau, n_subjects))
  #   cat(sprintf("Testing %d lambda values\n", length(lambda_seq)))
  #   if (parallel) cat(sprintf("Parallel: %d workers\n", if (!is.null(cl)) length(cl) else (n_cores %||% (detectCores() - 1))))
  #   cat("\n")
  # }

  sel_solver_control <- list(maxit = 50000, trace = 0, Tolerance = 1e-03)

  cv_one_lambda <- function(lam) {
    fold_errors <- numeric(nfolds)
    for (fold in 1:nfolds) {
      test_ids <- which(fold_ids == fold)
      train_ids <- which(fold_ids != fold)
      train_data <- data
      train_data$data_list <- data$data_list[train_ids]
      train_data$n <- length(train_ids)
      if (!is.null(data$tau_init)) {
        train_data$tau_init <- list(U = data$tau_init$U,
                                    Xi = data$tau_init$Xi[train_ids, , drop = FALSE])
      }
      train_fit <- tryCatch(
        fit_model(train_data, nbasis_val, npc_val, lam, tau,
                  basis_type, grid_size, init_coef, solver,
                  sel_solver_control, mu_nbasis),
        error = function(e) NULL
      )
      if (is.null(train_fit)) { fold_errors[fold] <- NA; next }
      U_train <- train_fit$U
      basis_data <- get_basis_terms(nbasis_val, basis_type, grid_size, data$data_list)
      squared_errors <- c()
      for (i in test_ids) {
        y_i <- basis_data$y_vec_list[[i]]
        orthB_i <- basis_data$orthB_list[[i]]
        BU <- orthB_i %*% U_train
        xi_hat <- as.vector(solve(crossprod(BU) + tau * diag(npc_val), crossprod(BU, y_i)))
        y_pred <- orthB_i %*% U_train %*% xi_hat
        squared_errors <- c(squared_errors, (y_i - y_pred)^2)
      }
      fold_errors[fold] <- mean(squared_errors)
    }
    fold_errors
  }

  if (parallel) {
    own_cluster <- is.null(cl)
    if (own_cluster) {
      if (is.null(n_cores)) n_cores <- max(1, detectCores() - 1)
      cl <- tryCatch(init_cluster(n_cores), error = function(e) {
        stop(sprintf("Failed to initialize cluster: %s", e$message))
      })
    }

    cv_errors_list <- tryCatch({
      clusterExport(cl, varlist = c(
        "data", "nbasis_val", "npc_val", "tau", "basis_type", "grid_size",
        "init_coef", "solver", "sel_solver_control", "mu_nbasis",
        "fold_ids", "nfolds"
      ), envir = environment())
      parLapply(cl, lambda_seq, cv_one_lambda)
    }, finally = {
      if (own_cluster) safe_stop_cluster(cl)
    })
  } else {
    cv_errors_list <- lapply(seq_along(lambda_seq), function(lam_idx) {
      lam <- lambda_seq[lam_idx]
      # if (verbose) cat(sprintf("Lambda %.2e (%d/%d): ", lam, lam_idx, length(lambda_seq)))
      errs <- cv_one_lambda(lam)
      # if (verbose) cat(sprintf("MSE = %.4f\n", mean(errs, na.rm = TRUE)))
      errs
    })
  }

  cv_errors <- do.call(cbind, cv_errors_list)
  colnames(cv_errors) <- sprintf("%.2e", lambda_seq)
  rownames(cv_errors) <- paste0("fold_", 1:nfolds)

  mean_cv_errors <- colMeans(cv_errors, na.rm = TRUE)
  se_cv_errors <- apply(cv_errors, 2, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

  cv_results <- data.frame(
    lambda = lambda_seq, mean_cv_error = mean_cv_errors,
    se_cv_error = se_cv_errors,
    n_converged = apply(cv_errors, 2, function(x) sum(!is.na(x)))
  )

  valid_idx <- which(!is.na(cv_results$mean_cv_error))
  if (length(valid_idx) == 0) stop("All lambdas failed CV")

  opt_idx <- valid_idx[which.min(cv_results$mean_cv_error[valid_idx])]
  lambda_opt <- cv_results$lambda[opt_idx]
  min_error <- cv_results$mean_cv_error[opt_idx]
  min_se <- cv_results$se_cv_error[opt_idx]
  threshold <- min_error + min_se
  candidates_1se <- which(cv_results$mean_cv_error <= threshold &
                            cv_results$lambda >= lambda_opt)
  lambda_1se <- max(cv_results$lambda[candidates_1se])

  # if (verbose) {
  #   cat(sprintf("\n--- CV Results ---\n"))
  #   cat(sprintf("Optimal lambda:   %.2e (CV error: %.4f +/- %.4f)\n", lambda_opt, min_error, min_se))
  #   cat(sprintf("1-SE rule lambda: %.2e\n\n", lambda_1se))
  # }

  list(lambda_opt = lambda_opt, lambda_1se = lambda_1se,
       cv_results = cv_results, cv_errors = cv_errors)
}

# Select optimal nbasis via k-fold cross-validation
select_nbasis_cv <- function(
    data, nbasis_seq, npc_val, lambda_val, tau, grid_size,
    basis_type, solver, solver_control, init_coef, mu_nbasis,
    nfolds = 5, seed = 123, verbose = TRUE,
    parallel = FALSE, n_cores = NULL, cl = NULL
) {

  n_subjects <- data$n
  if (n_subjects < nfolds) {
    nfolds <- n_subjects
    warning(sprintf("Using %d folds (n < %d)", nfolds, nfolds))
  }

  set.seed(seed)
  fold_ids <- sample(rep(1:nfolds, length.out = n_subjects))

  # if (verbose) {
  #   cat(sprintf("\n=== %d-Fold CV for Nbasis Selection ===\n", nfolds))
  #   cat(sprintf("npc=%d, lambda=%.2e, tau=%.3f, n=%d\n", npc_val, lambda_val, tau, n_subjects))
  #   cat(sprintf("Testing %d nbasis values: %s\n",
  #               length(nbasis_seq), paste(nbasis_seq, collapse = ", ")))
  #   if (parallel) cat(sprintf("Parallel: %d workers\n", if (!is.null(cl)) length(cl) else (n_cores %||% (detectCores() - 1))))
  #   cat("\n")
  # }

  sel_solver_control <- list(maxit = 300000, trace = 0, Tolerance = 1e-03)

  cv_one_nbasis <- function(nb) {
    fold_errors <- numeric(nfolds)
    for (fold in 1:nfolds) {
      test_ids <- which(fold_ids == fold)
      train_ids <- which(fold_ids != fold)
      train_data <- data
      train_data$data_list <- data$data_list[train_ids]
      train_data$n <- length(train_ids)
      if (!is.null(data$tau_init)) {
        train_data$tau_init <- list(U = data$tau_init$U,
                                    Xi = data$tau_init$Xi[train_ids, , drop = FALSE])
      }
      train_fit <- tryCatch(
        fit_model(train_data, nb, npc_val, lambda_val, tau,
                  basis_type, grid_size, init_coef, solver,
                  sel_solver_control, mu_nbasis),
        error = function(e) NULL
      )
      if (is.null(train_fit)) { fold_errors[fold] <- NA; next }
      U_train <- train_fit$U
      basis_data <- get_basis_terms(nb, basis_type, grid_size, data$data_list)
      squared_errors <- c()
      for (i in test_ids) {
        y_i <- basis_data$y_vec_list[[i]]
        orthB_i <- basis_data$orthB_list[[i]]
        BU <- orthB_i %*% U_train
        xi_hat <- as.vector(solve(crossprod(BU) + tau * diag(npc_val), crossprod(BU, y_i)))
        y_pred <- orthB_i %*% U_train %*% xi_hat
        squared_errors <- c(squared_errors, (y_i - y_pred)^2)
      }
      fold_errors[fold] <- mean(squared_errors)
    }
    fold_errors
  }

  if (parallel) {
    own_cluster <- is.null(cl)
    if (own_cluster) {
      if (is.null(n_cores)) n_cores <- max(1, detectCores() - 1)
      cl <- tryCatch(init_cluster(n_cores), error = function(e) {
        stop(sprintf("Failed to initialize cluster: %s", e$message))
      })
    }

    cv_errors_list <- tryCatch({
      clusterExport(cl, varlist = c(
        "data", "npc_val", "lambda_val", "tau", "basis_type", "grid_size",
        "init_coef", "solver", "sel_solver_control", "mu_nbasis",
        "fold_ids", "nfolds"
      ), envir = environment())
      parLapply(cl, nbasis_seq, cv_one_nbasis)
    }, finally = {
      if (own_cluster) safe_stop_cluster(cl)
    })
  } else {
    cv_errors_list <- lapply(seq_along(nbasis_seq), function(nb_idx) {
      nb <- nbasis_seq[nb_idx]
      if (verbose) cat(sprintf("Nbasis %d (%d/%d): ", nb, nb_idx, length(nbasis_seq)))
      errs <- cv_one_nbasis(nb)
      if (verbose) cat(sprintf("MSE = %.4f\n", mean(errs, na.rm = TRUE)))
      errs
    })
  }

  cv_errors <- do.call(cbind, cv_errors_list)
  colnames(cv_errors) <- as.character(nbasis_seq)
  rownames(cv_errors) <- paste0("fold_", 1:nfolds)

  mean_cv_errors <- colMeans(cv_errors, na.rm = TRUE)
  se_cv_errors <- apply(cv_errors, 2, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

  cv_results <- data.frame(
    nbasis = nbasis_seq, mean_cv_error = mean_cv_errors,
    se_cv_error = se_cv_errors,
    n_converged = apply(cv_errors, 2, function(x) sum(!is.na(x)))
  )

  valid_idx <- which(!is.na(cv_results$mean_cv_error))
  if (length(valid_idx) == 0) stop("All nbasis values failed CV")

  opt_idx <- valid_idx[which.min(cv_results$mean_cv_error[valid_idx])]
  nbasis_opt <- cv_results$nbasis[opt_idx]
  min_error <- cv_results$mean_cv_error[opt_idx]
  min_se <- cv_results$se_cv_error[opt_idx]
  candidates_1se <- which(cv_results$mean_cv_error <= min_error + min_se &
                            cv_results$nbasis <= nbasis_opt)
  nbasis_1se <- min(cv_results$nbasis[candidates_1se])

  # if (verbose) {
  #   cat(sprintf("\n--- CV Results ---\n"))
  #   cat(sprintf("Optimal nbasis:   %d (CV error: %.4f +/- %.4f)\n", nbasis_opt, min_error, min_se))
  #   cat(sprintf("1-SE rule nbasis: %d\n\n", nbasis_1se))
  # }

  list(nbasis_opt = nbasis_opt, nbasis_1se = nbasis_1se,
       cv_results = cv_results, cv_errors = cv_errors)
}
