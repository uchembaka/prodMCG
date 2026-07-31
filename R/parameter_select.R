# Initialize a parallel cluster with prodMCG loaded
init_cluster <- function(n_cores) {
  cl <- makeCluster(n_cores)
  tryCatch({
    clusterEvalQ(cl, { library(prodMCG) })
  }, error = function(e) {
    stopCluster(cl)
    stop(sprintf("Failed to initialize cluster: %s", e$message))
  })
  cl
}

# Safely stop a parallel cluster
safe_stop_cluster <- function(cl) {
  if (!is.null(cl)) tryCatch(stopCluster(cl), error = function(e) NULL)
}


# Select optimal tau via eigenvalue ordering + convergence
select_tau_eigen_order <- function(
    data, nbasis, npc, lambda, tau_seq, grid_size,
    basis_type, solver, solver_control, init_coef, mu_nbasis,
    verbose = FALSE, parallel = FALSE, n_cores = NULL, cl = NULL,
    scores_weights = NULL, seed = NULL
) {
  if (parallel) {
    if (is.null(n_cores)) n_cores <- max(1, detectCores() - 1)

    # fit_one_tau <- function(tau_val) {
    #   fit <- tryCatch(
    #     fit_model(data, nbasis, npc, lambda, tau_val, basis_type,
    #               grid_size, init_coef, solver, solver_control, mu_nbasis,
    #               scores_weights = scores_weights),
    #     error = function(e) NULL
    #   )
    #   if (is.null(fit)) {
    #     return(list(tau = tau_val, fval = NA, convergence = NA,
    #                 iter = NA, correct_order = NA, ord = NULL, fit = NULL))
    #   }
    #   list(tau = tau_val, fval = fit$fval, convergence = fit$convergence,
    #        iter = fit$iter,
    #        correct_order = identical(as.integer(fit$ord), 1:npc),
    #        ord = fit$ord, fit = fit)
    # }
    fit_one_tau <- function(tau_idx) {
      tau_val <- tau_seq[tau_idx]

      if (!is.null(seed)) {
        set.seed(seed + 10000L + tau_idx)
      }

      fit <- tryCatch(
        fit_model(data, nbasis, npc, lambda, tau_val, basis_type,
                  grid_size, init_coef, solver, solver_control, mu_nbasis,
                  scores_weights = scores_weights),
        error = function(e) NULL
      )

      if (is.null(fit)) {
        return(list(tau = tau_val, fval = NA, convergence = NA,
                    iter = NA, correct_order = NA, ord = NULL, fit = NULL))
      }

      list(tau = tau_val, fval = fit$fval, convergence = fit$convergence,
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
        "init_coef", "solver", "solver_control", "mu_nbasis", "scores_weights"
      ), envir = environment())
      # parLapply(cl, tau_seq, fit_one_tau)
      parLapply(cl, seq_along(tau_seq), fit_one_tau)
    }, finally = {
      if (own_cluster) safe_stop_cluster(cl)
    })

    tau_tab <- data.frame(
      tau           = sapply(results_list, `[[`, "tau"),
      fval          = sapply(results_list, `[[`, "fval"),
      convergence   = sapply(results_list, `[[`, "convergence"),
      iter          = sapply(results_list, `[[`, "iter"),
      correct_order = sapply(results_list, `[[`, "correct_order")
    )
    selected_fit_fn <- function(idx) results_list[[idx]]$fit

  } else {
    fits_list <- vector("list", length(tau_seq))
    tau_tab   <- data.frame()

    # for (i in seq_along(tau_seq)) {
    #   tau_val <- tau_seq[i]
    #   fit <- tryCatch(
    #     fit_model(data, nbasis, npc, lambda, tau_val, basis_type,
    #               grid_size, init_coef, solver, solver_control, mu_nbasis,
    #               scores_weights = scores_weights),
    #     error = function(e) NULL
    #   )
    #   fits_list[[i]] <- fit
    #   tau_tab <- rbind(tau_tab, data.frame(
    #     tau           = tau_val,
    #     fval          = if (!is.null(fit)) fit$fval          %||% NA else NA,
    #     convergence   = if (!is.null(fit)) fit$convergence   %||% NA else NA,
    #     iter          = if (!is.null(fit)) fit$iter          %||% NA else NA,
    #     correct_order = if (!is.null(fit))
    #                       identical(as.integer(fit$ord), 1:npc) else NA
    #   ))
    # }
    for (i in seq_along(tau_seq)) {
      fit <- fit_one_tau(i)
      fits_list[[i]] <- fit
      tau_tab <- rbind(tau_tab, data.frame(
        tau           = tau_val,
        fval          = if (!is.null(fit)) fit$fval          %||% NA else NA,
        convergence   = if (!is.null(fit)) fit$convergence   %||% NA else NA,
        iter          = if (!is.null(fit)) fit$iter          %||% NA else NA,
        correct_order = if (!is.null(fit))
                          identical(as.integer(fit$ord), 1:npc) else NA
      ))
    }
    selected_fit_fn <- function(idx) fits_list[[idx]]$fit
  }

  # First tau after the last incorrect-ordering index
  correct <- tau_tab$correct_order & tau_tab$convergence
  correct[is.na(correct)] <- FALSE
  incorrect_idx <- which(!correct)

  if (length(incorrect_idx) == 0) {
    selected_idx <- 1L
  } else {
    last_incorrect <- max(incorrect_idx)
    correct_after  <- which(correct & seq_along(correct) > last_incorrect)
    if (length(correct_after) == 0L) {
      warning("No correct ordering after last incorrect; falling back to tau = 0.1")
      return(list(tau = 0.1, fit = NULL, tau_results = tau_tab))
    }
    selected_idx <- correct_after[1]
  }

  list(tau         = tau_tab$tau[selected_idx],
       fit         = selected_fit_fn(selected_idx),
       tau_results = tau_tab)
}


# Select optimal npc via AIC_K (convergence-filtered)
select_npc_aic <- function(
    data, nbasis_val, npc_seq, lambda_val, tau, grid_size,
    basis_type, solver, solver_control, init_coef, mu_nbasis,
    verbose = FALSE, parallel = FALSE, n_cores = NULL, cl = NULL,
    scores_weights = NULL, seed = NULL
) {
  if (length(npc_seq) < 2) stop("npc_seq must contain at least 2 values for AIC selection")

  valid_npc <- npc_seq[npc_seq <= nbasis_val]
  if (length(valid_npc) == 0) stop("All npc candidates exceed nbasis")

  # fit_one_npc <- function(np) {
  #   fit <- tryCatch(
  #     fit_model(prep_data = data, nbasis_val = nbasis_val, npc_val = np,
  #               lambda_val = lambda_val, tau = tau, basis_type = basis_type,
  #               grid_size = grid_size, init_coef = init_coef, solver = solver,
  #               solver_control = solver_control, mu_nbasis = mu_nbasis,
  #               scores_weights = scores_weights),
  #     error = function(e) NULL
  #   )
  #   if (is.null(fit)) return(list(npc = np, AIC_K = NA, fval = NA,
  #                                 convergence = NA, iter = NA))
  #   list(npc = np, AIC_K = fit$AIC_K %||% NA, fval = fit$fval,
  #        convergence = fit$convergence, iter = fit$iter)
  # }

  fit_one_npc <- function(idx) {
    np <- valid_npc[idx]

    if (!is.null(seed)) {
      set.seed(seed + 20000L + idx)
    }

    fit <- tryCatch(
      fit_model(prep_data = data, nbasis_val = nbasis_val, npc_val = np,
                lambda_val = lambda_val, tau = tau, basis_type = basis_type,
                grid_size = grid_size, init_coef = init_coef, solver = solver,
                solver_control = solver_control, mu_nbasis = mu_nbasis,
                scores_weights = scores_weights),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      return(list(npc = np, AIC_K = NA, fval = NA,
                  convergence = NA, iter = NA))
    }

    list(npc = np, AIC_K = fit$AIC_K %||% NA, fval = fit$fval,
         convergence = fit$convergence, iter = fit$iter)
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
        "init_coef", "solver", "solver_control", "mu_nbasis", "scores_weights"
      ), envir = environment())
      # parLapply(cl, valid_npc, fit_one_npc)
      parLapply(cl, seq_along(valid_npc), fit_one_npc)
    }, finally = {
      if (own_cluster) safe_stop_cluster(cl)
    })
  } else {
    # results_list <- lapply(valid_npc, fit_one_npc)
    results_list <- lapply(seq_along(valid_npc), fit_one_npc)
  }

  AIC_tab <- data.frame(
    nbasis      = nbasis_val,
    npc         = sapply(results_list, `[[`, "npc"),
    lambda      = lambda_val,
    fn          = sapply(results_list, `[[`, "fval"),
    AIC_K       = sapply(results_list, `[[`, "AIC_K"),
    convergence = sapply(results_list, `[[`, "convergence"),
    iter        = sapply(results_list, `[[`, "iter")
  )

  if (nrow(AIC_tab) == 0) stop("No valid models found for npc selection")

  converged_idx <- which(!is.na(AIC_tab$convergence) & AIC_tab$convergence)
  if (length(converged_idx) == 0) {
    warning("No converged models for npc selection; using all finite AIC_K values")
    converged_idx <- seq_len(nrow(AIC_tab))
  }

  valid_idx <- converged_idx[is.finite(AIC_tab$AIC_K[converged_idx])]
  if (length(valid_idx) == 0) stop("All AIC_K values are non-finite among converged models")

  opt_idx <- valid_idx[which.min(AIC_tab$AIC_K[valid_idx])]
  list(npc_opt = AIC_tab$npc[opt_idx], AIC_tab = AIC_tab[order(AIC_tab$npc), ])
}


# Select optimal lambda via LRPsOCV
select_lambda_LRPsOCV <- function(
    data, nbasis_val, npc_val, lambda_seq, tau, grid_size,
    basis_type, solver, solver_control, init_coef, mu_nbasis,
    R = 1, h = 1, m_min = 3,
    seed = NULL, verbose = FALSE,
    parallel = FALSE, n_cores = NULL, cl = NULL,
    scores_weights = NULL
) {
  if (length(lambda_seq) < 2)
    stop("lambda_seq must contain at least 2 values for LRPsOCV")

  N          <- data$n
  basis_full <- get_basis_terms(nbasis_val, basis_type, grid_size, data$data_list)
  Gram_mat   <- basis_full$Gram_mat
  basis_obj  <- basis_full$basis

  LRPsOCV_one_lambda <- function(lam_idx) {
    lam        <- lambda_seq[lam_idx]
    rep_errors <- rep(NA_real_, R)

    for (r in seq_len(R)) {
      if (!is.null(seed)) set.seed(seed + r - 1L)
      data_reduced     <- data
      held_out_indices <- vector("list", N)

      for (i in seq_len(N)) {
        sub_i <- data$data_list[[i]]
        n_i   <- nrow(sub_i)
        if (n_i < m_min) next
        h_i     <- min(h, n_i - 2L)
        I_test  <- sort(sample.int(n_i, h_i))
        held_out_indices[[i]]       <- I_test
        data_reduced$data_list[[i]] <- sub_i[setdiff(seq_len(n_i), I_test), , drop = FALSE]
      }

      fit <- tryCatch(
        fit_model(data_reduced, nbasis_val, npc_val, lam, tau,
                  basis_type, grid_size, init_coef, solver,
                  solver_control, mu_nbasis,
                  scores_weights = scores_weights),
        error = function(e) NULL
      )
      if (is.null(fit)) next

      U_hat <- fit$U; rep_error <- 0; rep_count <- 0L

      for (i in seq_len(N)) {
        I_test <- held_out_indices[[i]]
        if (is.null(I_test) || length(I_test) == 0L) next
        sub_i  <- data$data_list[[i]]
        xi_hat <- fit$Xi[i, ]
        t_test     <- sub_i[I_test, "time", drop = TRUE]
        orthB_test <- fda::eval.basis(t_test, basis_obj) %*% Gram_mat
        y_pred     <- as.vector(orthB_test %*% U_hat %*% xi_hat)
        y_true     <- sub_i[I_test, "centered_value", drop = TRUE]
        rep_error  <- rep_error + sum((y_true - y_pred)^2)
        rep_count  <- rep_count + length(I_test)
      }

      if (rep_count > 0L) rep_errors[r] <- rep_error / rep_count
    }
    rep_errors
  }

  lam_indices <- seq_along(lambda_seq)

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
        "data", "nbasis_val", "npc_val", "lambda_seq", "tau",
        "basis_type", "grid_size", "init_coef", "solver", "solver_control",
        "mu_nbasis", "scores_weights", "R", "h", "m_min", "seed", "N",
        "Gram_mat", "basis_obj"
      ), envir = environment())
      parLapply(cl, lam_indices, LRPsOCV_one_lambda)
    }, finally = {
      if (own_cluster) safe_stop_cluster(cl)
    })
  } else {
    cv_errors_list <- lapply(lam_indices, LRPsOCV_one_lambda)
  }

  cv_errors      <- do.call(cbind, cv_errors_list)
  colnames(cv_errors) <- sprintf("%.2e", lambda_seq)
  rownames(cv_errors) <- paste0("rep_", seq_len(R))

  mean_cv_errors <- colMeans(cv_errors, na.rm = TRUE)
  se_cv_errors   <- apply(cv_errors, 2,
                          function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

  cv_results <- data.frame(
    lambda        = lambda_seq,
    mean_cv_error = mean_cv_errors,
    se_cv_error   = se_cv_errors,
    n_converged   = apply(cv_errors, 2, function(x) sum(!is.na(x)))
  )

  valid_idx  <- which(is.finite(cv_results$mean_cv_error))
  if (length(valid_idx) == 0L) stop("All lambda candidates failed LRPsOCV")

  opt_idx    <- valid_idx[which.min(cv_results$mean_cv_error[valid_idx])]
  lambda_opt <- cv_results$lambda[opt_idx]

  if (verbose)
    cat(sprintf("\n--- LRPsOCV Results ---\nOptimal lambda: %.2e (CV: %.4f +/- %.4f)\n",
                lambda_opt, cv_results$mean_cv_error[opt_idx],
                cv_results$se_cv_error[opt_idx]))

  list(lambda_opt = lambda_opt, cv_results = cv_results, cv_errors = cv_errors)
}


# Select optimal nbasis via LRPsOCV
select_nbasis_LRPsOCV <- function(
    data, nbasis_seq, npc_val, lambda_val, tau, grid_size,
    basis_type, solver, solver_control, init_coef, mu_nbasis,
    R = 1, h = 1, m_min = 3,
    seed = NULL, verbose = FALSE,
    parallel = FALSE, n_cores = NULL, cl = NULL,
    scores_weights = NULL
) {
  N <- data$n

  LRPsOCV_one_nbasis <- function(nb_idx) {
    nb         <- nbasis_seq[nb_idx]
    basis_full <- tryCatch(
      get_basis_terms(nb, basis_type, grid_size, data$data_list),
      error = function(e) NULL
    )
    if (is.null(basis_full)) return(rep(NA_real_, R))

    Gram_mat  <- basis_full$Gram_mat
    basis_obj <- basis_full$basis
    rep_errors <- rep(NA_real_, R)

    for (r in seq_len(R)) {
      if (!is.null(seed)) set.seed(seed + r - 1L)
      data_reduced     <- data
      held_out_indices <- vector("list", N)

      for (i in seq_len(N)) {
        sub_i <- data$data_list[[i]]
        n_i   <- nrow(sub_i)
        if (n_i < m_min) next
        h_i     <- min(h, n_i - 2L)
        I_test  <- sort(sample.int(n_i, h_i))
        held_out_indices[[i]]       <- I_test
        data_reduced$data_list[[i]] <- sub_i[setdiff(seq_len(n_i), I_test), , drop = FALSE]
      }

      fit <- tryCatch(
        fit_model(data_reduced, nb, npc_val, lambda_val, tau,
                  basis_type, grid_size, init_coef, solver,
                  solver_control, mu_nbasis,
                  scores_weights = scores_weights),
        error = function(e) NULL
      )
      if (is.null(fit)) next

      U_hat <- fit$U; rep_error <- 0; rep_count <- 0L

      for (i in seq_len(N)) {
        I_test <- held_out_indices[[i]]
        if (is.null(I_test) || length(I_test) == 0L) next
        sub_i  <- data$data_list[[i]]
        xi_hat <- fit$Xi[i, ]
        t_test     <- sub_i[I_test, "time", drop = TRUE]
        orthB_test <- fda::eval.basis(t_test, basis_obj) %*% Gram_mat
        y_pred     <- as.vector(orthB_test %*% U_hat %*% xi_hat)
        y_true     <- sub_i[I_test, "centered_value", drop = TRUE]
        rep_error  <- rep_error + sum((y_true - y_pred)^2)
        rep_count  <- rep_count + length(I_test)
      }

      if (rep_count > 0L) rep_errors[r] <- rep_error / rep_count
    }
    rep_errors
  }

  nb_indices <- seq_along(nbasis_seq)

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
        "data", "nbasis_seq", "npc_val", "lambda_val", "tau",
        "basis_type", "grid_size", "init_coef", "solver", "solver_control",
        "mu_nbasis", "scores_weights", "R", "h", "m_min", "seed", "N"
      ), envir = environment())
      parLapply(cl, nb_indices, LRPsOCV_one_nbasis)
    }, finally = {
      if (own_cluster) safe_stop_cluster(cl)
    })
  } else {
    cv_errors_list <- lapply(nb_indices, LRPsOCV_one_nbasis)
  }

  cv_errors <- do.call(cbind, cv_errors_list)
  colnames(cv_errors) <- as.character(nbasis_seq)
  rownames(cv_errors) <- paste0("rep_", seq_len(R))

  mean_cv_errors <- colMeans(cv_errors, na.rm = TRUE)
  se_cv_errors   <- apply(cv_errors, 2,
                          function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

  cv_results <- data.frame(
    nbasis        = nbasis_seq,
    mean_cv_error = mean_cv_errors,
    se_cv_error   = se_cv_errors,
    n_converged   = apply(cv_errors, 2, function(x) sum(!is.na(x)))
  )

  valid_idx <- which(!is.na(cv_results$mean_cv_error))
  if (length(valid_idx) == 0) stop("All nbasis values failed CV")

  opt_idx    <- valid_idx[which.min(cv_results$mean_cv_error[valid_idx])]
  nbasis_opt <- cv_results$nbasis[opt_idx]

  if (verbose)
    cat(sprintf("\n--- LRPsOCV Results ---\nOptimal nbasis: %d (CV: %.4f +/- %.4f)\n",
                nbasis_opt, cv_results$mean_cv_error[opt_idx],
                cv_results$se_cv_error[opt_idx]))

  list(nbasis_opt = nbasis_opt, cv_results = cv_results, cv_errors = cv_errors)
}


# Joint 2-D LRPsOCV grid search
#
# mode = "nbasis_npc"  : grid over nbasis_seq x npc_seq, lambda fixed at fixed_lambda (1e-6)
# mode = "npc_lambda"  : grid over npc_seq x lambda_seq, nbasis fixed at fixed_nbasis
select_grid_LRPsOCV <- function(
    data, npc_seq, vary_seq, tau, grid_size,
    basis_type, solver, solver_control, init_coef, mu_nbasis,
    mode = "npc_lambda",
    fixed_nbasis = NULL,
    fixed_lambda = 1e-6,
    R = 1, h = 1, m_min = 3,
    seed = NULL, verbose = FALSE,
    parallel = FALSE, n_cores = NULL, cl = NULL,
    scores_weights = NULL
) {
  stopifnot(mode %in% c("nbasis_npc", "npc_lambda"))
  if (mode == "npc_lambda" && is.null(fixed_nbasis))
    stop("select_grid_LRPsOCV: fixed_nbasis must be supplied for mode = 'npc_lambda'")

  N <- data$n

  # Build feasible grid
  if (mode == "nbasis_npc") {
    grid_pairs <- expand.grid(nbasis = vary_seq, npc = npc_seq,
                              KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    grid_pairs <- grid_pairs[grid_pairs$npc <= grid_pairs$nbasis, , drop = FALSE]
  } else {
    grid_pairs <- expand.grid(npc = npc_seq, lambda = vary_seq,
                              KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    grid_pairs <- grid_pairs[grid_pairs$npc <= fixed_nbasis, , drop = FALSE]
  }

  if (nrow(grid_pairs) == 0) stop("No feasible grid pairs")
  n_pairs <- nrow(grid_pairs)
  row.names(grid_pairs) <- NULL

  # Precompute fixed basis for npc_lambda (same nbasis every pair)
  if (mode == "npc_lambda") {
    basis_fixed <- tryCatch(
      get_basis_terms(fixed_nbasis, basis_type, grid_size, data$data_list),
      error = function(e) NULL
    )
    if (is.null(basis_fixed))
      stop("select_grid_LRPsOCV: failed to compute basis for fixed_nbasis")
  } else {
    basis_fixed <- NULL
  }

  LRPsOCV_one_pair <- function(pair_idx) {
    row <- grid_pairs[pair_idx, ]

    if (mode == "nbasis_npc") {
      nb  <- row$nbasis; np <- row$npc; lam <- fixed_lambda
      basis_eval <- tryCatch(
        get_basis_terms(nb, basis_type, grid_size, data$data_list),
        error = function(e) NULL
      )
    } else {
      nb  <- fixed_nbasis; np <- row$npc; lam <- row$lambda
      basis_eval <- basis_fixed
    }

    if (is.null(basis_eval) || np > nb) return(rep(NA_real_, R))
    Gram_mat  <- basis_eval$Gram_mat
    basis_obj <- basis_eval$basis

    rep_errors <- rep(NA_real_, R)

    for (r in seq_len(R)) {
      if (!is.null(seed)) set.seed(seed + r - 1L)
      data_reduced     <- data
      held_out_indices <- vector("list", N)

      for (i in seq_len(N)) {
        sub_i <- data$data_list[[i]]
        n_i   <- nrow(sub_i)
        if (n_i < m_min) next
        h_i     <- min(h, n_i - 2L)
        I_test  <- sort(sample.int(n_i, h_i))
        held_out_indices[[i]]       <- I_test
        data_reduced$data_list[[i]] <- sub_i[setdiff(seq_len(n_i), I_test), , drop = FALSE]
      }

      fit <- tryCatch(
        fit_model(data_reduced, nb, np, lam, tau, basis_type, grid_size,
                  init_coef, solver, solver_control, mu_nbasis,
                  scores_weights = scores_weights),
        error = function(e) NULL
      )
      if (is.null(fit)) next

      U_hat <- fit$U; rep_error <- 0; rep_count <- 0L

      for (i in seq_len(N)) {
        I_test <- held_out_indices[[i]]
        if (is.null(I_test) || length(I_test) == 0L) next
        sub_i  <- data$data_list[[i]]
        xi_hat <- fit$Xi[i, ]
        t_test     <- sub_i[I_test, "time", drop = TRUE]
        orthB_test <- fda::eval.basis(t_test, basis_obj) %*% Gram_mat
        y_pred     <- as.vector(orthB_test %*% U_hat %*% xi_hat)
        y_true     <- sub_i[I_test, "centered_value", drop = TRUE]
        rep_error  <- rep_error + sum((y_true - y_pred)^2)
        rep_count  <- rep_count + length(I_test)
      }

      if (rep_count > 0L) rep_errors[r] <- rep_error / rep_count
    }
    rep_errors
  }

  pair_indices <- seq_len(n_pairs)

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
        "data", "tau", "basis_type", "grid_size", "init_coef", "solver",
        "solver_control", "mu_nbasis", "scores_weights",
        "R", "h", "m_min", "seed", "N", "mode",
        "fixed_nbasis", "fixed_lambda", "grid_pairs", "basis_fixed"
      ), envir = environment())
      parLapply(cl, pair_indices, LRPsOCV_one_pair)
    }, finally = {
      if (own_cluster) safe_stop_cluster(cl)
    })
  } else {
    cv_errors_list <- lapply(pair_indices, LRPsOCV_one_pair)
  }

  mean_cv <- sapply(cv_errors_list, function(e) mean(e, na.rm = TRUE))
  se_cv   <- sapply(cv_errors_list, function(e) {
    v <- e[!is.na(e)]
    if (length(v) < 2) NA_real_ else sd(v) / sqrt(length(v))
  })
  n_conv <- sapply(cv_errors_list, function(e) sum(!is.na(e)))

  cv_results <- cbind(
    grid_pairs,
    data.frame(mean_cv_error = mean_cv, se_cv_error = se_cv,
               n_converged = n_conv),
    stringsAsFactors = FALSE
  )

  valid_idx <- which(is.finite(cv_results$mean_cv_error))
  if (length(valid_idx) == 0L) stop("All grid pairs failed LRPsOCV")

  opt_idx <- valid_idx[which.min(cv_results$mean_cv_error[valid_idx])]

  if (mode == "nbasis_npc") {
    opt         <- list(nbasis = cv_results$nbasis[opt_idx],
                        npc    = cv_results$npc[opt_idx])
    param_desc  <- sprintf("nbasis=%d, npc=%d", opt$nbasis, opt$npc)
  } else {
    opt         <- list(npc    = cv_results$npc[opt_idx],
                        lambda = cv_results$lambda[opt_idx])
    param_desc  <- sprintf("npc=%d, lambda=%.2e", opt$npc, opt$lambda)
  }

  if (verbose)
    cat(sprintf("\n--- Grid LRPsOCV (%s) ---\nOptimal: %s (CV: %.4f)\n",
                mode, param_desc, cv_results$mean_cv_error[opt_idx]))

  list(opt = opt, cv_results = cv_results)
}
