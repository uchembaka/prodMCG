# Validate and subset initial coefficients
init_coeff <- function(init_coef, n, npc, nbasis) {
  if (is.null(init_coef)) return(NULL)
  if (!is.list(init_coef) || is.null(init_coef$U) || is.null(init_coef$Xi)) {
    stop("init_coef must be NULL or a list with components $U and $Xi")
  }
  U <- init_coef$U
  Xi <- init_coef$Xi
  if (nrow(U) != nbasis) stop(sprintf("init_coef$U has %d rows but nbasis = %d", nrow(U), nbasis))
  if (ncol(U) < npc) stop(sprintf("init_coef$U has %d columns but npc = %d", ncol(U), npc))
  if (nrow(Xi) != n) stop(sprintf("init_coef$Xi has %d rows but n = %d", nrow(Xi), n))
  if (ncol(Xi) < npc) stop(sprintf("init_coef$Xi has %d columns but npc = %d", ncol(Xi), npc))
  list(U = U[, 1:npc, drop = FALSE], Xi = Xi[, 1:npc, drop = FALSE])
}

# Fit a single model via Riemannian optimization
fit_model <- function(
    prep_data, nbasis_val, npc_val, lambda_val, tau, basis_type, grid_size,
    init_coef = NULL, solver = "CG",
    solver_control = list(maxit = 50000, trace = 0, Tolerance = 1e-04),
    mu_nbasis = 15
) {
  if (npc_val > nbasis_val) {
    if (getOption("AIC_verbose", FALSE)) {
      cat(sprintf("  Skipping: npc=%d > nbasis=%d\n", npc_val, nbasis_val))
    }
    return(NULL)
  }
  tryCatch({
    basis_data <- get_basis_terms(nbasis_val, basis_type, grid_size, prep_data$data_list)
    if (!is.null(init_coef)) {
      validated_init <- init_coeff(init_coef, prep_data$n, npc_val, nbasis_val)
    } else if (!is.null(prep_data$tau_init)) {
      tau_U <- prep_data$tau_init$U
      tau_Xi <- prep_data$tau_init$Xi
      validated_init <- list(
        U = tau_U[, 1:min(npc_val, ncol(tau_U)), drop = FALSE],
        Xi = tau_Xi[, 1:min(npc_val, ncol(tau_Xi)), drop = FALSE]
      )
      if (nrow(validated_init$U) != nbasis_val || nrow(validated_init$Xi) != prep_data$n) {
        validated_init <- NULL
      }
    } else {
      validated_init <- NULL
    }
    comp_data <- list(
      prep_data = prep_data, basis_data = basis_data,
      K = npc_val, D = nbasis_val, lambda = lambda_val,
      tau = tau, init_coef = validated_init
    )
    solver_mapped <- solver %||% "RCG"
    solver_mapped <- switch(solver_mapped,
      "BFGS" = "RBFGS", "L-BFGS-B" = "LRBFGS", "CG" = "RCG", solver_mapped)
    fit <- Riemann(
      comp_data, method = solver_mapped,
      mani.params = get.manifold.params(),
      solver.params = get.solver.params(
        Max_Iteration = solver_control$maxit %||% 150000,
        Tolerance = solver_control$Tolerance %||% 1e-04,
        DEBUG = solver_control$trace %||% 0
      )
    )
    return(fit)
  }, error = function(e) {
    if (getOption("AIC_verbose", FALSE)) {
      cat(sprintf("  Error: nbasis=%d, npc=%d, lambda=%g: %s\n",
                  nbasis_val, npc_val, lambda_val, e$message))
    }
    return(NULL)
  })
}

# Route parameter selection to appropriate strategy
select_parameters <- function(
    data, nbasis_seq, npc_seq, lambda_seq, tau, grid_size, basis_type,
    init_coef = NULL, solver = "CG",
    solver_control = list(maxit = 150000, trace = 0, Tolerance = 1e-04),
    verbose = FALSE, mu_nbasis = 15, parallel = FALSE, n_cores = NULL,
    cv_nfolds = 5, tau_selection_method = "eigen_order", selection_method = "cv"
) {

  old_verbose <- getOption("AIC_verbose", FALSE)
  options(AIC_verbose = verbose)
  on.exit(options(AIC_verbose = old_verbose), add = TRUE)

  final_solver_control <- list(maxit = 50000, trace = 0, Tolerance = 1e-04)

  mi_vec <- sapply(data$data_list, nrow)
  max_mi <- max(mi_vec)
  npc_seq <- npc_seq[npc_seq <= max_mi]
  if (length(npc_seq) == 0) {
    npc_seq <- max_mi
    warning(sprintf("All npc candidates exceeded max(m_i)=%d; using npc=%d", max_mi, max_mi))
  }

  cl <- NULL
  if (parallel || tau_selection_method %in% c("aic", "cv")) {
    if (is.null(n_cores)) n_cores <- max(1, detectCores() - 1)
    cl <- tryCatch(init_cluster(n_cores), error = function(e) {
      warning(sprintf("Cluster init failed: %s. Falling back to sequential.", e$message))
      NULL
    })
  }
  on.exit(safe_stop_cluster(cl), add = TRUE)

  tau_results <- NULL

  if (length(tau) > 1) {
    if (max(npc_seq) > 1 &&
        npc_seq[length(npc_seq)] > nbasis_seq[length(nbasis_seq)]) {
      stop("Max npc is greater than max nbasis")
    }

    tau_nbasis <- min(10, nbasis_seq[length(nbasis_seq)])
    tau_npc <- npc_seq[length(npc_seq)]
    tau_lambda <- lambda_seq[1]

    if (tau_selection_method == "aic") {
      tau_fit <- select_tau_aic(
        data = data, nbasis = tau_nbasis, npc = tau_npc, lambda = tau_lambda,
        tau_seq = tau, grid_size = grid_size, basis_type = basis_type,
        solver = solver, solver_control = solver_control, init_coef = NULL,
        mu_nbasis = mu_nbasis, verbose = verbose, n_cores = n_cores, cl = cl
      )
    } else if (tau_selection_method == "cv") {
      tau_fit <- select_tau_cv(
        data = data, nbasis = tau_nbasis, npc = tau_npc, lambda = tau_lambda,
        tau_seq = tau, grid_size = grid_size, basis_type = basis_type,
        solver = solver, solver_control = solver_control, init_coef = NULL,
        mu_nbasis = mu_nbasis, nfolds = cv_nfolds, verbose = verbose,
        n_cores = n_cores, cl = cl
      )
    } else {
      tau_fit <- select_tau_eigen_order(
        data = data, nbasis = tau_nbasis, npc = tau_npc, lambda = tau_lambda,
        tau_seq = tau, grid_size = grid_size, basis_type = basis_type,
        solver = solver, solver_control = solver_control, init_coef = NULL,
        mu_nbasis = mu_nbasis, verbose = verbose, parallel = parallel,
        n_cores = n_cores, cl = cl
      )
    }

    tau <- tau_fit$tau
    tau_results <- tau_fit$tau_results
    if (!is.null(tau_fit$fit)) {
      data$tau_init <- list(U = tau_fit$fit$U, Xi = tau_fit$fit$Xi)
    }
  } else {
    tau <- tau[1]
  }

  nbasis_opt <- nbasis_seq[1]
  npc_opt <- npc_seq[1]
  lambda_opt <- lambda_seq[1]
  cv_results <- NULL
  AIC_tab <- data.frame()
  AIC_tab_param <- NULL

  multi_nbasis <- length(nbasis_seq) > 1
  multi_lambda <- length(lambda_seq) > 1
  multi_npc <- length(npc_seq) > 1

  if (multi_nbasis) {
    if (multi_lambda || (length(lambda_seq) == 1 && lambda_seq[1] > 1e-6)) {
      warning("Setting lambda = 1e-6 for nbasis selection")
    }
    lambda_fixed <- 1e-6
    npc_for_sel <- if (multi_npc) max(npc_seq) else npc_seq[1]

    if (selection_method == "aic") {
      # if (verbose) cat("=== Selecting nbasis via AIC ===\n")
      aic_nbasis_res <- select_nbasis_aic(
        data = data, nbasis_seq = nbasis_seq, npc_val = npc_for_sel,
        lambda_val = lambda_fixed, tau = tau, grid_size = grid_size,
        basis_type = basis_type, solver = solver, solver_control = solver_control,
        init_coef = init_coef, mu_nbasis = mu_nbasis, verbose = verbose,
        parallel = parallel, n_cores = n_cores, cl = cl
      )
      nbasis_opt <- aic_nbasis_res$nbasis_opt
      lambda_opt <- lambda_fixed
      AIC_tab_param <- aic_nbasis_res$AIC_tab
    } else {
      # if (verbose) cat("=== Selecting nbasis via CV ===\n")
      cv_res <- select_nbasis_cv(
        data = data, nbasis_seq = nbasis_seq, npc_val = npc_for_sel,
        lambda_val = lambda_fixed, tau = tau, grid_size = grid_size,
        basis_type = basis_type, solver = solver, solver_control = solver_control,
        init_coef = init_coef, mu_nbasis = mu_nbasis, nfolds = cv_nfolds,
        parallel = parallel, n_cores = n_cores, verbose = verbose, cl = cl
      )
      nbasis_opt <- cv_res$nbasis_opt
      lambda_opt <- lambda_fixed
      cv_results <- cv_res$cv_results
    }

    if (multi_npc) {
      # if (verbose) cat(sprintf("=== Selecting npc via AIC_K (nbasis=%d) ===\n", nbasis_opt))
      aic_res <- select_npc_aic(
        data = data, nbasis_val = nbasis_opt, npc_seq = npc_seq,
        lambda_val = lambda_opt, tau = tau, grid_size = grid_size,
        basis_type = basis_type, solver = solver, solver_control = solver_control,
        init_coef = init_coef, mu_nbasis = mu_nbasis, verbose = verbose,
        parallel = parallel, n_cores = n_cores, cl = cl
      )
      npc_opt <- aic_res$npc_opt
      AIC_tab <- aic_res$AIC_tab
    }

  } else if (multi_lambda) {
    if (nbasis_seq[1] < 5) warning(sprintf("Low nbasis (%d). Consider increasing", nbasis_seq[1]))
    npc_for_sel <- if (multi_npc) max(npc_seq) else npc_seq[1]

    if (selection_method == "aic") {
      # if (verbose) cat("=== Selecting lambda via AIC ===\n")
      aic_lambda_res <- select_lambda_aic(
        data = data, nbasis_val = nbasis_seq[1], npc_val = npc_for_sel,
        lambda_seq = lambda_seq, tau = tau, grid_size = grid_size,
        basis_type = basis_type, solver = solver, solver_control = solver_control,
        init_coef = init_coef, mu_nbasis = mu_nbasis, verbose = verbose,
        parallel = parallel, n_cores = n_cores, cl = cl
      )
      lambda_opt <- aic_lambda_res$lambda_opt
      AIC_tab_param <- aic_lambda_res$AIC_tab
    } else {
      # if (verbose) cat("=== Selecting lambda via CV ===\n")
      cv_res <- select_lambda_cv(
        data = data, nbasis_val = nbasis_seq[1], npc_val = npc_for_sel,
        lambda_seq = lambda_seq, tau = tau, grid_size = grid_size,
        basis_type = basis_type, solver = solver, solver_control = solver_control,
        init_coef = init_coef, mu_nbasis = mu_nbasis, nfolds = cv_nfolds,
        parallel = parallel, n_cores = n_cores, verbose = verbose, cl = cl
      )
      lambda_opt <- cv_res$lambda_opt
      cv_results <- cv_res$cv_results
    }

    if (multi_npc) {
      # if (verbose) cat(sprintf("=== Selecting npc via AIC_K (lambda=%.2e) ===\n", lambda_opt))
      aic_res <- select_npc_aic(
        data = data, nbasis_val = nbasis_seq[1], npc_seq = npc_seq,
        lambda_val = lambda_opt, tau = tau, grid_size = grid_size,
        basis_type = basis_type, solver = solver, solver_control = solver_control,
        init_coef = init_coef, mu_nbasis = mu_nbasis, verbose = verbose,
        parallel = parallel, n_cores = n_cores, cl = cl
      )
      npc_opt <- aic_res$npc_opt
      AIC_tab <- aic_res$AIC_tab
    }

  } else if (multi_npc) {
    # if (verbose) cat("=== Selecting npc via AIC_K ===\n")
    aic_res <- select_npc_aic(
      data = data, nbasis_val = nbasis_seq[1], npc_seq = npc_seq,
      lambda_val = lambda_seq[1], tau = tau, grid_size = grid_size,
      basis_type = basis_type, solver = solver, solver_control = solver_control,
      init_coef = init_coef, mu_nbasis = mu_nbasis, verbose = verbose,
      parallel = parallel, n_cores = n_cores, cl = cl
    )
    npc_opt <- aic_res$npc_opt
    AIC_tab <- aic_res$AIC_tab
  }

  # if (verbose) {
  #   cat(sprintf("\n=== Final Fit: nbasis=%d, npc=%d, lambda=%.2e ===\n",
  #               nbasis_opt, npc_opt, lambda_opt))
  # }

  final_fit <- fit_model(
    prep_data = data, nbasis_val = nbasis_opt, npc_val = npc_opt,
    lambda_val = lambda_opt, tau = tau, basis_type = basis_type,
    grid_size = grid_size, init_coef = init_coef, solver = solver,
    solver_control = final_solver_control, mu_nbasis = mu_nbasis
  )

  if (is.null(final_fit)) stop("Final model fitting failed with optimal parameters")

  if (nrow(AIC_tab) == 0) {
    AIC_tab <- data.frame(
      nbasis = nbasis_opt, npc = npc_opt, lambda = lambda_opt,
      fn = final_fit$fval, AIC = final_fit$AIC,
      AIC_K = if (!is.null(final_fit$AIC_K)) final_fit$AIC_K else NA,
      convergence = final_fit$convergence,
      aDF = final_fit$aDF, iter = final_fit$iter
    )
  }

  list(
    fit = final_fit, nbasis_opt = nbasis_opt, npc_opt = npc_opt,
    lambda_opt = lambda_opt, tau = tau, AIC_tab = AIC_tab,
    AIC_tab_param = AIC_tab_param, cv_results = cv_results,
    tau_results = tau_results
  )
}
