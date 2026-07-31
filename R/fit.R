# Validate and subset initial coefficients
init_coeff <- function(init_coef, n, npc, nbasis) {
  if (is.null(init_coef)) return(NULL)
  if (!is.list(init_coef) || is.null(init_coef$U) || is.null(init_coef$Xi))
    stop("init_coef must be NULL or a list with components $U and $Xi")
  U <- init_coef$U; Xi <- init_coef$Xi
  if (nrow(U)  != nbasis) stop(sprintf("init_coef$U has %d rows but nbasis = %d",  nrow(U),  nbasis))
  if (ncol(U)  <  npc)   stop(sprintf("init_coef$U has %d columns but npc = %d",   ncol(U),  npc))
  if (nrow(Xi) != n)     stop(sprintf("init_coef$Xi has %d rows but n = %d",        nrow(Xi), n))
  if (ncol(Xi) <  npc)   stop(sprintf("init_coef$Xi has %d columns but npc = %d",  ncol(Xi), npc))
  list(U = U[, 1:npc, drop = FALSE], Xi = Xi[, 1:npc, drop = FALSE])
}


# Fit a single model via Riemannian optimization
fit_model <- function(
    prep_data, nbasis_val, npc_val, lambda_val, tau, basis_type, grid_size,
    init_coef = NULL, solver = "CG",
    solver_control = list(maxit = 50000, trace = 0, Tolerance = 1e-04),
    mu_nbasis = 15, scores_weights = NULL
) {
  if (npc_val > nbasis_val) {
    if (getOption("AIC_verbose", FALSE))
      warning(sprintf("  Skipping: npc=%d > nbasis=%d\n", npc_val, nbasis_val))
    return(NULL)
  }

  tryCatch({
    basis_data <- get_basis_terms(nbasis_val, basis_type, grid_size, prep_data$data_list)

    if (!is.null(init_coef)) {
      validated_init <- init_coeff(init_coef, prep_data$n, npc_val, nbasis_val)
    } else if (!is.null(prep_data$tau_init)) {
      tau_U  <- prep_data$tau_init$U
      tau_Xi <- prep_data$tau_init$Xi
      validated_init <- list(
        U  = tau_U[,  1:min(npc_val, ncol(tau_U)),  drop = FALSE],
        Xi = tau_Xi[, 1:min(npc_val, ncol(tau_Xi)), drop = FALSE]
      )
      if (nrow(validated_init$U) != nbasis_val || nrow(validated_init$Xi) != prep_data$n)
        validated_init <- NULL
    } else {
      validated_init <- NULL
    }

    comp_data <- list(
      prep_data = prep_data, basis_data = basis_data,
      K = npc_val, D = nbasis_val, lambda = lambda_val,
      tau = tau, scores_weights = scores_weights, init_coef = validated_init
    )

    solver_mapped <- solver %||% "RCG"
    solver_mapped <- switch(solver_mapped,
                            "BFGS"     = "RBFGS",
                            "L-BFGS-B" = "LRBFGS",
                            "CG"       = "RCG",
                            solver_mapped)

    Riemann(
      comp_data, method = solver_mapped,
      mani.params   = get.manifold.params(),
      solver.params = get.solver.params(
        Max_Iteration = solver_control$maxit     %||% 150000,
        Tolerance     = solver_control$Tolerance  %||% 1e-04,
        DEBUG         = solver_control$trace      %||% 0
      )
    )
  }, error = function(e) {
    if (getOption("AIC_verbose", FALSE))
      cat(sprintf("  Error: nbasis=%d, npc=%d, lambda=%g: %s\n",
                  nbasis_val, npc_val, lambda_val, e$message))
    NULL
  })
}


# Route parameter selection to the appropriate strategy
#
# Solver-control policy:
#   tau_selection   -> user's solver_control  (cold start)
#   selection loops -> INTERNAL when tau_init warm start is available;
#                      user's solver_control when tau is fixed (cold)
#   final fit       -> INTERNAL when any search was done;
#                      user's solver_control when all params are scalars
#
# grid_search (logical in LRPsOCV_parameters):
#   TRUE  + vector nbasis -> joint nbasis × npc LRPsOCV ("nbasis_npc" mode)
#   TRUE  + vector lambda -> joint npc × lambda LRPsOCV ("npc_lambda" mode)
#   FALSE -> sequential selection (default)
select_parameters <- function(
    data, nbasis_seq, npc_seq, lambda_seq, tau, grid_size, basis_type,
    init_coef = NULL, solver = "CG",
    solver_control = list(maxit = 150000, trace = 0, Tolerance = 1e-04),
    verbose = FALSE, mu_nbasis = 15, parallel = TRUE, n_cores = NULL,
    LRPsOCV_parameters = list(R = 1, h = 1, m_min = 3),
    seed = NULL, scores_weights = NULL
) {
  INTERNAL_CTRL <- list(maxit = 50000, trace = 0, Tolerance = 1e-04)

  old_verbose <- getOption("AIC_verbose", FALSE)
  options(AIC_verbose = verbose)
  on.exit(options(AIC_verbose = old_verbose), add = TRUE)

  # Clip npc to maximum feasible per curve
  max_mi  <- max(sapply(data$data_list, nrow))
  npc_seq <- npc_seq[npc_seq <= max_mi]
  if (length(npc_seq) == 0) {
    npc_seq <- max_mi
    warning(sprintf("All npc candidates exceeded max(m_i)=%d; using npc=%d", max_mi, max_mi))
  }

  cl <- NULL
  if (parallel) {
    if (is.null(n_cores)) n_cores <- max(1, detectCores() - 1)
    cl <- tryCatch(init_cluster(n_cores), error = function(e) {
      warning(sprintf("Cluster init failed: %s. Falling back to sequential.", e$message))
      NULL
    })
  }
  on.exit(safe_stop_cluster(cl), add = TRUE)

  # --- Tau selection (always uses user's solver_control — cold start) ---
  tau_results      <- NULL
  tau_was_searched <- length(tau) > 1

  if (tau_was_searched) {
    if (max(npc_seq) > nbasis_seq[length(nbasis_seq)])
      stop("Max npc is greater than max nbasis")

    tau_fit <- select_tau_eigen_order(
      data = data,
      nbasis = min(10, nbasis_seq[length(nbasis_seq)]),
      npc    = npc_seq[length(npc_seq)],
      lambda = 1e-6, tau_seq = tau,
      grid_size = grid_size, basis_type = basis_type,
      solver = solver, solver_control = solver_control,
      init_coef = NULL, mu_nbasis = mu_nbasis,
      verbose = verbose, parallel = parallel, n_cores = n_cores,
      cl = cl, scores_weights = scores_weights
    )
    tau         <- tau_fit$tau
    tau_results <- tau_fit$tau_results
    if (!is.null(tau_fit$fit))
      data$tau_init <- list(U = tau_fit$fit$U, Xi = tau_fit$fit$Xi)
  } else {
    tau <- tau[1]
  }

  # --- Solver controls ---
  has_tau_warm <- tau_was_searched && !is.null(data$tau_init)
  loop_ctrl    <- if (has_tau_warm) INTERNAL_CTRL else solver_control

  multi_nbasis <- length(nbasis_seq) > 1
  multi_lambda <- length(lambda_seq)  > 1
  multi_npc    <- length(npc_seq)     > 1
  any_search   <- multi_nbasis || multi_lambda || multi_npc

  final_ctrl   <- if (any_search || tau_was_searched) INTERNAL_CTRL else solver_control

  # LRPsOCV options
  lrps_R    <- LRPsOCV_parameters$R      %||% 1
  lrps_h    <- LRPsOCV_parameters$h      %||% 1
  lrps_m    <- LRPsOCV_parameters$m_min  %||% 3
  use_grid  <- isTRUE(LRPsOCV_parameters$grid_search)

  nbasis_opt       <- nbasis_seq[1]
  npc_opt          <- npc_seq[1]
  lambda_opt       <- lambda_seq[1]
  AIC_K_tab        <- data.frame()
  lambda_selection <- NULL
  nbasis_selection <- NULL
  grid_cv_results  <- NULL

  # -----------------------------------------------------------------------
  # Both nbasis and lambda are vectors:
  #   always warn + fix lambda = 1e-6 + select nbasis only.
  #   grid_search is not supported in this case.
  # -----------------------------------------------------------------------
  if (multi_nbasis && multi_lambda) {
    warning(paste(
      "Both nbasis and lambda are vectors.",
      "Fixing lambda = 1e-6 and selecting nbasis only via LRPsOCV.",
      "Supply only one as a vector to enable lambda selection or grid search."
    ))
    lambda_opt  <- 1e-6
    npc_for_sel <- if (multi_npc) max(npc_seq) else npc_seq[1]

    cv_res <- select_nbasis_LRPsOCV(
      data = data, nbasis_seq = nbasis_seq, npc_val = npc_for_sel,
      lambda_val = lambda_opt, tau = tau, grid_size = grid_size,
      basis_type = basis_type, solver = solver, solver_control = loop_ctrl,
      init_coef = init_coef, mu_nbasis = mu_nbasis,
      R = lrps_R, h = lrps_h, m_min = lrps_m, seed = seed,
      parallel = parallel, n_cores = n_cores, verbose = verbose,
      cl = cl, scores_weights = scores_weights
    )
    nbasis_opt       <- cv_res$nbasis_opt
    nbasis_selection <- cv_res$cv_results

    if (multi_npc) {
      aic_k_res <- select_npc_aic(
        data = data, nbasis_val = nbasis_opt, npc_seq = npc_seq,
        lambda_val = lambda_opt, tau = tau, grid_size = grid_size,
        basis_type = basis_type, solver = solver, solver_control = loop_ctrl,
        init_coef = init_coef, mu_nbasis = mu_nbasis, verbose = verbose,
        parallel = parallel, n_cores = n_cores, cl = cl,
        scores_weights = scores_weights, seed = seed
      )
      npc_opt   <- aic_k_res$npc_opt
      AIC_K_tab <- aic_k_res$AIC_tab
    }

  # -----------------------------------------------------------------------
  # Only nbasis is a vector
  # -----------------------------------------------------------------------
  } else if (multi_nbasis) {
    lambda_opt <- 1e-6   # always use 1e-6 for nbasis selection

    if (use_grid && multi_npc) {
      # Joint nbasis × npc grid search via LRPsOCV
      if (verbose) cat("=== Grid LRPsOCV: nbasis x npc (lambda = 1e-6) ===\n")
      grid_res <- select_grid_LRPsOCV(
        data = data, npc_seq = npc_seq, vary_seq = nbasis_seq,
        tau = tau, grid_size = grid_size, basis_type = basis_type,
        solver = solver, solver_control = loop_ctrl, init_coef = init_coef,
        mu_nbasis = mu_nbasis, mode = "nbasis_npc",
        fixed_lambda = lambda_opt,
        R = lrps_R, h = lrps_h, m_min = lrps_m, seed = seed,
        verbose = verbose, parallel = parallel, n_cores = n_cores,
        cl = cl, scores_weights = scores_weights
      )
      nbasis_opt      <- grid_res$opt$nbasis
      npc_opt         <- grid_res$opt$npc
      grid_cv_results <- grid_res$cv_results

    } else {
      if (use_grid && !multi_npc)
        warning("grid_search = TRUE but npc is scalar; performing sequential nbasis selection.")

      npc_for_sel <- if (multi_npc) max(npc_seq) else npc_seq[1]
      cv_res <- select_nbasis_LRPsOCV(
        data = data, nbasis_seq = nbasis_seq, npc_val = npc_for_sel,
        lambda_val = lambda_opt, tau = tau, grid_size = grid_size,
        basis_type = basis_type, solver = solver, solver_control = loop_ctrl,
        init_coef = init_coef, mu_nbasis = mu_nbasis,
        R = lrps_R, h = lrps_h, m_min = lrps_m, seed = seed,
        parallel = parallel, n_cores = n_cores, verbose = verbose,
        cl = cl, scores_weights = scores_weights
      )
      nbasis_opt       <- cv_res$nbasis_opt
      nbasis_selection <- cv_res$cv_results

      if (multi_npc) {
        aic_k_res <- select_npc_aic(
          data = data, nbasis_val = nbasis_opt, npc_seq = npc_seq,
          lambda_val = lambda_opt, tau = tau, grid_size = grid_size,
          basis_type = basis_type, solver = solver, solver_control = loop_ctrl,
          init_coef = init_coef, mu_nbasis = mu_nbasis, verbose = verbose,
          parallel = parallel, n_cores = n_cores, cl = cl,
          scores_weights = scores_weights, seed = seed
        )
        npc_opt   <- aic_k_res$npc_opt
        AIC_K_tab <- aic_k_res$AIC_tab
      }
    }

  # -----------------------------------------------------------------------
  # Only lambda is a vector
  # -----------------------------------------------------------------------
  } else if (multi_lambda) {
    if (nbasis_seq[1] < 5)
      warning(sprintf("Low nbasis (%d) for lambda selection. Consider increasing.", nbasis_seq[1]))

    if (use_grid && multi_npc) {
      # Joint npc × lambda grid search via LRPsOCV
      if (verbose) cat("=== Grid LRPsOCV: npc x lambda ===\n")
      grid_res <- select_grid_LRPsOCV(
        data = data, npc_seq = npc_seq, vary_seq = lambda_seq,
        tau = tau, grid_size = grid_size, basis_type = basis_type,
        solver = solver, solver_control = loop_ctrl, init_coef = init_coef,
        mu_nbasis = mu_nbasis, mode = "npc_lambda",
        fixed_nbasis = nbasis_seq[1],
        R = lrps_R, h = lrps_h, m_min = lrps_m, seed = seed,
        verbose = verbose, parallel = parallel, n_cores = n_cores,
        cl = cl, scores_weights = scores_weights
      )
      npc_opt         <- grid_res$opt$npc
      lambda_opt      <- grid_res$opt$lambda
      grid_cv_results <- grid_res$cv_results

    } else {
      if (use_grid && !multi_npc)
        warning("grid_search = TRUE but npc is scalar; performing sequential lambda selection.")

      npc_for_sel <- if (multi_npc) max(npc_seq) else npc_seq[1]
      cv_res     <- select_lambda_LRPsOCV(
        data = data, nbasis_val = nbasis_seq[1], npc_val = npc_for_sel,
        lambda_seq = lambda_seq, tau = tau, grid_size = grid_size,
        basis_type = basis_type, solver = solver, solver_control = loop_ctrl,
        init_coef = init_coef, mu_nbasis = mu_nbasis,
        R = lrps_R, h = lrps_h, m_min = lrps_m, seed = seed,
        parallel = parallel, n_cores = n_cores, verbose = verbose,
        cl = cl, scores_weights = scores_weights
      )
      lambda_opt       <- cv_res$lambda_opt
      lambda_selection <- cv_res$cv_results

      if (multi_npc) {
        aic_k_res <- select_npc_aic(
          data = data, nbasis_val = nbasis_seq[1], npc_seq = npc_seq,
          lambda_val = lambda_opt, tau = tau, grid_size = grid_size,
          basis_type = basis_type, solver = solver, solver_control = loop_ctrl,
          init_coef = init_coef, mu_nbasis = mu_nbasis, verbose = verbose,
          parallel = parallel, n_cores = n_cores, cl = cl,
          scores_weights = scores_weights, seed = seed
        )
        npc_opt   <- aic_k_res$npc_opt
        AIC_K_tab <- aic_k_res$AIC_tab
      }
    }

  # -----------------------------------------------------------------------
  # Only npc is a vector
  # -----------------------------------------------------------------------
  } else if (multi_npc) {
    aic_k_res <- select_npc_aic(
      data = data, nbasis_val = nbasis_seq[1], npc_seq = npc_seq,
      lambda_val = lambda_seq[1], tau = tau, grid_size = grid_size,
      basis_type = basis_type, solver = solver, solver_control = loop_ctrl,
      init_coef = init_coef, mu_nbasis = mu_nbasis, verbose = verbose,
      parallel = parallel, n_cores = n_cores, cl = cl,
      scores_weights = scores_weights, seed = seed
    )
    npc_opt   <- aic_k_res$npc_opt
    AIC_K_tab <- aic_k_res$AIC_tab
  }

  if (verbose)
    cat(sprintf("\n=== Final Fit: nbasis=%d, npc=%d, lambda=%.2e ===\n",
                nbasis_opt, npc_opt, lambda_opt))

  final_fit <- fit_model(
    prep_data = data, nbasis_val = nbasis_opt, npc_val = npc_opt,
    lambda_val = lambda_opt, tau = tau, basis_type = basis_type,
    grid_size = grid_size, init_coef = init_coef, solver = solver,
    solver_control = final_ctrl, mu_nbasis = mu_nbasis,
    scores_weights = scores_weights
  )

  if (is.null(final_fit)) stop("Final model fitting failed with optimal parameters")

  list(
    fit              = final_fit,
    nbasis_opt       = nbasis_opt,
    npc_opt          = npc_opt,
    lambda_opt       = lambda_opt,
    tau              = tau,
    AIC_K_tab        = AIC_K_tab,
    lambda_selection = lambda_selection,
    nbasis_selection = nbasis_selection,
    grid_cv_results  = grid_cv_results,
    tau_results      = tau_results
  )
}
