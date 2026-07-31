#' Reconstruct sparse functional data via Riemannian optimization
#'
#' Main user-facing function that performs curve reconstruction and functional
#' principal component analysis of sparse longitudinal data using Riemannian
#' optimization on the product manifold Stiefel(D,K) x Euclidean(N,K).
#'
#' @param data Matrix or data frame with 3 columns: subject ID, observation
#'   time, observed value.
#' @param nbasis Integer or integer vector. Number of basis functions. If a
#'   vector, optimal nbasis is selected via LRPsOCV (lambda forced to 1e-6).
#' @param npc Integer or integer vector. Number of principal components.
#'   If a vector, optimal npc is selected via AIC (or jointly with nbasis/lambda
#'   via LRPsOCV when \code{LRPsOCV_parameters$grid_search = TRUE}).
#' @param grid_size Integer. Number of grid points for evaluation (default 101).
#' @param basis_type Character. Type of basis: \code{"bspline"} (default) or
#'   \code{"fourier"}.
#' @param lambda Numeric or numeric vector. Smoothing penalty. If a vector,
#'   optimal lambda is selected via LRPsOCV. When both \code{nbasis} and
#'   \code{lambda} are vectors, lambda is fixed to 1e-6 and only nbasis is
#'   selected (with a warning).
#' @param tau Numeric or numeric vector. Ridge penalty. If a vector, optimal
#'   tau is selected via eigenvalue ordering and convergence checking.
#' @param solver Character. Optimization method: \code{"CG"} (default),
#'   \code{"BFGS"}, or \code{"L-BFGS-B"}.
#' @param solver_control List with solver control parameters:
#'   \itemize{
#'     \item \code{maxit}: Maximum iterations (default 150000).
#'     \item \code{trace}: Debug level (default 0).
#'     \item \code{Tolerance}: Convergence tolerance (default 1e-4).
#'   }
#'   The user-supplied control is used for tau selection, parameter selection
#'   with fixed tau, and when all parameters are scalars. Internal fits that
#'   have warm starts from tau selection or parameter search use a faster
#'   internal control (maxit = 50000).
#' @param obs_range Numeric vector of length 2. Observation time range
#'   (default c(0,1)).
#' @param init_coef NULL or list with components U (D x K) and Xi (N x K).
#' @param plot_results Logical. Plot AIC and CV results (default FALSE).
#' @param verbose Logical. Print progress messages (default FALSE).
#' @param mu_nbasis Integer. Number of basis functions for mean estimation
#'   (default 15).
#' @param parallel Logical. Use parallel computation (default FALSE).
#' @param n_cores Integer or NULL. Number of cores; NULL uses detectCores()-1.
#'
#' @param seed Integer or NULL. Random seed set once at the start of the
#'   function, making all downstream RNG usage reproducible: Riemannian
#'   initialization of U and Xi, and point sampling in LRPsOCV replications.
#'   Default NULL (no seed set).
#'
#' @param LRPsOCV_parameters A list controlling Leave-Random-Points-Out CV:
#'   \describe{
#'     \item{\code{R}}{Number of random repetitions (default 1).}
#'     \item{\code{h}}{Observations held out per curve per repetition
#'       (actual hold-out is \code{min(h, n_i - 2)}; default 1).}
#'     \item{\code{m_min}}{Minimum observations for a curve to participate
#'       in the hold-out (default 3).}
#'     \item{\code{grid_search}}{Logical. When \code{TRUE}, performs a joint
#'       2-D LRPsOCV grid search instead of sequential selection. The grid
#'       axes are determined by which parameter is a vector: if \code{nbasis}
#'       is a vector the grid is nbasis x npc; if \code{lambda} is a vector
#'       the grid is npc x lambda. Requires \code{npc} to also be a vector.
#'       Default: \code{FALSE}.}
#'   }
#'
#' @param scores_weights Optional numeric vector of length \code{max(npc)}.
#'   Component-specific ridge-penalty weights. When NULL, weights are set
#'   automatically (1:max(npc) or adaptively from Stage 1 eigenvalues when
#'   \code{two_stage = TRUE}).
#'
#' @param two_stage Logical. When TRUE a two-stage procedure sets adaptive
#'   ridge weights from Stage 1 eigenvalues. Ignored when \code{scores_weights}
#'   is supplied. Default: TRUE.
#'
#' @param share_tau Logical. When TRUE (default) and \code{two_stage = TRUE},
#'   reuse the Stage 1 tau in Stage 2. Requires \code{two_stage = TRUE}.
#'
#' @param bootstrap_parameters List controlling bootstrap CIs:
#'   \itemize{
#'     \item \code{compute}: Logical (default FALSE).
#'     \item \code{B}: Number of replicates (default 200).
#'     \item \code{alpha}: Significance level (default 0.05).
#'     \item \code{method}: \code{"parametric"} (default) or
#'       \code{"residual"}.
#'   }
#'
#' @return A list with components:
#'   \describe{
#'     \item{U}{Orthonormal basis coefficient matrix (D x K).}
#'     \item{Xi}{Score matrix (N x K).}
#'     \item{xHat}{Reconstructed functions (N x G matrix).}
#'     \item{Phi_mat}{Eigenfunction evaluations on grid (G x K).}
#'     \item{eigenfunctions}{List of K eigenfunction fd objects.}
#'     \item{eigenValues}{Eigenvalues (variances of scores).}
#'     \item{mu_function}{Mean function (fd object).}
#'     \item{scores_weights}{Scores weights used.}
#'     \item{sig2}{Residual variance.}
#'     \item{npc}{Selected number of principal components.}
#'     \item{nbasis}{Selected number of basis functions.}
#'     \item{nbasis_selection}{LRPsOCV results for nbasis (if applicable).}
#'     \item{grid}{Evaluation grid.}
#'     \item{AIC}{AIC_K value of the final model.}
#'     \item{AIC_tab}{Data frame of AIC_K values across npc candidates.}
#'     \item{tau}{Selected tau value.}
#'     \item{lambda}{Selected lambda value.}
#'     \item{lambda_selection}{LRPsOCV results for lambda (if applicable).}
#'     \item{grid_cv_results}{2-D LRPsOCV results table (if grid_search used).}
#'     \item{tau_results}{Tau selection table (if tau was a vector).}
#'     \item{Fits}{List: fitted_values, residuals, var_explained.}
#'     \item{Data}{Full data object from the fit.}
#'     \item{manOptim_output}{Raw ManifoldOptim output.}
#'     \item{xHat_CI}{List of N matrices (G x 2: lower/upper bootstrap CIs),
#'       or NULL.}
#'     \item{call_args}{All arguments used.}
#'   }
#'
#' @export
reconstructSparseFuncs <- function(
    data,
    nbasis     = 10,
    npc        = 2:5,
    grid_size  = 101,
    basis_type = "bspline",
    lambda     = exp(-12:0),
    tau        = 10^(-6:-1),
    solver     = "CG",
    solver_control = list(maxit    = 50000,
                          trace    = 0,
                          Tolerance = 1e-04),
    obs_range  = c(0, 1),
    init_coef  = NULL,
    plot_results = TRUE,
    verbose    = FALSE,
    mu_nbasis  = 15,
    parallel   = TRUE,
    n_cores    = NULL,
    LRPsOCV_parameters = list(R = 1, h = 1, m_min = 3, grid_search = FALSE),
    seed       = NULL,
    scores_weights = NULL,
    two_stage  = TRUE,
    share_tau  = TRUE,
    bootstrap_parameters = list(compute = FALSE,
                                B       = 200,
                                alpha   = 0.05,
                                method  = "parametric")
) {
  if (!is.null(seed)) set.seed(seed)

  solver     <- match.arg(solver,     c("BFGS", "L-BFGS-B", "CG"))
  basis_type <- match.arg(basis_type, c("bspline", "fourier"))
  n_cores_org <- n_cores
  if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)

  p_max <- max(npc)

  if (!is.null(scores_weights)) {
    if (!is.numeric(scores_weights))
      stop("scores_weights must be a numeric vector")
    if (length(scores_weights) != p_max)
      stop(sprintf("scores_weights must have length %d (max of npc), got %d",
                   p_max, length(scores_weights)))
  }

  if (two_stage && !is.null(scores_weights)) {
    warning("Both two_stage=TRUE and scores_weights provided; using provided scores_weights (skipping Stage 1)")
    two_stage <- FALSE
  }

  if (share_tau && !two_stage) {
    warning("share_tau=TRUE requires two_stage=TRUE; ignoring share_tau")
    share_tau <- FALSE
  }

  if (!is.list(bootstrap_parameters))
    stop("bootstrap_parameters must be a list with elements: compute, B, alpha, method")

  boot_compute <- isTRUE(bootstrap_parameters$compute)
  boot_B       <- bootstrap_parameters$B      %||% 200
  boot_alpha   <- bootstrap_parameters$alpha  %||% 0.05
  boot_method  <- bootstrap_parameters$method %||% "parametric"

  if (!is.numeric(boot_B) || length(boot_B) != 1 || boot_B < 2)
    stop("bootstrap_parameters$B must be a positive integer >= 2")
  boot_B <- as.integer(boot_B)

  if (!is.numeric(boot_alpha) || length(boot_alpha) != 1 ||
      boot_alpha <= 0 || boot_alpha >= 1)
    stop("bootstrap_parameters$alpha must be a numeric value in (0, 1)")

  boot_method <- match.arg(boot_method, c("parametric", "residual"))

  if (basis_type == "fourier") {
    if (length(nbasis) > 1) {
      nbasis <- nbasis[nbasis %% 2 == 1]
    } else {
      nbasis <- ifelse(nbasis %% 2 == 1, nbasis, nbasis + 1)
    }
    if (verbose) message("nbasis: ", paste(nbasis, collapse = ", "))
  }

  if (!is.null(init_coef)) {
    if (!is.list(init_coef) || is.null(init_coef$U) || is.null(init_coef$Xi))
      stop("init_coef must be NULL or a list with components $U and $Xi")
  }

  prep_data       <- preprocess_data(data, basis_type, grid_size, obs_range, mu_nbasis)
  lambda_original <- lambda

  # --- Two-stage adaptive scores weights ---
  stage1_tau         <- NULL
  stage1_tau_results <- NULL

  if (two_stage) {
    stage1_nbasis <- min(nbasis)
    stage1_npc    <- p_max
    if (stage1_npc > stage1_nbasis) stage1_nbasis <- stage1_npc

    if (length(tau) > 1) {
      if (verbose)
        cat(sprintf("=== Stage 1: Tau selection (nbasis=%d, npc=%d, lambda=1e-6) ===\n",
                    stage1_nbasis, stage1_npc))
      tau_sel <- select_tau_eigen_order(
        data = prep_data, nbasis = stage1_nbasis, npc = stage1_npc,
        lambda = 1e-6, tau_seq = tau, grid_size = grid_size,
        basis_type = basis_type, solver = solver,
        solver_control = solver_control, init_coef = NULL,
        mu_nbasis = mu_nbasis, verbose = verbose,
        parallel = parallel, n_cores = n_cores, cl = NULL,
        scores_weights = seq_len(stage1_npc), seed = seed
      )
      stage1_tau         <- tau_sel$tau
      stage1_tau_results <- tau_sel$tau_results
      stage1_fit         <- tau_sel$fit
    } else {
      stage1_tau <- tau[1]
      stage1_fit <- NULL
    }

    if (is.null(stage1_fit)) {
      if (verbose)
        cat(sprintf("=== Stage 1: Initial fit (nbasis=%d, npc=%d, lambda=1e-6, tau=%.4g) ===\n",
                    stage1_nbasis, stage1_npc, stage1_tau))
      stage1_fit <- fit_model(
        prep_data      = prep_data,
        nbasis_val     = stage1_nbasis,
        npc_val        = stage1_npc,
        lambda_val     = 1e-6,
        tau            = stage1_tau,
        basis_type     = basis_type,
        grid_size      = grid_size,
        solver         = solver,
        solver_control = solver_control,
        mu_nbasis      = mu_nbasis,
        scores_weights = seq_len(stage1_npc)
      )
    }

    if (is.null(stage1_fit)) {
      warning("Stage 1 fit failed; falling back to default linear weights")
      scores_weights <- seq_len(p_max)
    } else {
      eig_floor      <- pmax(stage1_fit$eigenvalues, 1e-6 * max(stage1_fit$eigenvalues))
      adaptive_weights <- max(eig_floor) / eig_floor
      scores_weights   <- pmin(adaptive_weights, 1e6)
    }

    if (share_tau && !is.null(stage1_tau))
      tau <- stage1_tau
  }

  selection_result <- select_parameters(
    data           = prep_data,
    nbasis_seq     = nbasis,
    npc_seq        = npc,
    lambda_seq     = lambda,
    tau            = tau,
    grid_size      = grid_size,
    basis_type     = basis_type,
    init_coef      = init_coef,
    solver         = solver,
    solver_control = solver_control,
    verbose        = verbose,
    mu_nbasis      = mu_nbasis,
    parallel       = parallel,
    n_cores        = n_cores,
    scores_weights = scores_weights,
    LRPsOCV_parameters = LRPsOCV_parameters,
    seed           = seed
  )

  optimal_fit <- selection_result$fit

  if (!identical(obs_range, c(0, 1))) {
    fxns           <- obs_range_reset(optimal_fit, basis_type, obs_range, grid_size, mu_nbasis)
    grid           <- fxns$grid
    eigenfunctions <- fxns$eigenfunctions
    mu_function    <- fxns$mu_fd
  } else {
    grid           <- optimal_fit$data$prep_data$grid
    eigenfunctions <- optimal_fit$eigenfunctions
    mu_function    <- optimal_fit$data$prep_data$mean_estimation$fdobj$fd
  }

  output <- list(
    U              = optimal_fit$U,
    Xi             = optimal_fit$Xi,
    xHat           = optimal_fit$xHat,
    Phi_mat        = optimal_fit$Phi_mat,
    eigenfunctions = eigenfunctions,
    eigenValues    = optimal_fit$eigenvalues,
    mu_function    = mu_function,
    scores_weights = scores_weights,
    sig2           = optimal_fit$sig2,
    npc            = ncol(optimal_fit$Phi_mat),
    nbasis         = nrow(optimal_fit$U),
    grid           = grid,
    AIC            = optimal_fit$AIC_K,
    AIC_tab        = selection_result$AIC_K_tab,
    tau            = selection_result$tau,
    lambda_selection   = selection_result$lambda_selection,
    nbasis_selection   = selection_result$nbasis_selection,
    grid_cv_results    = selection_result$grid_cv_results,
    lambda             = optimal_fit$lambda,
    tau_results        = selection_result$tau_results,
    stage1_tau         = stage1_tau,
    stage1_tau_results = stage1_tau_results,
    original_lambda_order = optimal_fit$ord,
    Fits = list(
      fitted_values = optimal_fit$fitted_values,
      residuals     = optimal_fit$residuals,
      var_explained = optimal_fit$var_explained
    ),
    Data            = optimal_fit$data,
    manOptim_output = optimal_fit$optimization_result,
    call_args = list(
      nbasis             = selection_result$nbasis_opt,
      npc                = selection_result$npc_opt,
      lambda             = selection_result$lambda_opt,
      lambda_original    = lambda_original,
      tau                = tau,
      grid_size          = grid_size,
      basis_type         = basis_type,
      solver             = solver,
      solver_control     = solver_control,
      mu_nbasis          = mu_nbasis,
      parallel           = parallel,
      n_cores            = n_cores_org,
      seed               = seed,
      LRPsOCV_parameters = LRPsOCV_parameters,
      scores_weights     = scores_weights,
      two_stage          = two_stage,
      share_tau          = share_tau,
      bootstrap_parameters = bootstrap_parameters
    )
  )

  xHat_CI <- NULL
  if (boot_compute) {
    if (verbose) cat("\n=== Computing Bootstrap Confidence Intervals ===\n")
    xHat_CI <- bootstrap_xHat_CI(
      optimal_fit    = optimal_fit,
      prep_data      = prep_data,
      nbasis_val     = selection_result$nbasis_opt,
      npc_val        = selection_result$npc_opt,
      lambda_val     = selection_result$lambda_opt,
      tau            = selection_result$tau,
      basis_type     = basis_type,
      grid_size      = grid_size,
      solver         = solver,
      solver_control = solver_control,
      mu_nbasis      = mu_nbasis,
      B              = boot_B,
      alpha          = boot_alpha,
      method         = boot_method,
      parallel       = parallel,
      n_cores        = n_cores,
      verbose        = verbose,
      scores_weights = scores_weights,
      seed = seed
    )
  }
  output$xHat_CI <- xHat_CI

  # --- Plots ---
  if (plot_results) {
    # 2-D grid CV results (nbasis×npc or npc×lambda)
    if (!is.null(selection_result$grid_cv_results)) {
      gcv       <- selection_result$grid_cv_results
      stat_cols <- c("mean_cv_error", "se_cv_error", "n_converged")
      p_cols    <- setdiff(names(gcv), stat_cols)
      if (length(p_cols) == 2) {
        opt_vals <- lapply(setNames(p_cols, p_cols), function(p)
          selection_result[[paste0(p, "_opt")]])
        plot_grid_cv_results(gcv, opt_vals)
      }
    }

    # AIC_K plot for npc selection
    if (nrow(selection_result$AIC_K_tab) > 1) {
      plot_AIC(
        AIC_tab    = selection_result$AIC_K_tab,
        nbasis_opt = selection_result$nbasis_opt,
        npc_opt    = selection_result$npc_opt,
        lambda_opt = selection_result$lambda_opt,
        verbose    = verbose
      )
    }

    # 1-D CV plots for lambda and nbasis
    if (!is.null(selection_result$lambda_selection)) {
      plot_cv_results(selection_result$lambda_selection, selection_result$lambda_opt)
    }
    if (!is.null(selection_result$nbasis_selection)) {
      plot_cv_results(selection_result$nbasis_selection, selection_result$nbasis_opt)
    }
  }

  if (verbose)
    cat(sprintf("Optimal: nbasis=%d, npc=%d, lambda=%g\n",
                selection_result$nbasis_opt,
                selection_result$npc_opt,
                selection_result$lambda_opt))

  gc()
  output
}
