# Main user-facing wrapper for sparse functional PCA

#' Reconstruct sparse functional data via Riemannian FPCA
#'
#' Main user-facing function that performs functional principal component
#' analysis and curve reconstruction of sparse longitudinal data
#'  using Riemannian optimization on the product manifold
#' Stiefel(D,K) x Euclidean(N,K).
#'
#' @param data Matrix or data frame with 3 columns: subject ID, observation time, observed value.
#' @param nbasis Integer or integer vector. Number of basis functions. If a vector, optimal
#'   nbasis is selected via cross-validation (lambda forced to 1e-6).
#' @param npc Integer or integer vector. Number of principal components. If a vector,
#'   optimal npc is selected via AIC.
#' @param grid_size Integer. Number of grid points for evaluation (default 101).
#' @param basis_type Character. Type of basis: "bspline" (default) or "fourier".
#' @param lambda Numeric or numeric vector. Smoothing penalty. If a vector, optimal lambda
#'   is selected via AIC (default) or cross-validation.
#' @param tau Numeric or numeric vector. Ridge penalty for component ordering.
#'   If a vector, optimal tau is selected using the method specified by
#'   \code{tau_selection_method}.
#' @param solver Character. Optimization method: "CG" (default), "BFGS", or "L-BFGS-B".
#' @param solver_control List with solver control parameters:
#'   \itemize{
#'     \item maxit: Maximum iterations (default 3e+5)
#'     \item trace: Debug level (default 0)
#'     \item Tolerance: Convergence tolerance (default 1e-4)
#'   }
#' @param obs_range Numeric vector of length 2. Observation time range (default c(0,1)).
#' @param init_coef NULL or list with components U (D x K matrix) and Xi (N x K matrix).
#'   Initial coefficient values for optimization.
#' @param plot_results Logical. Plot AIC and CV results (default FALSE).
#' @param verbose Logical. Print progress messages (default FALSE).
#' @param mu_nbasis Integer. Number of basis functions for mean estimation (default 15).
#' @param parallel Logical. Use parallel computation for CV (default FALSE).
#' @param n_cores Integer or NULL. Number of cores for parallel computation.
#'   NULL uses detectCores() - 1.
#' @param cv_nfolds Integer. Number of folds for cross-validation (default 5).
#'   Must be a positive integer >= 2.
#' @param tau_selection_method Character. Method for selecting tau when \code{tau}
#'   is a vector. One of:
#'   \itemize{
#'     \item \code{"eigen_order"} (default): Searches for smallest tau that yields
#'       correct eigenvalue ordering (two consecutive correct orderings). When
#'       \code{parallel = TRUE}, fits ALL tau candidates in parallel, then applies
#'       the selection logic on collected results. When \code{parallel = FALSE},
#'       uses sequential search with early stopping.
#'     \item \code{"aic"}: Fit all tau candidates in parallel and select the one
#'       with minimum AIC. Always uses parallel computation.
#'     \item \code{"cv"}: k-fold cross-validation over tau candidates with 1-SE
#'       rule. Always uses parallel computation.
#'   }
#' @param selection_method Character. Method for selecting lambda and nbasis when
#'   multiple candidates are provided. One of:
#'   \itemize{
#'     \item \code{"cv"} (default): k-fold cross-validation with 1-SE rule.
#'     \item \code{"aic"}: AIC-based selection (minimum AIC).
#'   }
#' @param bootstrap_parameters List controlling bootstrap confidence intervals
#'   for the reconstructed trajectories (xHat). Elements:
#'   \itemize{
#'     \item \code{compute}: Logical. Whether to compute bootstrap CIs (default FALSE).
#'     \item \code{B}: Integer. Number of bootstrap replicates (default 200).
#'     \item \code{alpha}: Numeric. Significance level for (1-alpha)*100\% CIs (default 0.05).
#'     \item \code{method}: Character. Bootstrap method: \code{"parametric"} (default)
#'       keeps original scores fixed and samples fresh errors from N(0, sig2),
#'       or \code{"residual"} resamples residuals from the original fit.
#'   }
#'   Model parameters are fixed from the original fit for both methods.
#'
#' @return A list with components:
#'   \describe{
#'     \item{U}{Orthonormal basis coefficient matrix (D x K)}
#'     \item{Xi}{Score matrix (N x K)}
#'     \item{xHat}{Reconstructed functions (N x G matrix)}
#'     \item{Phi_mat}{Eigenfunction evaluations on grid (G x K)}
#'     \item{eigenfunctions}{List of K eigenfunction fd objects}
#'     \item{eigenValues}{Eigenvalues (variances of scores)}
#'     \item{mu_function}{Mean function (fd object)}
#'     \item{sig2}{Residual variance}
#'     \item{npc}{Selected number of principal components}
#'     \item{nbasis}{Selected number of basis functions}
#'     \item{grid}{Evaluation grid}
#'     \item{AIC}{AIC value of final model}
#'     \item{AIC_tab}{Data frame of AIC values across npc candidates}
#'     \item{AIC_tab_param}{Data frame of AIC values for lambda/nbasis selection
#'       (when \code{selection_method = "aic"}). NULL if CV was used.}
#'     \item{tau}{Selected/used tau value}
#'     \item{lambda}{Selected/used lambda value}
#'     \item{CV_results}{Cross-validation results (if applicable)}
#'     \item{tau_results}{Tau selection results table (if tau was a vector).
#'       For \code{"aic"}: data frame with tau, AIC, AIC_K, fval, aDF,
#'       convergence, iter, correct_order, delta_AIC.
#'       For \code{"cv"}: data frame with tau, mean_cv_error, se_cv_error,
#'       n_converged.
#'       For \code{"eigen_order"}: data frame with tau, AIC, AIC_K, fval,
#'       aDF, convergence, iter, correct_order.}
#'     \item{Fits}{List with fitted_values, residuals, var_explained}
#'     \item{Data}{Full data object from fit}
#'     \item{manOptim_output}{Raw ManifoldOptim output}
#'     \item{xHat_CI}{List of N matrices (each G x 2 with columns "lower" and
#'       "upper") for pointwise bootstrap confidence intervals. NULL if bootstrap
#'       not computed.}
#'     \item{call_args}{List of all arguments used}
#'   }
#'
#' @export
reconstructSparseFuncs <- function(
    data,
    nbasis = 10,
    npc = 2:5,
    grid_size = 101,
    basis_type = "bspline",
    lambda = exp(seq(-3, 2, length.out = 10)),
    tau = 10^(-6:1),
    solver = "CG",
    solver_control = list(maxit = 300000, trace = 1, Tolerance = 1e-04),
    obs_range = c(0,1),
    init_coef = NULL,
    plot_results = FALSE,
    verbose = FALSE,
    mu_nbasis = 15,
    parallel = TRUE,
    n_cores = NULL,
    cv_nfolds = 5,
    tau_selection_method = "eigen_order",
    selection_method = 'aic',
    bootstrap_parameters = list(compute = FALSE, B = 100, alpha = 0.05, method = 'residual')
) {

  solver     <- match.arg(solver, c("BFGS", "L-BFGS-B", "CG"))
  basis_type <- match.arg(basis_type, c("bspline", "fourier"))

  if (!is.numeric(cv_nfolds) || length(cv_nfolds) != 1 ||
      cv_nfolds != as.integer(cv_nfolds) || cv_nfolds < 2) {
    stop("cv_nfolds must be a positive integer >= 2")
  }
  cv_nfolds <- as.integer(cv_nfolds)

  tau_selection_method <- match.arg(tau_selection_method,
                                     c("eigen_order", "aic", "cv"))

  selection_method <- match.arg(selection_method, c("cv", "aic"))

  if (!is.list(bootstrap_parameters)) {
    stop("bootstrap_parameters must be a list with elements: compute, B, alpha, method")
  }
  boot_compute <- isTRUE(bootstrap_parameters$compute)
  boot_B       <- bootstrap_parameters$B %||% 200
  boot_alpha   <- bootstrap_parameters$alpha %||% 0.05
  boot_method  <- bootstrap_parameters$method %||% "parametric"

  if (!is.numeric(boot_B) || length(boot_B) != 1 || boot_B < 2) {
    stop("bootstrap_parameters$B must be a positive integer >= 2")
  }
  boot_B <- as.integer(boot_B)

  if (!is.numeric(boot_alpha) || length(boot_alpha) != 1 ||
      boot_alpha <= 0 || boot_alpha >= 1) {
    stop("bootstrap_parameters$alpha must be a numeric value in (0, 1)")
  }

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
    if (!is.list(init_coef) || is.null(init_coef$U) || is.null(init_coef$Xi)) {
      stop("init_coef must be NULL or a list with components $U and $Xi")
    }
  }

  prep_data <- preprocess_data(
    data       = data,
    basis_type = basis_type,
    G          = grid_size,
    obs_range  = obs_range,
    mu_nbasis  = mu_nbasis
  )

  lambda_original <- lambda

  if (verbose && (length(nbasis) > 1 || length(npc) > 1 || length(lambda) > 1)) {
    cat("Starting parameter selection...\n")
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
    n_cores              = n_cores,
    cv_nfolds            = cv_nfolds,
    tau_selection_method = tau_selection_method,
    selection_method     = selection_method
  )

  optimal_fit <- selection_result$fit

  if (!identical(obs_range, c(0, 1))) {
    fxns <- obs_range_reset(
      result     = optimal_fit,
      basis_type = basis_type,
      obs_range  = obs_range,
      G          = grid_size,
      mu_nbasis  = mu_nbasis
    )
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
    sig2           = optimal_fit$sig2,
    npc            = ncol(optimal_fit$Phi_mat),
    nbasis         = nrow(optimal_fit$U),
    grid           = grid,
    AIC            = optimal_fit$AIC,
    AIC_tab        = selection_result$AIC_tab,
    AIC_tab_param  = selection_result$AIC_tab_param,
    tau            = selection_result$tau,
    lambda         = optimal_fit$lambda,
    CV_results     = selection_result$cv_results,
    tau_results    = selection_result$tau_results,
    original_lambda_order = optimal_fit$ord,
    Fits = list(
      fitted_values = optimal_fit$fitted_values,
      residuals     = optimal_fit$residuals,
      var_explained = optimal_fit$var_explained
    ),
    Data             = optimal_fit$data,
    manOptim_output  = optimal_fit$optimization_result,
    call_args = list(
      nbasis          = selection_result$nbasis_opt,
      npc             = selection_result$npc_opt,
      lambda          = selection_result$lambda_opt,
      lambda_original = lambda_original,
      tau             = tau,
      grid_size       = grid_size,
      basis_type      = basis_type,
      solver          = solver,
      solver_controls = solver_control,
      mu_nbasis       = mu_nbasis,
      parallel             = parallel,
      cv_nfolds            = cv_nfolds,
      tau_selection_method  = tau_selection_method,
      selection_method      = selection_method,
      bootstrap_parameters  = bootstrap_parameters
    )
  )

  xHat_CI <- NULL
  if (boot_compute) {
    if (verbose) cat("\n=== Computing Bootstrap Confidence Intervals ===\n")
    # if (verbose) cat(sprintf("method=%s, B=%d, alpha=%.2f, parallel=%s\n",
    #                          boot_method, boot_B, boot_alpha, parallel))

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
      verbose        = verbose
    )
  }
  output$xHat_CI <- xHat_CI

  if (plot_results && nrow(selection_result$AIC_tab) > 1) {
    plot_AIC(
      AIC_tab    = selection_result$AIC_tab,
      nbasis_opt = selection_result$nbasis_opt,
      npc_opt    = selection_result$npc_opt,
      lambda_opt = selection_result$lambda_opt,
      verbose    = verbose
    )
  }
  if (plot_results && !is.null(selection_result$AIC_tab_param)) {
    plot_AIC(
      AIC_tab    = selection_result$AIC_tab_param,
      nbasis_opt = selection_result$nbasis_opt,
      npc_opt    = selection_result$npc_opt,
      lambda_opt = selection_result$lambda_opt,
      verbose    = verbose
    )
  }
  if (plot_results && !is.null(selection_result$cv_results)) {
    cv_res <- selection_result$cv_results
    if ("lambda" %in% names(cv_res)) {
      plot_cv_results(cv_res, selection_result$lambda_opt)
    } else if ("nbasis" %in% names(cv_res)) {
      plot_cv_results(cv_res, selection_result$nbasis_opt)
    }
  }

  if (verbose) {
    cat(sprintf("Optimal: nbasis=%d, npc=%d, lambda=%g\n",
                selection_result$nbasis_opt,
                selection_result$npc_opt,
                selection_result$lambda_opt))
    if (nrow(selection_result$AIC_tab) > 1) {
      cat(sprintf("Minimum AIC: %.2f\n", min(selection_result$AIC_tab$AIC, na.rm = TRUE)))
    }
  }

  gc()

  output
}
