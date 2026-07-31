#' Compute bootstrap confidence intervals for reconstructed functions
#'
#' Generates pointwise bootstrap confidence intervals for the estimated curves
#' (\eqn{\hat{x}_i(t)}) from a fitted \code{reconstructSparseFuncs} object.
#'
#' @param fit_obj A fitted object returned by \code{\link{reconstructSparseFuncs}}.
#' @param B Integer. Number of bootstrap replicates (default 200).
#' @param alpha Numeric. Significance level for (1-alpha)*100\% CIs (default 0.05).
#' @param method Character. Bootstrap method: \code{"parametric"} (default) or
#'   \code{"residual"}.
#' @param parallel Logical. Use parallel computation (default TRUE).
#' @param n_cores Integer or NULL. Number of cores; NULL uses detectCores() - 1.
#' @param verbose Logical. Print progress messages (default FALSE).
#'
#' @return A list of N matrices, each G x 2 with columns \code{"lower"} and
#'   \code{"upper"} for pointwise bootstrap confidence intervals.
#'
#' @seealso \code{\link{reconstructSparseFuncs}}
#' @export
compute_bootstrap <- function(
    fit_obj,
    B = 200,
    alpha = 0.05,
    method = "parametric",
    parallel = TRUE,
    n_cores = max(1, parallel::detectCores() - 1),
    verbose = FALSE
) {
  optimal_fit <- list(
    U              = fit_obj$U,
    Xi             = fit_obj$Xi,
    xHat           = fit_obj$xHat,          # named element (was unnamed before)
    fitted_values  = fit_obj$Fits$fitted_values,
    residuals      = fit_obj$Fits$residuals,
    sig2           = fit_obj$sig2
  )
  prep_data <- fit_obj$Data$prep_data

  Xhat_CI <- bootstrap_xHat_CI(
    optimal_fit    = optimal_fit,
    prep_data      = prep_data,
    nbasis_val     = fit_obj$nbasis,
    npc_val        = fit_obj$npc,
    lambda_val     = fit_obj$lambda,
    tau            = fit_obj$tau,
    basis_type     = fit_obj$call_args$basis_type,
    grid_size      = fit_obj$call_args$grid_size,
    solver         = fit_obj$call_args$solver,
    solver_control = fit_obj$call_args$solver_control,  # fixed: was solver_controls
    mu_nbasis      = fit_obj$call_args$mu_nbasis,
    B              = B,
    alpha          = alpha,
    method         = method,
    parallel       = parallel,
    n_cores        = n_cores,
    verbose        = verbose,
    scores_weights = fit_obj$scores_weights,
    seed           = fit_obj$seed
  )

  return(Xhat_CI)
}
