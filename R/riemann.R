# Riemannian optimization on the Stiefel x Euclidean product manifold
Riemann <- function(
    data, method = "RCG",
    mani.params = get.manifold.params(),
    solver.params = get.solver.params(Max_Iteration = 150000, Tolerance = 1e-04, DEBUG = 0)
) {
  prep_data <- data$prep_data
  basis_data <- data$basis_data
  orthB_list <- basis_data$orthB_list
  y_vec_list <- basis_data$y_vec_list
  mi_vec <- lengths(y_vec_list)
  data_list <- prep_data$data_list
  orthP <- basis_data$orthP
  Gram_mat <- basis_data$Gram_mat
  orthB_grid <- basis_data$orthB_grid
  mu_fd <- prep_data$mu_fxn
  grid <- prep_data$grid
  N <- prep_data$n
  D <- data$D
  K <- data$K
  G <- prep_data$G
  lambda <- data$lambda
  tau <- data$tau
  init_coef <- data$init_coef

  mani.defn <- get.product.defn(
    get.stiefel.defn(n = D, p = K),
    get.euclidean.defn(n = N, m = K)
  )

  stiefel_dim <- D * K
  total_dim <- stiefel_dim + N * K

  mod <- Rcpp::Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
  prob <- new(
    mod$RProblem,
    function(x) objective_manifold(x, D, K, N, stiefel_dim, total_dim,
                                   orthB_list, y_vec_list, orthP, lambda, tau),
    function(x) gradient_manifold(x, D, K, N, stiefel_dim, total_dim,
                                  orthB_list, y_vec_list, orthP, lambda, tau)
  )

  if (!is.null(init_coef)) {
    U_init <- init_coef$U
    Xi_init <- init_coef$Xi
  } else {
    set.seed(123)
    U_init <- qr.Q(qr(matrix(runif(D * K, -2, 2), nrow = D)))
    Xi_init <- matrix(0, nrow = N, ncol = K)
    for (k in 1:K) {
      set.seed(123)
      Xi_init[, k] <- rnorm(N, sd = 0.1 / sqrt(k))
    }
  }

  ord_eig_init <- eig_order(U_init, Xi_init)
  U_init <- ord_eig_init$U
  Xi_init <- ord_eig_init$Xi
  x0 <- c(as.vector(U_init), as.vector(Xi_init))

  result <- manifold.optim(
    prob, mani.defn, method = method,
    mani.params = mani.params, solver.params = solver.params, x0 = x0
  )

  U_final <- matrix(result$xopt[1:stiefel_dim], nrow = D, ncol = K)
  Xi_final <- matrix(result$xopt[(stiefel_dim + 1):total_dim], nrow = N, ncol = K)

  orthonormality_check <- max(abs(t(U_final) %*% U_final - diag(K)))

  U_original <- Gram_mat %*% U_final
  eigenvalues <- apply(Xi_final, 2, var)
  ord <- order(eigenvalues, decreasing = TRUE)
  eigenvalues <- eigenvalues[ord]
  U_original <- U_original[, ord, drop = FALSE]
  U_final <- U_final[, ord, drop = FALSE]
  Xi_final <- Xi_final[, ord, drop = FALSE]

  eigenfunctions <- lapply(1:K, function(k) fd(U_original[, k], basis_data$basis))
  Phi_mat <- orthB_grid %*% U_final

  fitted_values <- list()
  residuals <- list()
  xHat <- matrix(NA, N, G)

  for (i in 1:N) {
    fitted_values[[i]] <- orthB_list[[i]] %*% U_final %*% Xi_final[i, ]
    residuals[[i]] <- y_vec_list[[i]] - fitted_values[[i]]
    xHat[i, ] <- orthB_grid %*% U_final %*% Xi_final[i, ] + eval.fd(grid, mu_fd)
  }

  NN <- sum(mi_vec)
  sig2 <- mean(sapply(residuals, function(r) crossprod(r) / length(r)))

  B_grid <- basis_data$B_grid
  S <- t(B_grid) %*% B_grid
  DFu <- sum(diag(S %*% solve(S + lambda * orthP))) * K
  aDF <- DFu + sum(N / (1 + (0:(K - 1))^2 * tau))

  AIC <- NN * log(result$fval / NN) + aDF
  AIC_K <- NN * log(sig2) + NN + 2 * N * K

  total_var <- sum(sapply(y_vec_list, var))
  var_explained <- cumsum(eigenvalues) / total_var

  list(
    solver = method, U = U_final, U_original = U_original,
    Xi = Xi_final, eigenvalues = eigenvalues, eigenfunctions = eigenfunctions,
    Phi_mat = Phi_mat, sig2 = sig2, AIC = AIC, AIC_K = AIC_K,
    aDF = aDF, ord = ord, fitted_values = fitted_values, residuals = residuals,
    xHat = xHat, var_explained = var_explained, data = data,
    optimization_result = result, fval = result$fval, iter = result$iter,
    convergence = (result$normgfgf0 < solver.params$Tolerance),
    orthonormality_error = orthonormality_check, lambda = lambda
  )
}
