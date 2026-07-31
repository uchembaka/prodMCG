
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>



// [[Rcpp::export]]
double objective_manifold(const arma::vec& x,
                          const int D, const int K, const int N,
                          const int stiefel_dim, const int total_dim,
                          const Rcpp::List& orthB_list,
                          const Rcpp::List& y_vec_list,
                          const arma::mat& orthP,
                          const double lambda,
                          const double tau,
                          const arma::vec& ridge_weights) {

  if ((int)x.n_elem != total_dim) Rcpp::stop("x length != total_dim");

  arma::mat U(const_cast<double*>(x.memptr()), D, K, false, true);
  const arma::vec x_tail = x.subvec(stiefel_dim, total_dim - 1);
  arma::mat Xi(const_cast<double*>(x_tail.memptr()), N, K, false, true);

  double obj = 0.0;

  // Data fitting term
  for (int i = 0; i < N; ++i) {
    const arma::mat Bi = Rcpp::as<arma::mat>(orthB_list[i]);
    const arma::vec yi = Rcpp::as<arma::vec>(y_vec_list[i]);
    arma::vec pred = Bi * U * Xi.row(i).t();
    arma::vec r = yi - pred;
    obj += 0.5 * arma::dot(r, r);
  }

  // Roughness penalty: 0.5 * lambda * tr(U' P U)
  obj += 0.5 * lambda * arma::trace(U.t() * orthP * U);

  // Component-specific ridge penalty: 0.5 * tau * sum_k w_k * ||Xi_k||^2
  double penalty = 0.0;
  for (int i = 0; i < N; ++i) {
    for (int k = 0; k < K; ++k) {
      penalty += ridge_weights(k) * Xi(i, k) * Xi(i, k);
    }
  }
  obj += 0.5 * tau * penalty;

  return obj;
}

// [[Rcpp::export]]
arma::vec gradient_manifold(const arma::vec& x,
                            const int D, const int K, const int N,
                            const int stiefel_dim, const int total_dim,
                            const Rcpp::List& orthB_list,
                            const Rcpp::List& y_vec_list,
                            const arma::mat& orthP,
                            const double lambda,
                            const double tau,
                            const arma::vec& ridge_weights) {

  if ((int)x.n_elem != total_dim) Rcpp::stop("x length != total_dim");

  arma::mat U(const_cast<double*>(x.memptr()), D, K, false, true);
  const arma::vec x_tail = x.subvec(stiefel_dim, total_dim - 1);
  arma::mat Xi(const_cast<double*>(x_tail.memptr()), N, K, false, true);

  arma::mat grad_U(D, K, arma::fill::zeros);
  arma::mat grad_Xi(N, K, arma::fill::zeros);

  // Data fitting gradient
  for (int i = 0; i < N; ++i) {
    const arma::mat Bi = Rcpp::as<arma::mat>(orthB_list[i]);
    const arma::vec yi = Rcpp::as<arma::vec>(y_vec_list[i]);
    arma::vec pred = Bi * U * Xi.row(i).t();
    arma::vec r = yi - pred;

    arma::vec gxi = -(U.t() * Bi.t() * r);
    grad_Xi.row(i) = gxi.t();

    grad_U -= (Bi.t() * r) * Xi.row(i);
  }

  // Roughness penalty gradient
  grad_U += lambda * (orthP * U);

  // Ridge penalty gradient: tau * w_k * Xi_k
  for (int k = 0; k < K; ++k) {
    grad_Xi.col(k) += tau * ridge_weights(k) * Xi.col(k);
  }

  // Concatenate gradients in column-major order
  arma::vec out(total_dim);
  out.subvec(0, stiefel_dim - 1) = arma::vectorise(grad_U);
  out.subvec(stiefel_dim, total_dim - 1) = arma::vectorise(grad_Xi);
  return out;
}


