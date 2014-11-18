mvrNormfun <- '#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma, int seed) {
   int ncols = sigma.n_cols;
   arma_rng::set_seed(seed);
   arma::mat Y = arma::randn(n, ncols);
   return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}'
sigma <- matrix(c(1, 0.9, -0.3, 0.9, 1, -0.4, -0.3, -0.4, 1), ncol = 3)
mu <- c(10, 5, -3)
cor(mvrnormArma(100, mu, sigma, 456))