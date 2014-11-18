library(Rcpp)
library(RcppArmadillo)
sourceCpp( code = '
// [[Rcpp::depends("RcppArmadillo")]]

#include <RcppArmadillo.h>

using namespace arma ; 
using namespace Rcpp ;

// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma, int seed) {
		int ncols = sigma.n_cols;
		arma_rng::set_seed(seed);
		arma::mat Y = arma::randn(n, ncols);
		return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
		}
')
