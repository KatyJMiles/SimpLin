#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List SimpLinCpp(arma::vec x, arma::vec y) {
  mat xmat = join_rows(ones(x.size(), 1), x);
  colvec beta = solve(xmat, y);
  colvec pred_vals = xmat*beta;
  colvec resid = y - pred_vals;
  colvec se = sqrt(as_scalar(arma::trans(resid)*resid / (y.size() - 2)) * diagvec(inv(trans(xmat)*xmat)));
  colvec ci_upper = beta + se*R::qnorm(.975, 0, 1, 1, 0);
  colvec ci_lower = beta - se*R::qnorm(.975, 0, 1, 1, 0);
  mat ci_mat = join_rows(ci_lower, ci_upper);
  return List::create(beta, se, ci_upper, ci_lower, resid, pred_vals);
}
