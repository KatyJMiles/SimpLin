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
  colvec ci_upper = beta + se*R::qt(.975, x.size() - 1, 1, 0);
  colvec ci_lower = beta - se*R::qt(.975, x.size() - 1, 1, 0);
  mat ci_mat = join_rows(ci_lower, ci_upper);
  return List::create(Named("coefficients") = beta, 
                      Named("se") = se, 
                      Named("CI_upper") = ci_upper,
                      Named("CI_lower") = ci_lower, 
                      Named("Residuals") = resid, 
                      Named("Fitted_Values") = pred_vals);
}
