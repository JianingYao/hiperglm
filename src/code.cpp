#include <Rcpp.h>
#include "hiperglm_types.h"
using namespace Rcpp;

// [[Rcpp::export]]
VectorXd qr_solve_rcpp(const MatrixXd& A, const VectorXd& y) {
  if (A.rows() != y.size()) {
    Rcpp::stop("Incompatible matrix dimensions.");
  }
  HouseholderQR<MatrixXd> qr(A);
  qr = qr.compute(A);
  VectorXd beta_hat = qr.solve(y);
  return beta_hat;
}
