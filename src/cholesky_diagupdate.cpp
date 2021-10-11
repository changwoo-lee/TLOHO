//' @name cholesky_diagonalcpp
//' @title Cholesky diagonal update
//' 
//' @description It calculates the Cholesky factor of (A + diag(diagadd)) directly from chol(A) and diagadd vector. See Algorithm 2 in Appendix A1 of Lee et al.(2021)
//' 
//' @param R n by n matrix, right triangular Cholesky factor R, where t(R)*R = A 
//' @param diagadd n by 1 vector, diagonal adding vector 
//' @export
//' @references Lee, C. J., Luo, Z. T., & Sang, H. (2021). T-LoHo: A Bayesian Regularization Model for Structured Sparsity and Smoothness on Graphs. \emph{Advances in Neural Information Processing Systems 34 (NeurIPS 2021)}.
//' 
//' @examples
//' 
//' X = matrix(rnorm(1000*200),1000, 200)
//' A = crossprod(X) # positive definite A
//' R = chol(A) 
//' diagadd = runif(200)
//' 
//' Rnew = cholesky_diagonalcpp(R, diagadd) 
//' all.equal(Rnew, chol(crossprod(X)+diag(diagadd)))
//'
#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
arma::mat cholesky_diagonalcpp(arma::mat& R, arma::vec& diagadd){
  int n = R.n_rows;
  if (n == 1){
    arma::mat Rnew = R;
    Rnew = sqrt(square(Rnew) + diagadd);
    return Rnew;
  }
  arma::mat Rnew = R;
  arma::vec Lii = R.diag();
  arma::vec diag_original = trans(sum(square(R)));
  arma::vec dnew = diag_original + diagadd;
  for (int i = 0; i < n; i++){
    if (i == 0){
      Rnew(0, 0) = sqrt(dnew(0));
      Rnew(i, span(1,n-1)) = R(i, span(1,n-1)) * Lii(0)/Rnew(0,0);
    } else if(i < n-1){
      Rnew(i, i) = sqrt(dnew(i)-sum( square(Rnew(span(0,i-1), i)) ) );
      Rnew(i, span(i+1,n-1)) = (R(i, span(i+1,n-1)) * Lii(i) + trans(R(span(0,i-1), i))*R(span(0,i-1),span(i+1,n-1)) - trans(Rnew(span(0,i-1), i))*Rnew(span(0,i-1),span(i+1,n-1)) )/Rnew(i,i);
    }else{
      Rnew(n-1,n-1) = sqrt(dnew(n-1)-sum( square(Rnew(span(0,n-2), n-1)) ) );
    }
  }
  return Rnew;
}


