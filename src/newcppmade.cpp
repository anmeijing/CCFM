#include <RcppArmadillo.h>
#include <omp.h>
#include <bigmemory/BigMatrix.h>

// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

// Function to compute covariance in C++
vec covarcpp(const arma::mat& XX, const arma::colvec& yy) {
  const int r = XX.n_rows;
  const int c = XX.n_cols;
  arma::mat XY(r, c);
  arma::colvec F(c), FF(c), covar(c), fir(c);

  if (r != yy.n_rows) {
    stop("The number of rows in XX must be equal to the number of rows in yy.");
  }

#pragma omp parallel for
  for(int i = 0; i < r; i++) {
    XY.row(i) = XX.row(i) * yy(i);
  }

#pragma omp parallel for
  for(int i = 0; i < c; i++) {
    F(i) = accu(XY.col(i));
    FF(i) = mean(XX.col(i));
  }

  fir = FF * accu(yy);
  covar = (F - fir) / (r - 1);
  return covar;
}

// Function to compute variance in C++
template<class T>
double varcpp(T x, const bool std = false, const bool na_rm = false) {
  int n = 0;
  double v = 0.0, sum1 = 0.0, sum2 = 0.0;
  for(auto val : x) {
    if (na_rm && R_IsNA(val)) {
      continue;
    } else {
      sum1 += val * val;
      sum2 += val;
      ++n;
    }
  }
  v = (n > 1) ? (sum1 - sum2 * sum2 / n) / (n - 1) : 0.0;
  return std ? std::sqrt(v) : v;
}

// Function to compute column variances
arma::colvec colVarscpp(const arma::mat& x, const bool std = false, const bool na_rm = false, const bool parallel = false) {
  arma::colvec f(x.n_cols);
  if(parallel) {
#pragma omp parallel for
    for(unsigned int i = 0; i < x.n_cols; ++i) {
      f[i] = varcpp(x.col(i), std, na_rm);
    }
  } else {
    for(unsigned int i = 0; i < x.n_cols; ++i) {
      f[i] = varcpp(x.col(i), std, na_rm);
    }
  }
  return f;
}

// Function to compute row means
arma::colvec row_means(const arma::mat& x) {
  return mean(x, 1);
}

// Main function to perform CC calculation
// [[Rcpp::export]]
Rcpp::List CCmadeM(arma::mat X_train, arma::mat X_test, const arma::colvec& y_train, int chunk, int cc_num) {
  int N = y_train.n_rows;
  int M = X_test.n_rows;
  arma::mat S1M(N, cc_num + 1, fill::zeros);
  S1M.col(0) = ones(N, 1);
  arma::mat S2M(M, cc_num, fill::zeros);
  vec U2, yt, r, sx, beta;
  uvec id;
  arma::mat CV, Xtr, Xte;
  int i, s, bb, nn;

  if (chunk == 1) {
    for (i = 0; i < cc_num; ++i) {
      CV = S1M.submat(0, 0, N - 1, i);
      U2 = solve(CV, y_train);
      yt = y_train - CV * U2;
      r = covarcpp(X_train, yt);
      sx = colVarscpp(X_train, false, false, true);
      id = find(sx == 0);
      if(id.n_elem > 0) {
        sx.shed_rows(id);
        r.shed_rows(id);
        beta = r / sx;
        bb = beta.n_elem;
        X_train.shed_cols(id);
        X_test.shed_cols(id);
        S1M.col(i + 1) = (X_train * beta) / bb;
        S2M.col(i) = (X_test * beta) / bb;
      } else {
        beta = r / sx;
        bb = beta.n_elem;
        S1M.col(i + 1) = (X_train * beta) / bb;
        S2M.col(i) = (X_test * beta) / bb;
      }
    }
  } else {
    for (i = 0; i < cc_num; ++i) {
      CV = S1M.submat(0, 0, N - 1, i);
      for (s = 0; s < chunk; ++s) {
        Xtr = X_train.cols(s * X_train.n_cols / chunk, (s + 1) * X_train.n_cols / chunk - 1);
        Xte = X_test.cols(s * X_test.n_cols / chunk, (s + 1) * X_test.n_cols / chunk - 1);
        U2 = solve(CV, y_train);
        yt = y_train - CV * U2;
        r = covarcpp(Xtr, yt);
        sx = colVarscpp(Xtr, false, false, true);
        id = find(sx == 0);
        if(id.n_elem > 0) {
          sx.shed_rows(id);
          r.shed_rows(id);
          beta = r / sx;
          bb = beta.n_elem;
          Xtr.shed_cols(id);
          Xte.shed_cols(id);
          S1M.col(i + 1) += (Xtr * beta) / bb;
          S2M.col(i) += (Xte * beta) / bb;
        } else {
          beta = r / sx;
          bb = beta.n_elem;
          S1M.col(i + 1) += (Xtr * beta) / bb;
          S2M.col(i) += (Xte * beta) / bb;
        }
      }
      S1M.col(i + 1) /= chunk;
      S2M.col(i) /= chunk;
    }
  }
  S1M.shed_col(0);
  return (List::create(Named("CCtrain") = S1M, Named("CCtest") = S2M));
}

// Function to handle bigmemory matrices
// [[Rcpp::export]]
Rcpp::List CCmadeBM(SEXP bX_train, SEXP bX_test, arma::colvec by_train, int chunk, int cc_num) {
  XPtr<BigMatrix> bmX_train(bX_train);
  XPtr<BigMatrix> bmX_test(bX_test);
  Mat<double> Xtrain((double *)bmX_train->matrix(), bmX_train->nrow(), bmX_train->ncol(), false);
  Mat<double> Xtest((double *)bmX_test->matrix(), bmX_test->nrow(), bmX_test->ncol(), false);
  Rcpp::List CC = CCmadeM(Xtrain, Xtest, by_train, chunk, cc_num);
  return CC;
}
