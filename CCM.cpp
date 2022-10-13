#include <RcppArmadillo.h>
#include <omp.h>
#include <algorithm>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>

// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;
using namespace fill;

#ifdef _OPENMP
#include <omp.h>
#endif
arma::vec covarcpp(const arma::mat& XX,
                   const arma::colvec& yy){
  const int r = XX.n_rows;
  const int c = XX.n_cols;
  arma::mat XY(size(XX));
  arma::colvec F(c),FF(c),covar(c),fir(c);
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i=0;i<r;i++){
    XY.row(i) = XX.row(i) * yy(i);
  }
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i=0;i<c;i++){
    F(i)=accu(XY.col(i));
    FF(i) = mean(XX.col(i));
  }
  fir = FF*accu(yy);
  covar = (F - fir)/(r-1);
  return covar;
}
template<class T>
double var(T x,const bool std=false,const bool na_rm=false){
  int n = 0;
  double v=0.0,sum1=0,sum2=0;
  if(na_rm){
    for(auto v : x){
      if(!R_IsNA(v)){
        sum1+=v*v;
        sum2+=v;
        ++n;
      }
    }
  }else{
    n=x.size();
    double *xx=&x[0],*end=xx+n;
    for(;end-xx;++xx){
      v=*xx;
      sum1+=v*v;
      sum2+=v;
    }
  }
  v=(sum1-sum2*sum2/n)/(n-1);
  return std ? sqrt(v) : v;
}

colvec colVarscpp(mat x,const bool std=false,const bool na_rm=false,const bool parallel=false){
  colvec f(x.n_cols);
  if(parallel){
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(unsigned int i=0;i<x.n_cols;++i){
      f[i]=var<colvec>(x.col(i),std,na_rm);
    }
  }else{
    for(unsigned int i=0;i<x.n_cols;++i){
      f[i]=var<colvec>(x.col(i),std,na_rm);
    }
  }
  return f;
}

colvec row_means(mat x){
  mat X = mat(x.begin(), x.n_rows, x.n_cols, false);
  return mean(X, 1);
}

// [[Rcpp::export]]
Rcpp::List CCmadeM(  arma::mat& X_train,
                     arma::mat& X_test,
                    const arma::colvec& y_train,
                    const int chunk,
                    const int cc_num)
{ int N = y_train.n_rows;
  int M = X_test.n_rows;
  mat S1M(N,cc_num+1,fill::zeros);
  S1M.col(0)= vec(N,fill::ones);
  mat S2M(M,cc_num,fill::zeros);
  arma::vec U2,yt,r,sx,beta;
  arma::uvec id;
  arma::mat CV,Xtr,Xte;
  int i,s,bb,nn;
  if (chunk == 1){
    for ( i = 0; i < cc_num; ++i){
      CV = S1M.submat( 0, 0, N-1, i );
    //  U1 = CV.t() * y_train;
   //   U2 = solve(CV.t()*CV, U1);
      U2 = solve(CV, y_train);
      yt = y_train - CV * U2;
      r = covarcpp(X_train,yt);
      sx = colVarscpp(X_train,false,false,true);
      id = find( sx == 0 );
      if( id.n_elem > 0 ){
        sx.shed_rows(id);
        r.shed_rows(id);
        beta = r/sx;
        bb = beta.n_elem;
        X_train.shed_cols(id);
        X_test.shed_cols(id);
        S1M.col(i+1) = (X_train*beta)/bb;
        S2M.col(i) = (X_test*beta)/bb;
      }else{
        beta = r/sx;
        bb = beta.n_elem;
        S1M.col(i+1) = (X_train*beta)/bb;
        S2M.col(i) = (X_test*beta)/bb;
      }
    }
  }else{
    mat S1(N,chunk,fill::zeros);
    mat S2(M,chunk,fill::zeros);
    nn = X_train.n_cols;
    arma::colvec cc = linspace(0,nn,chunk+1);
    for ( i = 0; i < cc_num; ++i){
      CV = S1M.submat( 0, 0, N-1, i );
      for ( s = 0; s < chunk; ++s){
        Xtr = X_train.cols(cc(s),cc(s+1)-1);
        Xte = X_test.cols(cc(s),cc(s+1)-1);
        U2 = solve(CV, y_train);
        yt = y_train - CV * U2;
        r = covarcpp(Xtr,yt);
        sx = colVarscpp(Xtr,false,false,true);
        id = find( sx == 0 );
        if( id.n_elem > 0 ){
          sx.shed_rows(id);
          r.shed_rows(id);
          beta = r/sx;
          bb = beta.n_elem;
          Xtr.shed_cols(id);
          Xte.shed_cols(id);
          S1.col(s) = (Xtr*beta)/bb;
          S2.col(s) = (Xte*beta)/bb;
        }else{
          beta = r/sx;
          bb = beta.n_elem;
          S1.col(s) = (Xtr*beta)/bb;
          S2.col(s) = (Xte*beta)/bb;
        }
        S1M.col(i+1) = row_means(S1);
        S2M.col(i) = row_means(S2);
      }
    }
  }
  S1M.shed_col(0);
  return (List::create( Named("CCtrain")= S1M,
                        Named("CCtest")=S2M));
}

// [[Rcpp::export]]
Rcpp::List CCmadeBM(SEXP bX_train,
                    SEXP bX_test,
                    arma::colvec by_train,
                    int chunk,
                    int cc_num)
  {
  XPtr<BigMatrix> bmX_train(bX_train);
  XPtr<BigMatrix> bmX_test(bX_test);
  arma::Mat<double> Xtrain ((double *)bmX_train->matrix(), bmX_train->nrow(), bmX_train->ncol(), false);
  arma::Mat<double> Xtest ((double *)bmX_test->matrix(), bmX_test->nrow(), bmX_test->ncol(), false);
  Rcpp::List CC = CCmadeM(Xtrain,Xtest,by_train,chunk,cc_num);
  return CC;
  }
