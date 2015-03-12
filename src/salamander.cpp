/*------------------- created 03/10/2014 by Erlis Ruli ------------------------
#################### iLaplace for salamander data  ############################
-----------------------------------------------------------------------------*/

#include <RcppArmadillo.h>
using namespace Rcpp;


double ouflow(double xb) {
  return -log(exp(0.0) + exp(xb));
}

arma::vec vouflow(arma::vec xb) {
  int n = xb.n_elem;
  arma::vec ans(n);
  for(int i=0; i<n; i++) {
    ans(i) = ouflow(xb(i));
  }
  return ans;
}


//[[Rcpp::export]]
double nlogH(arma::vec randPar, arma::vec fixPar, List data)
{
  //arma::vec y = data["y"];
  arma::vec y = data["y"];
  arma::mat X = data["x"];
  arma::mat Z = data["z"];
  int nData = Z.n_rows;
  arma::vec theta = fixPar(arma::span(0, 3));
  double sigA = fixPar(4), sigB = fixPar(5);
  
  arma::vec a = randPar(arma::span(0, 19));
  arma::vec b = randPar(arma::span(20, 39));
  
  double prob = 0.0, eta = 0.0, ll = 0.0, ans = 0.0;
  for(int i = 0; i< nData; i++) {
    eta = sum(X.row(i)*theta) + sum(Z.row(i)*randPar);

    ll += y(i)*ouflow(-eta) + (1-y(i))*ouflow(eta);
  }
  ans = ll - 0.5*(dot(a, a)/sigA + dot(b, b)/sigB) - 10.0*(log(sigA) + log(sigB));
  return -ans;
}

//[[Rcpp::export]]
arma::vec grad_nlogH(arma::vec randPar, arma::vec fixPar, List data) {
  arma::vec y = data["y"];
  arma::mat X = data["x"];
  arma::mat Z = data["z"];
  int nData = Z.n_rows;
  int nRand = Z.n_cols;
  arma::vec theta = fixPar(arma::span(0, 3));
  double sigA = fixPar(4), sigB = fixPar(5);
  arma::vec a = randPar(arma::span(0, nRand/2 -1));
  arma::vec b = randPar(arma::span(nRand/2, nRand - 1));
  arma::vec eta(nData), prob(nData);
  arma::vec ans(nRand);
  eta = X*theta + Z*randPar;
  //prob = vpran(1/(1+exp(-eta)));
  arma::vec tmpa(nRand/2);
  arma::vec tmpb(nRand/2);
  tmpb = Z.cols(nRand/2 , nRand -1).t()*(y - exp(vouflow(-eta)));
  tmpa = Z.cols(0 , nRand/2 -1).t()*(y - exp(vouflow(-eta)));
  ans(arma::span(0, nRand/2 -1)) = tmpa - a/sigA;
  ans(arma::span(nRand/2, nRand -1)) = tmpb - b/sigB;
  return -ans;
}

//[[Rcpp::export]]
arma::mat hess_nlogH(arma::vec randPar, arma::vec fixPar, List data) {
  arma::mat X = data["x"];
  arma::mat Z = data["z"];
  int nData = Z.n_rows;
  int nRand = Z.n_cols;
  arma::vec theta = fixPar(arma::span(0, 3));
  double sigA = fixPar(4), sigB = fixPar(5);
  
  arma::vec eta(nData), prob(nData);
  arma::mat ans(nRand, nRand, arma::fill::zeros);
  arma::vec dda(nRand/2);
  arma::vec ddb(nRand/2);
  dda.fill(1.0/sigA);
  ddb.fill(1.0/sigB); 
  eta = X*theta + Z*randPar;
  //prob = vpran(1/(1+exp(-eta)));
  prob = exp(vouflow(-eta))%exp(vouflow(eta));
  ans(arma::span(0, nRand/2 -1), arma::span(0, nRand/2 -1)) = arma::diagmat(dda);
  ans(arma::span(nRand/2, nRand-1), arma::span(nRand/2, nRand -1)) = arma::diagmat(ddb);
  ans += Z.t()*(diagmat(prob)*Z);
  return ans;
}


//[[Rcpp::export]]
NumericVector  ldetHess_nlogH(arma::vec randPar, arma::vec fixPar, List data) {
  int nRand = randPar.n_elem;
  arma::mat fullMat(nRand, nRand);
  fullMat = hess_nlogH(randPar, fixPar, data);
  NumericVector ans(nRand);
  double val = 0.0, sign = 1.0;
  for(int i = 0; i < nRand-1;i++){
    log_det(val, sign, fullMat(arma::span(i, nRand-1), arma::span(i, nRand-1)));
    ans[i] = val;
  }
  ans[nRand-1] = log(fullMat(nRand-1, nRand-1));
  return ans;
}
