#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
double nlpost_gomp(NumericVector param, NumericVector data) {
  int n = data.size();
  double ll = 0.0, xpar = 0.0;
  for(int i=0; i<n; i++){
    xpar = exp(param[1])*data[i];
    ll += param[0] + xpar - exp(param[0]-param[1])*(exp(xpar)-1.0);
  }
    ll += sum(Rcpp::dnorm(param, 0.0, 10.0, 1));
  return -ll;
}

//[[Rcpp::export]]
NumericVector grad_gomp(NumericVector param, NumericVector data) {
  int n = data.size();
  double ll = 0.0, xpar = 0.0;
  NumericVector ans(2);
  ans[0] = ans[1] = 0.0;
  for(int i=0; i<n; i++){
    xpar = exp(param[1])*data[i];
    ans[0] += 1.0 - exp(param[0] - param[1])*(-1.0 + exp(xpar));
    ans[1] += exp(param[0] - param[1])*(-1.0 + exp(xpar)) - exp(xpar + param[0])*data[i] + xpar;
  }
  ans[0] += -param[0]/100.0;
  ans[1] += -param[1]/100.0;
  return (-1.0)*ans;
}

//[[Rcpp::export]]
NumericMatrix hess_gomp(NumericVector param, NumericVector data) {
  int n = data.size();
  double ll = 0.0, xpar = 0.0;
  NumericMatrix ans(2, 2);
  for(int i=0; i<n; i++){
    xpar = exp(param[1])*data[i];
    ans(0,0) += -exp(param[0] - param[1])*(-1.0 + exp(xpar));
    ans(1,1) += -exp(param[0] - param[1])*(-1.0 + exp(xpar)) + exp(xpar + param[0])*data[i] +
    xpar - exp(xpar + param[0] + param[1])*pow(data[i], 2.0);
    ans(0,1) += exp(param[0] - param[1])*(-1.0 + exp(xpar)) - exp(xpar + param[0])*data[i];
  }
  ans(1,0) = ans(0,1);
  ans(0,0) += -1.0/100.0;
  ans(1,1) += -1.0/100.0;
  ans(0,0) = ans(0,0)* -1.0;
  ans(1,0) = ans(1,0) * -1.0;
  ans(0,1) = ans(0,1) * -1.0;
  ans(1,1) = ans(1,1) * -1.0;
  return ans;
}