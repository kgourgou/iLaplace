#include <RcppArmadillo.h>
using namespace Rcpp;

//[[Rcpp::export]]
double nlDen_mvt(arma::vec x, int m, double df){
  double kr = 1 + dot(x,x)/df;
  return -Rf_lgammafn(0.5*(df+m)) + Rf_lgammafn(0.5*df) + 0.5*m*log(df*M_PI) + 0.5*(df+m)*log(kr);
}
//[[Rcpp::export]]
double nlDen_mvskt(arma::vec x, double a, double c, int m, double df){
  double kr = 1 + dot(x,x)/df;
  double ans = -Rf_lgammafn(0.5*(df + m)) + Rf_lgammafn(0.5*(df +1)) + Rf_lbeta(a, c) + 0.5*log(a + c);
  ans += (a+c -1.0)*log(2.0) + 0.5*(m - 1)*log(df*M_PI) - 0.5*(df + 1)*log(1+pow(x(0),2.0)/df);
  ans += -(a + 0.5)*log(1 + x(0)/pow(a + c + pow(x(0), 2.0), 0.5)) - (c + 0.5)*log(1 - x(0)/pow(a + c+pow(x(0), 2.0), 0.5));
  ans += 0.5*(df+m)*log(kr);
  return ans;
}
//[[Rcpp::export]]
arma::vec grad_mvt(arma::vec x, int m, double df){
  double kr = 1 + dot(x,x)/df;
  return (m+df)*x/(df*kr);
}

//[[Rcpp::export]]
arma::vec grad_mvskt(arma::vec x, double a, double c, int m, double df){
  double gr1 = -(1+df)*x(0)/(df + pow(x(0), 2.0));
  gr1 += -0.5*(1+2*a)*(-x(0) + sqrt(a+c+pow(x(0),2.0)))/(a+c+pow(x(0),2.0));
  gr1 += 0.5*(1+2*c)*(x(0) + sqrt(a+c+pow(x(0),2.0)))/(a+c+pow(x(0),2.0));
  double kr = 1 + dot(x,x)/df;
  arma::vec ans = (m+df)*x/(df*kr);
  ans(0) += gr1; 
  return ans;
}

//[[Rcpp::export]]
arma::mat hess_mvt(arma::vec x, int m, double df){
  int i, j;
  double kr = 1+ dot(x,x)/df;
  arma::mat ll(m, m, arma::fill::zeros);
   for(i=0; i<(m-1); i++){
     for(j = i+1; j<m; j++){
       ll(i,j) = -(df+m)*2*x(i)*x(j)/pow(df*kr, 2.0);
   }
 }
 arma::vec dd = (df+m)*(kr - 2*pow(x,2)/df)/(df*pow(kr, 2.0));
 ll.diag() = dd;
 ll = ll + trans(ll) - diagmat(dd);
 return ll;
}

//[[Rcpp::export]]
arma::mat hess_mvskt(arma::vec x, double a, double c, int m, double df){
  int i, j;
  double kr = 1+ dot(x,x)/df;
  arma::mat ll(m, m, arma::fill::zeros);
   for(i=0; i<(m-1); i++){
     for(j = i+1; j<m; j++){
       ll(i,j) = -(df+m)*2*x(i)*x(j)/pow(df*kr, 2.0);
   }
 }
 arma::vec dd = (df+m)*(kr - 2*pow(x,2)/df)/(df*pow(kr, 2.0));
 ll.diag() = dd;
 ll = ll + trans(ll) - diagmat(dd);
 ll(0,0) += (1/(pow(a+c+x(0)*x(0), 2.0)*pow(df + x(0)*x(0), 2.0)))*(-(a+c)*(a+c-df)*df + 
           (-df*df -3*c*df*df + a*a*(1 + 3*df) + c*c*(1 + 3*df) + a*(-3*df*df + c*(2 + 6*df)))*x(0)*x(0) + 
           (a*a + 3*c + c*c + a*(3 + 2*c) - df*(3 + df))*pow(x(0),4.0) - (a + c - df)*pow(x(0),6.0) + 
           (a-c)*df*df*x(0)*sqrt(a + c + x(0)*x(0)) + 2*(a-c)*df*pow(x(0),3.0)*sqrt(a+c+ x(0)*x(0)) +
           (a-c)*pow(x(0),5.0)*sqrt(a + c + x(0)*x(0)));
 return ll;
}