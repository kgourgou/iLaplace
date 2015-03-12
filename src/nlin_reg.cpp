#include <RcppArmadillo.h>

using namespace Rcpp;


//[[Rcpp::export]]
double dmvt(arma::vec x, arma::vec mu, arma::mat S, int p, double df, bool lg)
{
  double lans = 0.0, ldet = 0.0, sign = 1.0;
  lans += Rf_lgammafn(0.5*(df+p)) - Rf_lgammafn(0.5*df) - 0.5*p*(log(df) + 
  log(arma::datum::pi));
  arma::log_det(ldet, sign, S);
  lans += -(0.5*arma::as_scalar(ldet) + 0.5*(df+p)*log(1.0 + arma::as_scalar(trans(x-mu)
  *inv( S )*(x-mu))/df));
  if(lg == true){
    return lans;
    } else {
      return exp(lans);
  }
}


// [[Rcpp::export]]
double 
dhalfCauchy(double x, double scale, bool lg = false) 
{
  if (scale <= 0) {
    throw std::range_error("The scale parameter must be positive in dhalfCauchy.\n");
    }
  double dens = log(2 * scale) - log(M_PI*(x*x + scale*scale));
    if (lg == false) 
        dens = exp(dens);
    return dens;
}


//[[Rcpp::export]]
arma::vec rmvnorm(arma::vec mu, arma::mat S, int p)
{
  /*
  Multivariate normal random variates Y = mu + AX
  X is iid standard normal
  mu mean vector and A is such that AA^T = S
   */
  arma::mat A;
  chol(A, S); // note: this give A^TA = S, with A upper-triangular
  RNGScope scope;
  bool status = A.is_empty();
  if (status==true) Rprintf("\nCholesky decomposition in rmvnorm failed!");
  arma::vec X = arma::randn(p);
  return mu + A.t()*X;
}

//[[Rcpp::export]]
arma::vec rmvt(arma::vec mu, arma::mat S, int p, double df)
{
  RNGScope scope;
  double chi = Rf_rchisq(df);
  arma::vec mm(p, arma::fill::zeros);
  arma::vec z = rmvnorm(mm, S, p);
  return mu + z*sqrt(df/chi);
}

// [[Rcpp::export]]
double prJeff(double nu, bool lg = false)
{
  // Prior from Fonseca et al. (2008), Eq. 5
  double ans=0.0;
  ans = 0.5*(log(nu) - log(nu+3.0) + log(Rf_trigamma(0.5*nu) -
        Rf_trigamma(0.5*(nu+1)) - 2*(nu+3.0)/(nu*(nu+1)*(nu+1))));
  if(lg == false){
    ans = exp(ans);
  }
  return ans;
}

double fbX(arma::vec beta, double x1i, double x2i)
{
  double eta = 100*beta(0)/(10*beta(1) + x1i) + beta(2)*x2i 
        	+ beta(3)*x2i*x2i + beta(4)*pow(x2i,3.0) 
        	+ (beta(5) + beta(6)*x2i*x2i)*x2i*exp(-x1i/(beta(7) 
        	+ beta(8)*x2i*x2i));
    return eta;
}

arma::vec 
fbXprime(arma::vec beta, double x1i, double x2i)
{
  arma::vec ans(9, arma::fill::zeros);
	ans(0) = 100.0/(x1i + 10*beta(1));
  ans(1) = -10*beta(0)*ans(0)/(x1i + 10*beta(1));
  ans(2) = x2i;
  ans(3) = ans(2)*x2i;
  ans(4) = ans(3)*x2i;
  ans(5) = exp(-x1i/(beta(7) + x2i*x2i*beta(8)))*x2i;
  ans(6) = ans(5)*x2i*x2i;
  ans(7) = exp(-x1i/(beta(7) + x2i*x2i*beta(8)))*x1i*x2i*(beta(5) + x2i*x2i*
            beta(6))/pow(beta(7) + x2i*x2i*beta(8), 2.0);
  ans(8) = ans(7)*x2i*x2i;
  return ans;
}


arma::mat 
fbXsec(arma::vec beta, double x1i, double x2i)
{
  arma::mat ans(9,9, arma::fill::zeros);
//  ans(0,0) = 0.0;
//  ans(0,arma::span(2,8)) = ans(arma::span(2,8),0) = 0.0; 
  ans(0,1) = ans(1,0) = -1000/pow(x1i + 10*beta(1),2.0);
  ans(1,1) = 20000*beta(0)/pow(x1i + 10*beta(1),3.0);
//  ans(1,arma::span(2,8)) = ans(arma::span(2,8),1) = arma::zeros<arma::vec>(7); 
//  ans(2,arma::span(0,8)) = ans(arma::span(0,8),2) = arma::zeros<arma::vec>(7); 
//  ans(3,arma::span(0,8)) = ans(arma::span(0,8),3) = arma::zeros<arma::vec>(7);
//  ans(4,arma::span(0,8)) = ans(arma::span(0,8),4) = arma::zeros<arma::vec>(7);
//  ans(5,arma::span(0,6)) = ans(arma::span(0,6)),5) = arma::zeros<arma::vec>(7);
  ans(5,7) = ans(7,5) = exp(-x1i/(beta(7) + x2i*x2i*beta(8)))*x1i*x2i/pow(beta(7) +
            x2i*x2i*beta(8), 2.0);
  ans(5,8) = ans(8,5) = ans(5,7)*x2i*x2i;
//  ans(6,arma::span(0,6)) = ans(arma::span(0,6),6) = arma::zeros<arma::vec>(7);
  ans(6,7) = ans(7,6) = ans(8,5) = ans(5,7)*x2i*x2i;
  ans(6,8) = ans(8,6) = ans(8,5)*pow(x2i,2.0);
  ans(7,7) = exp(-x1i/(beta(7) + x2i*x2i*beta(8)))*x1i*x1i*x2i*(beta(5) + x2i*x2i*
            beta(6))/pow(beta(7) +
            x2i*x2i*beta(8), 4.0) - 2*exp(-x1i/(beta(7) + x2i*x2i*beta(8)))*x1i*x2i*(beta(5) + x2i*x2i*
            beta(6))/pow(beta(7) + x2i*x2i*beta(8), 3.0);
  ans(7,8) = ans(8,7) = exp(-x1i/(beta(7) + x2i*x2i*beta(8)))*x1i*x1i*pow(x2i,3.0)*(beta(5) + x2i*x2i*
            beta(6))/pow(beta(7) +x2i*x2i*beta(8), 4.0) -
            2*exp(-x1i/(beta(7) + x2i*x2i*beta(8)))*x1i*pow(x2i,3.0)*(beta(5) + x2i*x2i*
            beta(6))/pow(beta(7) + x2i*x2i*beta(8), 3.0);
  return ans;
}

//[[Rcpp::export]]
arma::vec grldmvt(arma::vec x, arma::vec mu, arma::mat S, int p, double df)
{
  double cc = 0.0;
  arma::vec gr(p, arma::fill::zeros);
  arma::mat Smin1 = inv(S);
  cc = (df + p)/(df + arma::as_scalar(trans(x-mu)*Smin1*(x-mu)));
  gr = Smin1*(x-mu);
  return -cc*gr;
}

//[[Rcpp::export]]
arma::mat hessldmvt(arma::vec x, arma::vec mu, arma::mat S, int p, double df)
{
  double cc = 0.0, tbsb = 0.0;
  arma::mat hh(p, p, arma::fill::zeros);
  arma::mat Smin1 = inv(S);
  tbsb = arma::as_scalar(trans(x-mu)*Smin1*(x-mu));
  cc = (df + p);
  hh = cc*(df + tbsb)*Smin1 - 2*cc*Smin1*(x-mu)*trans(x-mu)*Smin1;
  return -hh/(pow(df + tbsb, 2.0));
}

// [[Rcpp::export]]
double nlpost_lub(arma::vec beta,
   			 double lsig,
 				 arma::vec y,
 				 arma::vec x1,
 				 arma::vec x2,
 				 int n,
 				 arma::vec muBeta,
 				 arma::mat SigBeta,
 				 double sigScale)
{
  double ll = 0.0, eta = 0.0, df=2;
  int p = muBeta.n_elem + 1;
  for(int i=0; i<n; i++){
  	eta = fbX(beta, x1(i), x2(i));
//      zi = (y(i) - eta)/exp(lsig);
    ll += Rf_dnorm4(y(i), eta , exp(lsig), true);
  }
  ll += dmvt(beta, muBeta, SigBeta, p, df, true);
  ll += dhalfCauchy(exp(lsig), sigScale, true) + lsig;
  return -ll;
}

// [[Rcpp::export]]
double nlpostT_lub(arma::vec beta,
   			 double lsig,
 				 double lnu,
 				 arma::vec y,
 				 arma::vec x1,
 				 arma::vec x2,
 				 int n,
 				 arma::vec muBeta,
 				 arma::mat SigBeta,
 				 double sigScale)
{
  double ll = 0.0, eta = 0.0, zi = 0.0, df=2;
  int p = muBeta.n_elem + 2;
  for(int i=0; i<n; i++){
  	eta = fbX(beta, x1(i), x2(i));
    zi = (y(i) - eta)/exp(lsig);
    ll += Rf_dt(zi, exp(lnu), true) - lsig;
  }
  ll += dmvt(beta, muBeta, SigBeta, p, df, true) +
        dhalfCauchy(exp(lsig), sigScale, true) +
        lsig + prJeff(exp(lnu), true) + lnu;
  return -ll;
}

//[[Rcpp::export]]
arma::vec 
grad_lub(arma::vec beta,
        double lsig,
 			  arma::vec y,
 			  arma::vec x1,
 			  arma::vec x2,
 			  int n,
 			  arma::vec muBeta,
 			  arma::mat SigBeta,
 			  double sigScale) 
{
  arma::vec ans(10, arma::fill::zeros),
  etaprime(9, arma::fill::zeros);
  double eta = 0.0, Df = 0.0, zi = 0.0, df = 2.0;
  for(int i=0; i<n; i++){
    eta = fbX(beta, x1(i), x2(i));
    zi = (y(i)-eta)/exp(lsig);
    etaprime = fbXprime(beta, x1(i), x2(i));
    Df = exp(-lsig)*zi;
    ans(arma::span(0,8)) += Df*etaprime;
    ans(9) += -1 + zi*zi;
  }
  arma::vec grb = grldmvt(beta, muBeta, SigBeta, 9, df);
  ans(arma::span(0,8)) += grb;
  ans(9) += 1 - 2*exp(2*lsig)/(exp(2*lsig) + sigScale*sigScale);
  return -ans;
}

//[[Rcpp::export]]
arma::vec 
gradT_lub(arma::vec beta,
     	double lsig,
      double lnu,
 			arma::vec y,
 			arma::vec x1,
 			arma::vec x2,
 			int n,
 			arma::vec muBeta,
 			arma::mat SigBeta,
 			double sigScale) 
{
  arma::vec ans(11, arma::fill::zeros),
  etaprime(9, arma::fill::zeros);
  double eta = 0.0, Df = 0.0, zi = 0.0, df = 2.0;
  for(int i=0; i<n; i++){
    eta = fbX(beta, x1(i), x2(i));
    zi = (y(i)-eta)/exp(lsig);
    etaprime = fbXprime(beta, x1(i), x2(i));
    Df = exp(-lsig)*(1+exp(lnu))*zi/(exp(lnu) + zi*zi);
    
    ans(arma::span(0,8)) += Df*etaprime;
    ans(9) += -1.0 + (1+exp(lnu))*zi*zi/(exp(lnu) + zi*zi);
    ans(10) += 0.5*exp(lnu) - 0.5*exp(lnu)*(1+exp(lnu))/(exp(lnu) + zi*zi)
           + 0.5*exp(lnu)*lnu - 0.5*exp(lnu)*log(exp(lnu)+zi*zi)
           - 0.5*exp(lnu)*Rf_digamma(0.5*exp(lnu)) + 0.5*exp(lnu)
           *Rf_digamma(0.5*(1.0+exp(lnu)));
  }
  arma::vec grb = grldmvt(beta, muBeta, SigBeta, 9, df);
  ans(arma::span(0,8)) += grb;
  ans(9) += 1 - 2*exp(2*lsig)/(exp(2*lsig) + sigScale*sigScale);
  ans(10) += 1 + 0.5*
          (1 - exp(lnu)/(3+exp(lnu)) + 
          (2*exp(-lnu)*(3+9*exp(lnu)+2*exp(2*lnu))/pow(1+exp(lnu), 3.0) 
          + 0.5*exp(lnu)*Rf_tetragamma(0.5*exp(lnu))
          - 0.5*exp(lnu)*Rf_tetragamma(0.5*(1+exp(lnu))))
          /(-2.0*(3.0+exp(lnu))*exp(-lnu)/(pow(1+exp(lnu), 2.0)) +
          Rf_trigamma(0.5*exp(lnu)) -
          Rf_trigamma(0.5*(1+exp(lnu))))
          );
  return -ans;
}

//[[Rcpp::export]]
arma::mat 
hess_lub(arma::vec beta,
        double lsig,
   		  arma::vec y,
 			  arma::vec x1,
 			  arma::vec x2,
 			  int n,
 			  arma::vec muBeta,
 			  arma::mat SigBeta,
 			  double sigScale) 
{
  double zi=0.0, df = 2;
  arma::mat hh(10,10, arma::fill::zeros),
  etas(9,9,arma::fill::zeros);
  arma::vec etap(9, arma::fill::zeros);
  for(int i=0; i<n; i++){
    zi = (y(i) - fbX(beta, x1(i), x2(i)))/exp(lsig);
    etap = fbXprime(beta, x1(i), x2(i));
    etas = fbXsec(beta, x1(i), x2(i));
    hh(arma::span(0,8),arma::span(0,8)) += -exp(-2*lsig)*
                                          (etap*trans(etap) - 
                                          exp(lsig)*zi*etas);
    hh(arma::span(0,8),9) +=  -2*exp(-lsig)*zi*etap;
    hh(9,9) += -2*zi*zi;
    }
  hh(9,arma::span(0,8)) = trans(hh(arma::span(0,8),9));
  hh(arma::span(0,8),arma::span(0,8)) += hessldmvt(beta, muBeta, SigBeta, 9, df);
  hh(9,9) += -4*exp(2*lsig)*sigScale*sigScale/pow(exp(2*lsig) + sigScale*sigScale, 2.0);
  return -hh;
}

//[[Rcpp::export]]
arma::mat 
hessT_lub(arma::vec beta,
        double lsig,
        double lnu,
     	  arma::vec y,
 			  arma::vec x1,
 			  arma::vec x2,
 			  int n,
 			  arma::vec muBeta,
 			  arma::mat SigBeta,
 			  double sigScale) 
{
  double zi=0.0, df = 2;
  arma::mat hh(11,11, arma::fill::zeros),
  etas(9,9,arma::fill::zeros);
  arma::vec etap(9, arma::fill::zeros);
  for(int i=0; i<n; i++){
    zi = (y(i) - fbX(beta, x1(i), x2(i)))/exp(lsig);
    etap = fbXprime(beta, x1(i), x2(i));
    etas = fbXsec(beta, x1(i), x2(i));
    
//    for(int g=0; g<9; g++){
//      for(int h = 0; h<9; h++){
//        hh(g,h) += (1+exp(lnu))*(-exp(-2*lsig)*etap(g)*etap(h)/(exp(lnu)+zi*zi)+
//                        2*exp(-2*lsig)*zi*zi*etap(g)*etap(h)/pow(exp(lnu)+zi*zi, 2.0)+
//                        exp(-lsig)*zi*etas(g,h)/(exp(lnu)+zi*zi));
//      }
//    }
    hh(arma::span(0,8),arma::span(0,8)) += -(1+exp(lnu))*(exp(-2*lsig)*etap*trans(etap)/(exp(lnu)+zi*zi)-
                        2*exp(-2*lsig)*zi*zi*etap*trans(etap)/pow(exp(lnu)+zi*zi, 2.0)-exp(-lsig)*zi*etas/
                        (exp(lnu)+zi*zi));
                        
    hh(arma::span(0,8),9) +=  -2*exp(-lsig)*(1+exp(lnu))*zi*etap/(exp(lnu)+zi*zi)+
                            2*exp(-lsig)*(1+exp(lnu))*pow(zi,3.0)*etap/pow(exp(lnu)+
                            zi*zi,2.0);
                            
    hh(9,9) += -2*(1+exp(lnu))*zi*zi/(exp(lnu)+zi*zi) + 2*(1+exp(lnu))*pow(zi,4.0)/
              pow(exp(lnu)+zi*zi,2);
              
    hh(arma::span(0,8),10) +=  -exp(-lsig+lnu)*(1+exp(lnu))*zi*etap/pow(exp(lnu)+zi*zi,2.0)+
                            exp(-lsig+lnu)*zi*etap/(exp(lnu)+zi*zi);
                            
    hh(9,10) += -exp(lnu)*(1+exp(lnu))*zi*zi/pow(exp(lnu)+zi*zi,2.0)+
                            exp(lnu)*zi*zi/(exp(lnu)+zi*zi);
                            
    hh(10,10) += exp(lnu) + 0.5*exp(2*lnu)*(1+exp(lnu))/pow(exp(lnu)+zi*zi,2.0) -
               exp(2*lnu)/(exp(lnu)+zi*zi) - 0.5*exp(lnu)*(1+exp(lnu))/(exp(lnu)+zi*zi) +
               0.5*exp(lnu)*lnu - 0.5*exp(lnu)*log(exp(lnu)+zi*zi) - 0.5*exp(lnu)
               *Rf_digamma(0.5*exp(lnu)) + 0.5*exp(lnu)*Rf_digamma(0.5*(exp(lnu)+1)) 
               -0.25*exp(2*lnu)*Rf_trigamma(0.5*exp(lnu)) + 0.25*exp(2*lnu)*
               Rf_trigamma(0.5*(exp(lnu)+1));
    }
  hh(9,arma::span(0,8)) = trans(hh(arma::span(0,8),9));
  hh(10,arma::span(0,8)) = trans(hh(arma::span(0,8),10));
  hh(10,9) = hh(9,10);
  hh(arma::span(0,8),arma::span(0,8)) += hessldmvt(beta, muBeta, SigBeta, 9, df);
  hh(9,9) += -4*exp(2*lsig)*sigScale*sigScale/pow(exp(2*lsig) + sigScale*sigScale, 2.0);
  hh(10,10) += 0.5*exp(2*lnu)/pow(3+exp(lnu),2.0) - 0.5*exp(lnu)/(3+exp(lnu)) -pow(
               -2/pow(1+exp(lnu),2.0) + 4*(3+exp(lnu))/pow(1+exp(lnu),3.0) + 2*exp(-lnu)
               *(3+exp(lnu))/pow(1+exp(lnu),2.0) + 0.5*exp(lnu)*Rf_tetragamma(0.5*exp(lnu))
               - 0.5*exp(lnu)*Rf_tetragamma(0.5*(1+exp(lnu))),2.0)/(2*pow(-2*exp(-lnu)
               *(3+exp(lnu))/pow(1+exp(lnu),2.0) + Rf_trigamma(0.5*exp(lnu)) 
               - Rf_trigamma(0.5*(1+exp(lnu))), 2.0)) + (8*exp(lnu)/pow(1+exp(lnu),3.0) 
               + 2/pow(1+exp(lnu), 2.0) -12*exp(lnu)*(3+exp(lnu))/pow(1+exp(lnu), 4.0) 
               - 4*(3+exp(lnu))/pow(1+exp(lnu),3.0) - 2*exp(-lnu)*(3+exp(lnu))/pow(1+exp(lnu),2.0)
               + 0.5*exp(lnu)*Rf_tetragamma(0.5*exp(lnu)) - 0.5*exp(lnu)*Rf_tetragamma(0.5*
               (1+exp(lnu))) + 0.25*exp(2*lnu)*Rf_pentagamma(0.5
               *exp(lnu)) - 0.25*exp(2*lnu)*Rf_pentagamma(0.5*(1+exp(lnu))))/(2*(-2*exp(-lnu)
               *(3+exp(lnu))/pow(1+exp(lnu),2.0)+Rf_trigamma(0.5*exp(lnu)) - 
               Rf_trigamma(0.5*(1+exp(lnu)))));
  return -hh;
}


double user_fun_eval(SEXP fun, arma::vec theta, SEXP myframe) {
  SEXP R_fcall;
  if(!Rf_isFunction(fun)) throw std::range_error("`fun' must be a function");
  if(!Rf_isEnvironment(myframe)) throw std::range_error("myframe must be an environment");
  PROTECT(R_fcall = Rf_lang2(fun, R_NilValue));
  SETCADR(R_fcall, wrap(theta));
  SEXP funval;
  PROTECT(funval = Rf_eval(R_fcall, myframe));

//  if (!Rf_isReal(funval)) throw std::range_error("`fun' must return a double");
  if (!Rf_isReal(funval))  Rf_warning("WW`fun' must return a double");
  double fv = REAL(funval)[0];
//  if (fv == R_PosInf) throw std::range_error("`fun' returned +Inf");
  if (fv == R_PosInf) Rf_warning("WW`fun' returned +Inf");
//  if (R_IsNaN(fv) || R_IsNA(fv)) throw std::range_error("`fun' returned NaN or NA");
  if (R_IsNaN(fv) || R_IsNA(fv)) Rf_warning("WW`fun' returned NaN or NA");
  UNPROTECT(2);
  return fv;
}

// [[Rcpp::export]]
arma::mat
MCMCmetrop_cpp(Function lfun,
          arma::vec start,
          arma::mat V,
          int mcmc,
          int burnin,
          double df,
          int verbose,
          Environment myframe)
{
  int p = start.n_elem, nacc=1;
  int pV = V.n_cols;
  double logAccRat = 0.0, rUnif = 0.0;
  if (p!=pV) {
    throw std::range_error("start and V must have the same # rows.\n");
  }
//  if (nthin!=mcmc) {
//    throw std::range_error("thinvec must have mcmc elements.\n");
//  }
  arma::vec propPar(p, arma::fill::zeros);
  arma::mat ans(p, mcmc+burnin, arma::fill::zeros);
  ans.col(0) = start;
  for(int i = 1; i < (mcmc+burnin); i++){
    // draw from the jumping density
    propPar = rmvt(ans.col(i-1), V, p, df);
    // acceptance ratio of prior and jumping density
    logAccRat = user_fun_eval(lfun, propPar, myframe)-
                user_fun_eval(lfun, ans.col(i-1), myframe);
    
    if(R_finite(logAccRat) != true){
      ans.col(i) = ans.col(i-1);
    } else {          
      rUnif = as_scalar(arma::randu(1));
    
      ans.col(i) = ans.col(i-1);
      if (log(rUnif) <= logAccRat) {
        nacc += 1;
        ans.col(i) = propPar;
        }
    // print some information on the simulation process
      if(i%verbose == 0 ) Rprintf("\rMCMC: iteration %d of %d, accpted %d, acc. ratio %f\r",
      i, mcmc+burnin, nacc, static_cast<double>(nacc)/static_cast<double>(i+1));
    
    // for stopping R computations
      }
      R_CheckUserInterrupt();
    }
  Rprintf("\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
  Rprintf("MCMC accepteance ratio was: %3.5f", 
    static_cast<double>(nacc) / static_cast<double>(mcmc+burnin));
  Rprintf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
  
//  arma::mat ans2 = ans(arma::span::all,arma::span(burnin-1,mcmc-1));
  
  return ans;
}

// [[Rcpp::export]]
double 
mlChib_cpp(Function lfun,
          arma::vec hatPar,
          arma::mat V,
          arma::mat mcmc,
          double df,
          int verbose,
          Environment myframe)
{
  int nNum = mcmc.n_cols, p = mcmc.n_rows;
  double den = 0.0, num = 0.0;
  arma::vec minNum(2, arma::fill::zeros), minDen(2, arma::fill::zeros);
  minNum(0) = minDen(0)= 1.0;
  
  arma::vec propPar(p, arma::fill::zeros);
  double lfun_hatpar = user_fun_eval(lfun, hatPar, myframe);
  
  for(int i=0; i<nNum; i++){
    minNum(1) = exp(lfun_hatpar - 
                user_fun_eval(lfun, mcmc.col(i), myframe));
    num += min(minNum)*dmvt(mcmc.col(i),hatPar, V, p, df, false);
    
    propPar = rmvt(hatPar, V, p, df);
    minDen(1) = exp(user_fun_eval(lfun, propPar, myframe)-lfun_hatpar);
    den += min(minDen);
        // print some information on the simulation process
    if(i%verbose == 0 ) Rprintf("\rChib's marg. lik at iteration %d of %d, is %f",
    i, nNum, log(num)-log(den));
    
    // for stopping R computations
    R_CheckUserInterrupt();
  }
  return log(num) - log(den);
}

//[[Rcpp::export]]
double
mlIS_cpp(Function lfun,
         arma::vec hatPar,
         arma::mat V,
         int nsim,
         int p,
         double df,
         int verbose,
         Environment myframe)
{
  double malik = 0.0, wsamp = 0.0;
  arma::vec simpar(p, arma::fill::zeros);
  
  for (int i=0; i<nsim; i++){
    simpar = rmvt(hatPar, V, p, df);
    wsamp = user_fun_eval(lfun, simpar, myframe) - dmvt(simpar, hatPar, V, p, df, true);
    if(R_finite(wsamp) != true){
      wsamp = -1000.0;
    }
    malik += exp(wsamp);
    if(i%verbose == 0 ) Rprintf("\rImportance sampling marg. lik at iteration %d of %d, is %f",
    i, nsim, malik/static_cast<double>(i));
  }
  return log(malik) - log(nsim);
}

//[[Rcpp::export]]
double 
nlpost_bod2(arma::vec beta,
    		  double lsig,
				  arma::vec y,
				  arma::vec x,
          int n,
				  arma::vec muBeta,
				  arma::mat SigBeta,
				  double sigScale) {
  double ll = 0.0, eta = 0.0, df = 2.0;
  
  for(int i=0; i<n; i++){
    eta = beta(0)*(1-exp(-x(i)/beta(1)));
    ll += Rf_dnorm4(y(i), eta, exp(lsig), true);
  }
  ll += dmvt(beta, muBeta, SigBeta, 2, df, true);
  ll += dhalfCauchy(exp(lsig), sigScale, true) + lsig;
  return -ll;
}

// [[Rcpp::export]]
double nlpostT_bod2(arma::vec beta,
   			 double lsig,
 				 double lnu,
 				 arma::vec y,
 				 arma::vec x,
         int n,
 				 arma::vec muBeta,
 				 arma::mat SigBeta,
 				 double sigScale)
{
  	double ll = 0.0, eta = 0.0, zi = 0.0, df = 2.0;
    // likelihood
  	for(int i=0; i<n; i++){
  		eta = beta(0)*(1-exp(-x(i)/beta(1)));
      zi = (y(i) - eta)/exp(lsig);
    ll += Rf_dt(zi, exp(lnu), true) - lsig;
  }
  //prior for beta
  ll += dmvt(beta, muBeta, SigBeta, 2, df, true) 
  		+ dhalfCauchy(exp(lsig), sigScale, true)+ lsig//prior for sigma
      + prJeff(exp(lnu), true) + lnu;// prior for nu
  return -ll;
}

// [[Rcpp::export]]
arma::vec 
grad_bod2(arma::vec beta,
    		  double lsig,
				  arma::vec y,
				  arma::vec x,
          int n,
				  arma::vec muBeta,
				  arma::mat SigBeta,
				  double sigScale)
{
  double df = 2, zi = 0.0;
  arma::vec gr(3, arma::fill::zeros);

  //update with data
  for(int i=0; i<n; i++){
    zi = (y(i) - beta(0)*(1-exp(-x(i)/beta(1))))/exp(lsig);
    gr(0) += (1.0 - exp(-x(i)/beta(1)))*zi/exp(lsig);
    gr(1) += -beta(0)*exp(-x(i)/beta(1))*x(i)*zi/
            (exp(lsig)*beta(1)*beta(1));
    gr(2) += -1.0+ zi*zi;
  }
  //update with prior
  arma::vec grb = grldmvt(beta, muBeta, SigBeta, 2, df);
  gr(0) += grb(0);
  gr(1) += grb(1);
  gr(2) += 1 - 2*exp(2*lsig)/(exp(2*lsig) + sigScale*sigScale);
  return -gr;
}

// [[Rcpp::export]]
arma::mat 
hess_bod2(arma::vec beta,
      	  double lsig,
				  arma::vec y,
				  arma::vec x,
          int n,
				  arma::vec muBeta,
				  arma::mat SigBeta,
				  double sigScale)
{
  double zi=0.0, df = 2;
  arma::mat hh(3,3, arma::fill::zeros);
  for(int i=0; i<n; i++){
    zi = (y(i) - beta(0)*(1-exp(-x(i)/beta(1))))/exp(lsig); 
    hh(0,1) += - beta(0)*exp(-x(i)/beta(1))*x(i)*(-1.0 + exp(-x(i)/beta(1)))/
              (pow(beta(1), 2.0)*exp(lsig)*exp(lsig)) - exp(-x(i)/beta(1))*x(i)
              *zi/(pow(beta(1), 2.0)*exp(lsig));
    hh(0,2) += 2*exp(-2*lsig)*(-1.0 + exp(-x(i)/beta(1)))*zi*exp(lsig);
    hh(0,0) += -pow(-1.0 + exp(-x(i)/beta(1)), 2.0)/(exp(lsig)*exp(lsig));
    hh(1,2) += 2*beta(0)*exp(-2*lsig -x(i)/beta(1))*x(i)*zi*exp(lsig)/pow(beta(1), 2.0);
    hh(1,1) += - beta(0)*beta(0)*exp(-2*x(i)/beta(1))*x(i)*x(i)/
              (pow(beta(1), 4.0)*exp(lsig)*exp(lsig)) +
              2*beta(0)*exp(-x(i)/beta(1))*x(i)*zi/(pow(beta(1), 3.0)*exp(lsig))
              - beta(0)*exp(-x(i)/beta(1))*x(i)*x(i)*zi/
              (pow(beta(1), 4.0)*exp(lsig)); 
//    hh(2,2) += 1.0/pow(exp(lsig),2.0) - 3.0*zi*zi/pow(exp(lsig),2.0);
    hh(2,2) += -2*exp(-2*lsig)*zi*zi*exp(lsig)*exp(lsig);
  }
  hh = hh + trans(hh) - diagmat(hh);
  
  hh(arma::span(0,1),arma::span(0,1)) += hessldmvt(beta, muBeta, SigBeta, 2, df = 2); 
  hh(2,2) += -4*exp(2*lsig)*sigScale*sigScale/pow(exp(2*lsig) + sigScale*sigScale, 2.0);
  return -hh;
}

// [[Rcpp::export]]
arma::vec 
gradT_bod2(arma::vec beta,
     		 double lsig,
 				 double lnu,
 				 arma::vec y,
 				 arma::vec x,
         int n,
 				 arma::vec muBeta,
 				 arma::mat SigBeta,
 				 double sigScale)
{
  double df = 2, zi = 0.0, zi2 = 0.0;
  arma::vec gr(4, arma::fill::zeros);
  
    //update with data
  for(int i=0; i<n; i++){
    zi = (y(i) - beta(0)*(1-exp(-x(i)/beta(1))))/exp(lsig);
    zi2 = zi*zi;
    gr(0) += -(-1.0 + exp(-x(i)/beta(1)))*zi*(1+exp(lnu))
            /(exp(lsig)*(zi*zi + exp(lnu)));
            
    gr(1) += -beta(0)*exp(-x(i)/beta(1))*x(i)*zi*(1+exp(lnu))/
            (beta(1)*beta(1)*exp(lsig)*(zi2 + exp(lnu)));
            
    gr(2) += -1.0 + zi2*(1.0+exp(lnu))/(zi2 + exp(lnu));
    
    gr(3) += 0.5*exp(lnu) - 0.5*exp(lnu)*(1+exp(lnu))/(exp(lnu) + zi2)
           + 0.5*exp(lnu)*lnu - 0.5*exp(lnu)*log(exp(lnu)+zi2)
           - 0.5*exp(lnu)*Rf_digamma(0.5*exp(lnu)) + 0.5*exp(lnu)
           *Rf_digamma(0.5*(1.0+exp(lnu)));
  }
  //update with prior
  arma::vec grb = grldmvt(beta, muBeta, SigBeta, 2, df);
  gr(0) += grb(0);
  gr(1) += grb(1);
  
  gr(2) +=  1 - 2*exp(2*lsig)/(exp(2*lsig) + sigScale*sigScale);

  gr(3) += 1 + 0.5*
          (1 - exp(lnu)/(3+exp(lnu)) + 
          (2*exp(-lnu)*(3+9*exp(lnu)+2*exp(2*lnu))/pow(1+exp(lnu), 3.0) 
          + 0.5*exp(lnu)*Rf_tetragamma(0.5*exp(lnu))
          - 0.5*exp(lnu)*Rf_tetragamma(0.5*(1+exp(lnu))))
          /(-2.0*(3.0+exp(lnu))*exp(-lnu)/(pow(1+exp(lnu), 2.0)) +
          Rf_trigamma(0.5*exp(lnu)) -
          Rf_trigamma(0.5*(1+exp(lnu))))
          );
  return -gr;
}

// [[Rcpp::export]]
arma::mat 
hessT_bod2(arma::vec beta,
       	 double lsig,
 				 double lnu,
 				 arma::vec y,
 				 arma::vec x,
         int n,
 				 arma::vec muBeta,
 				 arma::mat SigBeta,
 				 double sigScale)
{
  double zi=0.0, df = 2, mEb2 = 0.0;
  arma::mat hh(4,4, arma::fill::zeros);
  for(int i=0; i<n; i++){
    mEb2 = 1-exp(-x(i)/beta(1));
    zi = (y(i) - beta(0)*mEb2)/exp(lsig); 
    
    hh(0,0) += 2*exp(-2*lsig)*mEb2*mEb2*(1+exp(lnu))*zi*zi/
              pow(exp(lnu) + zi*zi, 2.0) - exp(-2*lsig)*mEb2*mEb2*(1+exp(lnu))
              /(exp(lnu)+zi*zi);
              
    hh(1,1) += 2*pow(beta(0)*x(i)*zi, 2.0)*exp(-2.0*x(i)/beta(1) -2*lsig)*(1+exp(lnu))
              /pow(beta(1)*beta(1)*(exp(lnu)+zi*zi), 2.0) 
              - exp(-2*x(i)/beta(1) - 2*lsig)*(1+exp(lnu))*pow(beta(0)*x(i),2.0)/(pow(beta(1), 4.0)*(exp(lnu)+zi*zi))
              + 2*beta(0)*exp(-x(i)/beta(1) -lsig)*(1+exp(lnu))*x(i)*zi/(pow(beta(1), 3.0)*(exp(lnu)+zi*zi)) 
              - beta(0)*exp(-x(i)/beta(1)-lsig)*(1+exp(lnu))*x(i)*x(i)*zi/(pow(beta(1), 4.0)*(exp(lnu)+zi*zi));
              
    hh(2,2) += 2*(1+exp(lnu))*pow(zi,4.0)/pow(exp(lnu)+zi*zi,2.0) - 2*(1+exp(lnu))*zi*zi
              /(exp(lnu)+zi*zi);
              
    hh(3,3) += exp(lnu) + 0.5*exp(2*lnu)*(1+exp(lnu))/pow(exp(lnu)+zi*zi,2.0) -
               exp(2*lnu)/(exp(lnu)+zi*zi) - 0.5*exp(lnu)*(1+exp(lnu))/(exp(lnu)+zi*zi) +
               0.5*exp(lnu)*lnu - 0.5*exp(lnu)*log(exp(lnu)+zi*zi) - 0.5*exp(lnu)
               *Rf_digamma(0.5*exp(lnu)) + 0.5*exp(lnu)*Rf_digamma(0.5*(exp(lnu)+1)) 
               -0.25*exp(2*lnu)*Rf_trigamma(0.5*exp(lnu)) + 0.25*exp(2*lnu)*Rf_trigamma(0.5*(exp(lnu)+1));
               
    hh(0,1) += 2*beta(0)*exp(-x(i)/beta(1) - 2.0*lsig)*(-mEb2)*(1+exp(lnu))*x(i)*zi*zi/
              pow(beta(1)*(exp(lnu) + zi*zi), 2.0) - beta(0)*exp(-x(i)/beta(1) -2.0*lsig)*(-mEb2)*(1+exp(lnu))*x(i)/
              (beta(1)*beta(1)*(exp(lnu)+zi*zi)) - exp(-x(i)/beta(1)-lsig)*(1.0+exp(lnu))*x(i)*zi
              /(beta(1)*beta(1)*(exp(lnu)+zi*zi));
    hh(0,2) += -2*exp(-lsig)*(-mEb2)*(1+exp(lnu))*pow(zi,3.0)/pow(exp(lnu)+zi*zi, 2.0) 
              + 2*exp(-lsig)*(-mEb2)*(1+exp(lnu))*zi/(exp(lnu)+zi*zi);
    hh(0,3) += exp(-lsig+lnu)*(-mEb2)*(1+exp(lnu))*zi/pow(exp(lnu)+zi*zi, 2.0) - 
               exp(-lsig+lnu)*(-mEb2)*zi/(exp(lnu)+zi*zi);
    hh(1,2) += -2*beta(0)*exp(-x(i)/beta(1) - lsig)*(1+exp(lnu))*x(i)*pow(zi,3.0)
              /pow(beta(1)*(exp(lnu)+zi*zi), 2.0) + 2*beta(0)*exp(-x(i)/beta(1) -lsig)*(1+exp(lnu))*x(i)*zi/
              (pow(beta(1),2.0)*(exp(lnu)+zi*zi));
    hh(1,3) += beta(0)*exp(-x(i)/beta(1) -lsig + lnu)*(1+exp(lnu))*x(i)*zi/pow(beta(1)*(exp(lnu)+zi*zi), 2.0)
              - beta(0)*exp(-x(i)/beta(1) -lsig + lnu)*x(i)*zi/(beta(1)*beta(1)*(exp(lnu)+zi*zi));
    hh(2,3) += -exp(lnu)*(1+exp(lnu))*zi*zi/pow(exp(lnu)+zi*zi, 2.0) + exp(lnu)*zi*zi/(exp(lnu)+zi*zi);
  }
    hh = hh + trans(hh) - diagmat(hh);
    hh(arma::span(0,1),arma::span(0,1)) += hessldmvt(beta, muBeta, SigBeta, 2, df = 2);
    
    hh(2,2) += -4*exp(2*lsig)*sigScale*sigScale/pow(exp(2*lsig) + sigScale*sigScale, 2.0);
    
    hh(3,3) += 0.5*exp(2*lnu)/pow(3+exp(lnu),2.0) - 0.5*exp(lnu)/(3+exp(lnu)) 
              -pow(
               -2/pow(1+exp(lnu),2.0) + 4*(3+exp(lnu))/pow(1+exp(lnu),3.0) + 2*exp(-lnu)
               *(3+exp(lnu))/pow(1+exp(lnu),2.0) + 0.5*exp(lnu)*Rf_tetragamma(0.5*exp(lnu))
               - 0.5*exp(lnu)*Rf_tetragamma(0.5*(1+exp(lnu))),2.0)/(2*pow(-2*exp(-lnu)
               *(3+exp(lnu))/pow(1+exp(lnu),2.0) + Rf_trigamma(0.5*exp(lnu)) 
               - Rf_trigamma(0.5*(1+exp(lnu))), 2.0)) + (8*exp(lnu)/pow(1+exp(lnu),3.0) 
               + 2/pow(1+exp(lnu), 2.0) -12*exp(lnu)*(3+exp(lnu))/pow(1+exp(lnu), 4.0) 
               - 4*(3+exp(lnu))/pow(1+exp(lnu),3.0) - 2*exp(-lnu)*(3+exp(lnu))/pow(1+exp(lnu),2.0)
               + 0.5*exp(lnu)*Rf_tetragamma(0.5*exp(lnu)) - 0.5*exp(lnu)*Rf_tetragamma(0.5*
               (1+exp(lnu))) + 0.25*exp(2*lnu)*Rf_pentagamma(0.5
               *exp(lnu)) - 0.25*exp(2*lnu)*Rf_pentagamma(0.5*(1+exp(lnu))))/(2*(-2*exp(-lnu)
               *(3+exp(lnu))/pow(1+exp(lnu),2.0)+Rf_trigamma(0.5*exp(lnu)) - 
               Rf_trigamma(0.5*(1+exp(lnu)))));
  return -hh;
}


arma::vec seq_C(double a, double cent, double b, int lengthOut) {
  arma::vec out(lengthOut);
  out.fill(a);
  out[lengthOut-1] = b;
  double width = (b-a)/(lengthOut-1);
  for(int i=1; i<(lengthOut -1); i++){
    out(i) = out(i-1) + width;
  }
  out(lengthOut/2) = cent;
  return out;
}

// [[Rcpp::export]]
List seqMat(arma::vec par, arma::vec se, int lengthOut, int q, double delta){
  arma::vec av = par - delta*se;
  arma::vec bv = par + delta*se;
  arma::mat out(lengthOut, q, arma::fill::zeros);
  for(int i=0; i<q; i++){
    out.col(i) = seq_C(av(i), par(i), bv(i), lengthOut);
  }
  return List::create(Named("parVal") = out,
                      Named("lo") = av,
                      Named("up") = bv);
}
