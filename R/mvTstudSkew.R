nlDen_mvt <- function(x, m, df) {
    .Call('iLaplace_nlDen_mvt', PACKAGE = 'iLaplace', x, m, df)
}

nlDen_mvskt <- function(x, a, c, m, df) {
    .Call('iLaplace_nlDen_mvskt', PACKAGE = 'iLaplace', x, a, c, m, df)
}

grad_mvt <- function(x, m, df) {
    .Call('iLaplace_grad_mvt', PACKAGE = 'iLaplace', x, m, df)
}

grad_mvskt <- function(x, a, c, m, df) {
    .Call('iLaplace_grad_mvskt', PACKAGE = 'iLaplace', x, a, c, m, df)
}

hess_mvt <- function(x, m, df) {
    .Call('iLaplace_hess_mvt', PACKAGE = 'iLaplace', x, m, df)
}

hess_mvskt <- function(x, a, c, m, df) {
    .Call('iLaplace_hess_mvskt', PACKAGE = 'iLaplace', x, a, c, m, df)
}

iLap_mvt = function(m, df) {
  # full optimization and hessian
  fullOpt = nlminb(rep(0,m), objective = nlDen_mvt, gradient=grad_mvt, m=m, df=df)
  fullOpt$hessian = hess_mvt(fullOpt$par, m, df)
  fullOpt$ldet = as.double(determinant(fullOpt$hessian)$mod)
  
  lo = fullOpt$par - 100*sqrt(diag(solve(fullOpt$hessian)))
  up = fullOpt$par + 100*sqrt(diag(solve(fullOpt$hessian)))
  lnc = lden = 0.0
  
  # approximate marginal denstiy for the first component
  marg = function(y) {
    tmp = function(xx) nlDen_mvt(c(y, xx), m=m, df=df)
    gr.tmp = function(xx) grad_mvt(c(y, xx), m=m, df=df)[-1]
    tmpOpt = nlminb(fullOpt$par[2:m], objective = tmp, gradient=gr.tmp)
    tmpOpt$hessian = hess_mvt(c(y, tmpOpt$par), m=m, df=df)[-1,-1]
    ldetc = determinant(tmpOpt$hessian)$mod
    ans = -0.5*log(2*pi) + 0.5*(fullOpt$ldet - ldetc) - tmpOpt$obj + fullOpt$obj
    return(exp(ans))
  }
  
  # normalizes the marginal density and computes the density at the mode
  nc.marg = integrate(Vectorize(function(x) marg(x), "x"), lower=lo[1], upper=up[1])$value
  lnc = log(nc.marg) + lnc
  lden = log(marg(fullOpt$par[1])) + lden
  
  # conditionals form the second to the p-1 th component
  middleConds <- function(y, index){
    # index = 2, ..., m -1.
    no = c(1:index)
    tmp = function(xx) nlDen_mvt(c(fullOpt$par[1:(index-1)], y, xx), m=m, df=df)
    gr.tmp = function(xx) grad_mvt(c(fullOpt$par[1:(index-1)], y, xx), m=m, df=df)[-no]
    tmpOpt = nlminb(c(fullOpt$par[-no]), objective = tmp, gradient=gr.tmp)
    tmpOpt$hessian = hess_mvt(c(fullOpt$par[1:(index-1)], y, tmpOpt$par), m=m, df=df)[-no,-no]
    if(index==(m-1)) {
      ldetc = log(tmpOpt$hessian)
    } else {
      ldetc = determinant(tmpOpt$hessian)$mod
    }
    ans = -0.5*log(2*pi) + 0.5*(fullOpt$ldet - ldetc)- tmpOpt$obj + fullOpt$obj
    return(exp(ans))
  }
  
  # normalizes the conditionals and evaluates their density at the mode
  for(i in 2:(m-1)){
    nc.midcond = integrate(Vectorize(function(x) middleConds(x, index = i), "x"), lower=lo[i], upper=up[i])$value
    lnc = lnc + log(nc.midcond)
    lden = log(middleConds(fullOpt$par[i], index=i)) + lden
  }
  
  # The conditional of pth given the rest
  lastCond = function(y) {
    tt = -nlDen_mvt(c(fullOpt$par[1:(m-1)], y), m=m, df=df) + fullOpt$objective
    return(exp(tt))
  }
  
  # Normalized the last conditional and evaluates the density at the mode
  nc.lastcond = integrate(Vectorize(function(x) lastCond(x), "x"), lower=lo[m], upper=up[m])$value
  
  lnc = log(nc.lastcond) + lnc
  lden = log(lastCond(fullOpt$par[m])) + lden
  
  ans = c(-fullOpt$obj - lden + lnc, 0.5*m*log(2*pi) - fullOpt$obj - 0.5*fullOpt$ldet)
  return(ans)
}

iLap_mvskt = function(a, c, m, df) {
  # full optimization and hessian
  fullOpt = nlminb(rep(0,m), objective = nlDen_mvskt, gradient=grad_mvskt, a=a, c=c, m=m, df=df)
  fullOpt$hessian = hess_mvskt(fullOpt$par, a=a, c=c, m=m, df=df)
  fullOpt$ldet = as.double(determinant(fullOpt$hessian)$mod)
  
  lo = fullOpt$par - 100*sqrt(diag(solve(fullOpt$hessian)))
  up = fullOpt$par + 100*sqrt(diag(solve(fullOpt$hessian)))
  lnc = lden = 0.0
  
  # approximate marginal denstiy for the first component
  marg = function(y) {
    tmp = function(xx) nlDen_mvskt(c(y, xx), a=a, c=c, m=m, df=df)
    gr.tmp = function(xx) grad_mvskt(c(y, xx), a=a, c=c, m=m, df=df)[-1]
    tmpOpt = nlminb(fullOpt$par[2:m], objective = tmp, gradient=gr.tmp)
    tmpOpt$hessian = hess_mvskt(c(y, tmpOpt$par), a=a, c=c, m=m, df=df)[-1,-1]
    ldetc = determinant(tmpOpt$hessian)$mod
    ans = -0.5*log(2*pi) + 0.5*(fullOpt$ldet - ldetc) - tmpOpt$obj + fullOpt$obj
    return(exp(ans))
  }
  
  # normalizes the marginal density and computes the density at the mode
  nc.marg = integrate(Vectorize(function(x) marg(x), "x"), lower=lo[1], upper=up[1])$value
  lnc = log(nc.marg) + lnc
  lden = log(marg(fullOpt$par[1])) + lden
  
  # conditionals form the second to the p-1 th component
  middleConds <- function(y, index){
    # index = 2, ..., m -1.
    no = c(1:index)
    tmp = function(xx) nlDen_mvskt(c(fullOpt$par[1:(index-1)], y, xx), a=a, c=c, m=m, df=df)
    gr.tmp = function(xx) grad_mvskt(c(fullOpt$par[1:(index-1)], y, xx), a=a, c=c, m=m, df=df)[-no]
    tmpOpt = nlminb(c(fullOpt$par[-no]), objective = tmp, gradient=gr.tmp)
    tmpOpt$hessian = hess_mvskt(c(fullOpt$par[1:(index-1)], y, tmpOpt$par), a=a, c=c, m=m, df=df)[-no,-no]
    if(index==(m-1)) {
      ldetc = log(tmpOpt$hessian)
    } else {
      ldetc = determinant(tmpOpt$hessian)$mod
    }
    ans = -0.5*log(2*pi) + 0.5*(fullOpt$ldet - ldetc)- tmpOpt$obj + fullOpt$obj
    return(exp(ans))
  }
  
  # normalizes the conditionals and evaluates their density at the mode
  for(i in 2:(m-1)){
    nc.midcond = integrate(Vectorize(function(x) middleConds(x, index = i), "x"), lower=lo[i], upper=up[i])$value
    lnc = lnc + log(nc.midcond)
    lden = log(middleConds(fullOpt$par[i], index=i)) + lden
  }
  
  # The conditional of pth given the rest
  lastCond = function(y) {
    tt = -nlDen_mvskt(c(fullOpt$par[1:(m-1)], y), a=a, c=c, m=m, df=df) + fullOpt$objective
    return(exp(tt))
  }
  
  # Normalized the last conditional and evaluates the density at the mode
  nc.lastcond = integrate(Vectorize(function(x) lastCond(x), "x"), lower=lo[m], upper=up[m])$value
  
  lnc = log(nc.lastcond) + lnc
  lden = log(lastCond(fullOpt$par[m])) + lden
  
  ans = c(-fullOpt$obj - lden + lnc, 0.5*m*log(2*pi) - fullOpt$obj - 0.5*fullOpt$ldet)
  return(ans)
}

