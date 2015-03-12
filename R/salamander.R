nlogH_salam <- function(randPar, fixPar, data) {
  if(length(randPar)<40) stop("Random effects must be 40-dimensional vector!")
  if(length(fixPar)<6) stop("Fixed parameters must be 6-dimensional vector")
  if(length(data) != 3 | is.list(data) != TRUE) stop("data must be a list with elements y, x and z")
    .Call('iLaplace_nlogH', PACKAGE = 'iLaplace', randPar, fixPar, data)
}

grad_salam <- function(randPar, fixPar, data) {
  if(length(randPar)<40) stop("Random effects must be 40-dimensional vector!")
  if(length(fixPar)<6) stop("Fixed parameters must be 6-dimensional vector")
  if(length(data) != 3 | is.list(data) != TRUE) stop("data must be a list with elements y, x and z")
    .Call('iLaplace_grad_nlogH', PACKAGE = 'iLaplace', randPar, fixPar, data)
}

hess_salam <- function(randPar, fixPar, data) {
  if(length(randPar)<40) stop("Random effects must be 40-dimensional vector!")
  if(length(fixPar)<6) stop("Fixed parameters must be 6-dimensional vector")
  if(length(data) != 3 | is.list(data) != TRUE) stop("data must be a list with elements y, x and z")
    .Call('iLaplace_hess_nlogH', PACKAGE = 'iLaplace', randPar, fixPar, data)
}

ldetHess_nlogH <- function(randPar, fixPar, data) {
  if(length(randPar)<40) stop("Random effects must be 40-dimensional vector!")
  if(length(fixPar)<6) stop("Fixed parameters must be 6-dimensional vector")
  if(length(data) != 3 | is.list(data) != TRUE) stop("data must be a list with elements y, x and z")
    .Call('iLaplace_ldetHess_nlogH', PACKAGE = 'iLaplace', randPar, fixPar, data)
}

iLap_sal <- function(fixParm, data, sp.points, delta){
  fixPar = c(fixParm[1:4], exp(fixParm[5:6]))
  fullOpt = nlminb(rep(0, 40), nlogH_salam, gradient=grad_salam, fixPar = fixPar, data=data)
  fullOpt$hessian = hess_salam(fullOpt$par, fixPar = fixPar, data=data)
  fullOpt$ldblocks = ldetHess_nlogH(fullOpt$par, fixPar = fixPar, data=data)
  m = ncol(data$z)
  se = rep(0,m)
  for(i in 1:m){
    se[i] = sqrt(diag(solve(fullOpt$hessian[i:m,i:m])))[1]
  }
  oo = seqMat(par=fullOpt$par, se, lengthOut=sp.points, q=m, delta=delta)
  lo = oo$lo
  up = oo$up
  #up[5] = 0.033
  par.val <- oo$parVal
  
  #   lo = fullOpt$par - 20*sqrt(diag(solve(fullOpt$hessian)))
  #   up = fullOpt$par + 20*sqrt(diag(solve(fullOpt$hessian)))
  lnc = lden = 0.0
  
  # marginal density for the first prameter
  marg = function(par1) {
    tmp = function(x) nlogH_salam(c(par1, x), fixPar=fixPar, data=data)
    gr.tmp = function(x) grad_salam(c(par1, x), fixPar = fixPar, data=data)[-1]
    optTmp = nlminb(fullOpt$par[-1], tmp, gradient=gr.tmp)
    optTmp$hessian = hess_salam(c(par1,optTmp$par), fixPar = fixPar, data=data)[-1,-1]
    out  = -0.5*log(2*pi) + 0.5*(fullOpt$ldblock[1] - determinant(optTmp$hessian)$mod) -   optTmp$obj + fullOpt$obj
    return(exp(out))
  }
  #   nc.marg = integrate(Vectorize(function(x) marg(x), "x"), lower=lo[1], upper=up[1])$value
  #cat("\rMarginal normlizing constant!", nc.marg)
  fun.val <- sapply(par.val[,1], marg)
  marg.sp <- splinefun(x=par.val[,1], y=fun.val)
  #   plot(marg.sp, lo[1], up[1])
  nc.marg = integrate(marg.sp, lower=lo[1], upper=up[1])$value
  lnc = log(nc.marg) + lnc
  lden = log(marg(fullOpt$par[1])) + lden
  
  # conditionals form the second to the p-1 th component
  middleConds <- function(y, index){
    # index = 2, ..., m -1.
    no = c(1:index)
    tmp = function(xx) nlogH_salam(c(fullOpt$par[1:(index-1)], y, xx), fixPar = fixPar, data=data)
    gr.tmp = function(xx) grad_salam(c(fullOpt$par[1:(index-1)], y, xx), fixPar = fixPar, data = data)[-no]
    tmpOpt = nlminb(c(fullOpt$par[-no]), obj=tmp, gradient=gr.tmp)
    tmpOpt$hessian = hess_salam(c(fullOpt$par[1:(index-1)], y, tmpOpt$par), fixPar = fixPar, data = data)[-no,-no]
    if(index==(m-1)) {
      ldetc = log(tmpOpt$hessian)
    } else {
      ldetc = as.double(determinant(tmpOpt$hessian)$mod)
    }
    ans = -0.5*log(2*pi) + 0.5*(fullOpt$ldblock[1] - ldetc)- tmpOpt$obj + fullOpt$obj
    return(exp(ans))
  }
  
  cfun.val <- rep(0, sp.points)
  
  registerDoMC(detectCores()-1)
  # normalizes the conditionals and evaluates their density at the mode
  cfun.val = cbind(cfun.val, foreach(i=2:(m-1), .combine = cbind) %dopar% {
    #     log(integrate(Vectorize(function(x) middleConds(x, index = i), "x"), lower=lo[i], upper=up[i])$value)
    cfun.val <- sapply(par.val[,i], middleConds, index=i)
  }
  )
  for(i in 2:(m-1)){
    cond.sp <- splinefun(x=par.val[,i], y=cfun.val[,i])
    #     plot(cond.sp, lo[i], up[i])
    nc.midcond = integrate(cond.sp, lower=lo[i], upper=up[i])$value
    #nc.midcond = integrate(Vectorize(function(y) middleConds(y, index=i), "y"), lower=lo[i], upper  =up[i])$value
    #nc.midcond = integrate(Vectorize(function(x) middleConds(x, index = i), "x"), lower=lo[i], upper=up[i])$value
    lnc = lnc + log(nc.midcond)
    lden = log(middleConds(fullOpt$par[i], index=i)) + lden
    #cat(paste("\r conditional", i, "of", m-1, "\r"))
  }
  
  #last conditional posterior
  lastCond = function(par40) {
    tmp = function(x) nlogH_salam(c(fullOpt$par[c(1:(m-1))], x), fixPar = fixPar, data=data)
    out  = -0.5*log(2*pi) + 0.5*fullOpt$ldblock[m] - tmp(par40) + fullOpt$obj
    return(exp(out))
  }
  
  lcfun.val <- sapply(par.val[,m], lastCond)
  lcond.sp <- splinefun(x=par.val[,m], y=lcfun.val)
  #   plot(lcond.sp, lo[m], up[m])
  nc.lastcond = integrate(lcond.sp, lower=lo[m], upper=up[m])$value
  #cat("Experiment:", 1, "\r", sep=" ")
  lnc = log(nc.lastcond) + lnc
  lden = log(lastCond(fullOpt$par[m])) + lden
  
  return(as.double(-fullOpt$obj - lden + lnc))           
}
