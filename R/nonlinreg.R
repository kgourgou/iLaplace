# This file was generated by Rcpp::compileAttributes
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

dmvt <- function(x, mu, S, p, df, lg) {
    .Call('iLaplace_dmvt', PACKAGE = 'iLaplace', x, mu, S, p, df, lg)
}

dhalfCauchy <- function(x, scale, lg = FALSE) {
    .Call('iLaplace_dhalfCauchy', PACKAGE = 'iLaplace', x, scale, lg)
}

rmvnorm <- function(mu, S, p) {
    .Call('iLaplace_rmvnorm', PACKAGE = 'iLaplace', mu, S, p)
}

rmvt <- function(mu, S, p, df) {
    .Call('iLaplace_rmvt', PACKAGE = 'iLaplace', mu, S, p, df)
}

prJeff <- function(nu, lg = FALSE) {
    .Call('iLaplace_prJeff', PACKAGE = 'iLaplace', nu, lg)
}

grldmvt <- function(x, mu, S, p, df) {
    .Call('iLaplace_grldmvt', PACKAGE = 'iLaplace', x, mu, S, p, df)
}

hessldmvt <- function(x, mu, S, p, df) {
    .Call('iLaplace_hessldmvt', PACKAGE = 'iLaplace', x, mu, S, p, df)
}

nlpost_lub <- function(beta, lsig, y, x1, x2, n, muBeta, SigBeta, sigScale) {
    .Call('iLaplace_nlpost_lub', PACKAGE = 'iLaplace', beta, lsig, y, x1, x2, n, muBeta, SigBeta, sigScale)
}

nlpostT_lub <- function(beta, lsig, lnu, y, x1, x2, n, muBeta, SigBeta, sigScale) {
    .Call('iLaplace_nlpostT_lub', PACKAGE = 'iLaplace', beta, lsig, lnu, y, x1, x2, n, muBeta, SigBeta, sigScale)
}

grad_lub <- function(beta, lsig, y, x1, x2, n, muBeta, SigBeta, sigScale) {
    .Call('iLaplace_grad_lub', PACKAGE = 'iLaplace', beta, lsig, y, x1, x2, n, muBeta, SigBeta, sigScale)
}

gradT_lub <- function(beta, lsig, lnu, y, x1, x2, n, muBeta, SigBeta, sigScale) {
    .Call('iLaplace_gradT_lub', PACKAGE = 'iLaplace', beta, lsig, lnu, y, x1, x2, n, muBeta, SigBeta, sigScale)
}

hess_lub <- function(beta, lsig, y, x1, x2, n, muBeta, SigBeta, sigScale) {
    .Call('iLaplace_hess_lub', PACKAGE = 'iLaplace', beta, lsig, y, x1, x2, n, muBeta, SigBeta, sigScale)
}

hessT_lub <- function(beta, lsig, lnu, y, x1, x2, n, muBeta, SigBeta, sigScale) {
    .Call('iLaplace_hessT_lub', PACKAGE = 'iLaplace', beta, lsig, lnu, y, x1, x2, n, muBeta, SigBeta, sigScale)
}

nlpost_bod2 <- function(beta, lsig, y, x, n, muBeta, SigBeta, sigScale) {
    .Call('iLaplace_nlpost_bod2', PACKAGE = 'iLaplace', beta, lsig, y, x, n, muBeta, SigBeta, sigScale)
}

nlpostT_bod2 <- function(beta, lsig, lnu, y, x, n, muBeta, SigBeta, sigScale) {
    .Call('iLaplace_nlpostT_bod2', PACKAGE = 'iLaplace', beta, lsig, lnu, y, x, n, muBeta, SigBeta, sigScale)
}

grad_bod2 <- function(beta, lsig, y, x, n, muBeta, SigBeta, sigScale) {
    .Call('iLaplace_grad_bod2', PACKAGE = 'iLaplace', beta, lsig, y, x, n, muBeta, SigBeta, sigScale)
}

hess_bod2 <- function(beta, lsig, y, x, n, muBeta, SigBeta, sigScale) {
    .Call('iLaplace_hess_bod2', PACKAGE = 'iLaplace', beta, lsig, y, x, n, muBeta, SigBeta, sigScale)
}

gradT_bod2 <- function(beta, lsig, lnu, y, x, n, muBeta, SigBeta, sigScale) {
    .Call('iLaplace_gradT_bod2', PACKAGE = 'iLaplace', beta, lsig, lnu, y, x, n, muBeta, SigBeta, sigScale)
}

hessT_bod2 <- function(beta, lsig, lnu, y, x, n, muBeta, SigBeta, sigScale) {
    .Call('iLaplace_hessT_bod2', PACKAGE = 'iLaplace', beta, lsig, lnu, y, x, n, muBeta, SigBeta, sigScale)
}

seqMat <- function(par, Hess, lengthOut, q, delta) {
  .Call('iLaplace_seqMat', PACKAGE = 'iLaplace', par, Hess, lengthOut, q, delta)
}

# iLap_bod2 <- function(fullOpt, ff, ff.gr, ff.hess, lb, ub){
#   
#   # marginal density for the first prameter
#   marg <- function(x, fullOpt) {
#     tmp = function(xx) ff(c(x, xx))
#     gr.tmp = function(xx) ff.gr(c(x, xx))[-1]
#     optTmp = nlminb(fullOpt$par[-1], objective = tmp, gradient=gr.tmp)
#     optTmp$hessian = ff.hess(c(x,optTmp$par))[-1,-1]
#     ldetc = determinant(optTmp$hessian)$mod
#     out  = -0.5*log(2*pi) + 0.5*(determinant(fullOpt$hessian)$mod - ldetc) - optTmp$obj + fullOpt$obj
#     return(exp(out))
#   }
#   margv <- Vectorize(function(x) marg(x, fullOpt = fullOpt), "x")
#   
#   # conditional densities
#   middleCond1 <- function(x, fullOpt){
#     tmp = function(xx) ff(c(fullOpt$par[1], x, xx))
#     gr.tmp = function(xx) ff.gr(c(fullOpt$par[1], x, xx))[-c(1:2)]
#     tmpOpt = nlminb(c(fullOpt$par[-c(1:2)]), objective = tmp, gradient=gr.tmp)
#     tmpOpt$hessian = ff.hess(c(fullOpt$par[1], x, tmpOpt$par))[-c(1:2),-c(1:2)]
#     ldetc = determinant(tmpOpt$hessian)$mod
#     ans = -0.5*log(2*pi) + 0.5*(determinant(fullOpt$hessian[-1,-1])$mod - ldetc)- tmpOpt$obj + fullOpt$obj
#     return(exp(ans))
#   }
#   middleCond1v <- Vectorize(function(x) middleCond1(x, fullOpt = fullOpt), "x")
#   
#   middleCond2 <- function(x, fullOpt){
#     tmp = function(xx) ff(c(fullOpt$par[1:2], x, xx))
#     gr.tmp = function(xx) ff.gr(c(fullOpt$par[1:2], x, xx))[-c(1:3)]
#     tmpOpt = nlminb(c(fullOpt$par[-c(1:3)]), objective = tmp, gradient=gr.tmp)
#     tmpOpt$hessian = ff.hess(c(fullOpt$par[1:2], x, tmpOpt$par))[-c(1:3),-c(1:3)]
#     ldetc = log(tmpOpt$hessian)
#     ans = -0.5*log(2*pi) + 0.5*(determinant(fullOpt$hessian[-c(1:2),-c(1:2)])$mod - ldetc)- tmpOpt$obj + fullOpt$obj
#     return(exp(ans))
#   }
#   middleCond2v <- Vectorize(function(x) middleCond2(x, fullOpt = fullOpt), "x")
#   
#   #last conditional posterior
#   lastCond <- function(x, fullOpt) {
#     tmp = function(x) ff(c(fullOpt$par[c(1:3)], x))
#     out  = -0.5*log(2*pi) + 0.5*log(fullOpt$hessian[4,4]) - tmp(x) + fullOpt$obj
#     return(exp(out))
#   }
#   lastCondv <- Vectorize(function(x) lastCond(x, fullOpt=fullOpt), "x")
#   
#   # normalizing constants of the univariate densities
#   log.ncost = log(integrate(margv, lower=lb[1], upper=ub[1])$value)+
#     log(integrate(middleCond1v, lower=lb[2], upper=ub[2])$value)+
#     log(integrate(middleCond2v, lower=lb[3], upper=ub[3])$value)+
#     log(integrate(lastCondv, lower=lb[4], upper=ub[4])$value)
#   
#   # non-normalized univariate densities at the mode
#   log.den = log(margv(fullOpt$par[1]))+
#     log(middleCond1v(fullOpt$par[2]))+
#     log(middleCond2v(fullOpt$par[3]))+
#     log(lastCondv(fullOpt$par[4]))
#   
#   return(list(iLap = -fullOpt$obj - log.den + log.ncost, Lap = 0.5*3*log(2*pi) - fullOpt$obj -0.5*fullOpt$ldblock[1]))
# }

# iLaplace and Laplace approximation for lubricant data
iLap_an <- function(fullOpt, ff, ff.gr, ff.hess, sp.points, se, delta)
{ 
  m = length(fullOpt$par)
  oo = seqMat(par=fullOpt$par, se, lengthOut=sp.points, q=m, delta=delta)
  lo = oo$lo
  up = oo$up
  #up[5] = 0.033
  par.val <- oo$parVal
  fullOpt$ldblock =  c()
  for(i in 1:(m-1)) fullOpt$ldblock[i] = as.double(determinant(fullOpt$hessian[c(i:(m)), c(i:(m))])$mod)
  #   lo = fullOpt$par - 11*sqrt(diag(solve(fullOpt$hessian)))
  #   up = fullOpt$par + 11*sqrt(diag(solve(fullOpt$hessian)))
  # #   up[5] = .033
  log.ncost = log.den = 0.0
  
  # marginal density for the first prameter
  marg <- function(x) {
    tmp = function(xx) ff(c(x, xx))
    gr.tmp = function(xx) ff.gr(c(x, xx))[-1]
    optTmp = nlminb(start = fullOpt$par[-1], objective = tmp, gradient = gr.tmp)
    optTmp$hessian = ff.hess(c(x,optTmp$par))[-1,-1]
    ldetc = determinant(optTmp$hessian)$mod
    out  = -0.5*log(2*pi) + 0.5*(fullOpt$ldblock[1] - ldetc) - optTmp$obj + fullOpt$obj
    return(exp(as.double(out)))
  }
  
#   pdf("Rplots.pdf")
  fun.val <- sapply(par.val[,1], marg)
  marg.sp <- splinefun(x=par.val[,1], y=fun.val)
#   plot(marg.sp, lo[1], up[1])
  nc.marg = integrate(marg.sp, lower=lo[1], upper=up[1])$value
  
  #   nc.marg = integrate(Vectorize(marg, "x"), lower=lb[1], upper=ub[1])$value
  cat("\rMarginal normlizing constant!", nc.marg)
  log.ncost = log(nc.marg) + log.ncost
  log.den = log(marg(fullOpt$par[1])) + log.den
  
  # conditionals form the second to the p-1 th component
  middleConds <- function(x, index){
    # index = 2, ..., m -1.
    no = c(1:index)
    tmp = function(xx) ff(c(fullOpt$par[1:(index-1)], x, xx))
    gr.tmp = function(xx) ff.gr(c(fullOpt$par[1:(index-1)], x, xx))[-no]
    tmpOpt = nlminb(start = c(fullOpt$par[-no]), objective = tmp, gradient = gr.tmp)
    tmpOpt$hessian = ff.hess(c(fullOpt$par[1:(index-1)], x, tmpOpt$par))[-no,-no]
    if(index==(m-1)) {
      ldetc = log(tmpOpt$hessian)
    } else {
      ldetc = as.double(determinant(tmpOpt$hessian)$mod)
    }
    ans = -0.5*log(2*pi) + 0.5*(fullOpt$ldblock[index] - ldetc)- tmpOpt$obj + fullOpt$obj
    return(exp(ans))
  }
  
  # normalizes the conditionals and evaluates their density at the mode
  for(i in 2:(m-1)){
    cfun.val <- sapply(par.val[,i], middleConds, index=i)
    cond.sp <- splinefun(x=par.val[,i], y=cfun.val)
#     plot(cond.sp, lo[i], up[i])
    nc.midcond = integrate(cond.sp, lower=lo[i], upper=up[i])$value
    #nc.midcond = integrate(Vectorize(function(y) middleConds(y, index=i), "y"), lower=lo[i], upper  =up[i])$value
    log.ncost = log.ncost + log(nc.midcond)
    log.den = log.den + log(middleConds(fullOpt$par[i], index=i))
    #     log.den = log(cond.sp(fullOpt$par[i])) + lden
    cat(paste("\r\n @@@@@@@ Finished with cond.", i, "of", m-1, "@@@@@@@\r"))
    cat(paste("\n @@@@@@@ Cond. norm. const", i, "was", nc.midcond, "@@@@@@@"))
  }
  
  
  #   for(i in 2:(m-1)){
  #     nc.midcond = integrate(Vectorize(function(x) middleConds(x, index = i), "x"), lower=lb[i], upper=ub[i])$value
  #     log.ncost = log.ncost + log(nc.midcond)
  #     log.den = log(middleConds(fullOpt$par[i], index=i)) + log.den
  #     cat(paste("\r conditional", i, "of", m-1, "\r"))
  #   }
  
  #last conditional posterior
  lastCond = function(x) {
    tmp = function(xx) ff(c(fullOpt$par[c(1:(m-1))], xx))
    out  = -0.5*log(2*pi) + 0.5*log(fullOpt$hessian[m,m]) - tmp(x) + fullOpt$obj
    return(exp(out))
  }
  
  lcfun.val <- sapply(par.val[,m], lastCond)
  lcond.sp <- splinefun(x=par.val[,m], y=lcfun.val)
#   plot(lcond.sp, lo[m], up[m])
  nc.lastcond = integrate(lcond.sp, lower=lo[m], upper=up[m])$value
  #   nc.lastcond = integrate(Vectorize(lastCond, "x"), lower=lb[m], upper=ub[m])$value
  #cat("Experiment:", 1, "\r", sep=" ")
  log.ncost = log(nc.lastcond) + log.ncost
  log.den = log(lastCond(fullOpt$par[m])) + log.den
#   dev.off()
  
  return(list(iLap = -fullOpt$obj - log.den + log.ncost, Lap = 0.5*m*log(2*pi) - fullOpt$obj -0.5*fullOpt$ldblock[1]))  
}