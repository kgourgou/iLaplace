MHmcmcR <- function(logfun, burnin, mcmc, thin, tune, V, df, theta.init, verbose){
  my.env <- environment(fun = logfun)
  tune <- vector.tune(tune, length(theta.init))
  
  if (nrow(V) != ncol(V) || nrow(V) != length(theta.init)){
    cat("V not of appropriate dimension.\n")
    stop("Check V and theta.init and call MCMCmetrop1R() again. \n",
         call.=FALSE)     
  }
  CC <- NULL
  try(CC <- chol(V), silent=TRUE)
  if (is.null(CC)){
    cat("V not positive definite.\n")
    stop("Check V and call MCMCmetrop1R() again. \n",
         call.=FALSE)     
  }
  V <- tune %*% V %*% tune
  
  sample <- .Call('iLaplace_MCMCmetrop_cpp',
                  PACKAGE = 'iLaplace',
                  logfun,
                  theta.init,
                  V,
                  mcmc,
                  burnin,
                  df,
                  verbose,
                  my.env) 
  sample <- mcmc(data=t(sample[,-c(1:burnin)]), start=1, end=mcmc, thin=thin)
  return(sample)
}