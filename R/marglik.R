ChibML <- function(logfun, theta.star, tune, V, mcmcsamp, df, verbose){
  
  my.env <- environment(fun = logfun)
  
  tune <- vector.tune(tune, length(theta.star))
  
  if (nrow(V) != ncol(V) || nrow(V) != length(theta.star)){
    cat("V not of appropriate dimension.\n")
    stop("Check V and theta.star and call ChibML() again. \n",
         call.=FALSE)     
  }
  CC <- NULL
  try(CC <- chol(V), silent=TRUE)
  if (is.null(CC)){
    cat("V not positive definite.\n")
    stop("Check V and call ChibML() again. \n",
         call.=FALSE)     
  }
  V <- tune %*% V %*% tune
  
  ans = .Call('iLaplace_mlChib_cpp',
              PACKAGE = 'iLaplace',
              logfun,
              theta.star,
              V,
              mcmcsamp,
              df,
              verbose,
              my.env)
  return(ans)
}

isML <- function(logfun, nsim, theta.hat, tune, V, df, verbose){
  my.env <- environment(fun = logfun)
  
  tune <- vector.tune(tune, length(theta.hat))
  if (nrow(V) != ncol(V) || nrow(V) != length(theta.hat)){
    cat("V not of appropriate dimension.\n")
    stop("Check V and theta.star and call ChibML() again. \n",
         call.=FALSE)     
  }
  CC <- NULL
  try(CC <- chol(V), silent=TRUE)
  if (is.null(CC)){
    cat("V not positive definite.\n")
    stop("Check V and call ChibML() again. \n",
         call.=FALSE)     
  }
  V <- tune %*% V %*% tune
  p <- ncol(V)
  ans = .Call('iLaplace_mlIS_cpp', 
              PACKAGE = 'iLaplace',
              logfun,
              theta.hat,
              V,
              nsim, 
              p, 
              df,
              verbose,
              my.env)
  return(ans)
}