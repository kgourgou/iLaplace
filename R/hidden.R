## setup diagonal tuning matrix for vector parameters

# return name of the calling function
"calling.function" <-
  function(parentheses=TRUE) {
    calling.function <- strsplit(toString(sys.call(which=-3)),",")[[1]][1]
    if (parentheses){
      calling.function <- paste(calling.function, "()", sep="")
    }
    return(calling.function)
  }

"vector.tune" <- function(mcmc.tune, K){
  if (max(is.na(mcmc.tune))){
    cat("Error: Vector tuning parameter cannot contain NAs.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)    
  }
  if (length(mcmc.tune) == 1){
    mcmc.tune <- rep(mcmc.tune, K)
  }
  if (length(mcmc.tune) != K){
    cat("Error: length(vector tuning parameter) != length(theta) or 1.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (sum(mcmc.tune <= 0) != 0) {
    cat("Error: Vector tuning parameter cannot contain negative values.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)
  }
  if (length(mcmc.tune)==1){
    return(matrix(mcmc.tune, 1, 1))
  }
  else{
    return(diag(as.double(mcmc.tune)))
  }
}