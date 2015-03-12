nlpost_gomp <- function(param, data) {
    .Call('iLaplace_nlpost_gomp', PACKAGE = 'iLaplace', param, data)
}

grad_gomp <- function(param, data) {
    .Call('iLaplace_grad_gomp', PACKAGE = 'iLaplace', param, data)
}

hess_gomp <- function(param, data) {
    .Call('iLaplace_hess_gomp', PACKAGE = 'iLaplace', param, data)
}
