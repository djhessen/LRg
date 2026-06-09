#' Backtracking line search
#'
#' @param x parameter vector
#' @param X design matrix
#' @param Z response dummies
#' @param dx search direction
#' @param u score vector
#' @param alpha step size
#' @param rho contraction factor
#' @param c Armijo parameter
#' @param maxit maximum number of iterations
#' @return new parameter vector
#' @export
btls_stn <- function(x, X, Z, q, dx, u, alpha = 1, rho = 0.5, c = 1e-4, maxit = 20) {

  ll_old <- logLik_stn(x, X, Z, q)

  for (i in 1:maxit) {

    x_new <- x + alpha * dx
    ll_new <- logLik_stn(x_new, X, Z, q)

    if (is.na(ll_new) || is.infinite(ll_new)) {
      alpha <- rho * alpha
      next
    }

    # Armijo condition
    if (ll_new >= ll_old + c * alpha * sum(u * dx)) {
      return(x_new)
    }

    # reduce step size
    alpha <- rho * alpha

  }

  # if nothing works, return the smallest step
  x_new <- x + alpha * dx
  x_new
}
