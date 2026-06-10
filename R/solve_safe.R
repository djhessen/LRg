#' Solving Idx=u
#'
#' @param I information matrix or Hessian
#' @param u score or gradient
#' @param lambda ridge penalty parameter
#' @return ridge information matrix or Hessian
#' @export
solve_safe <- function(I, u, lambda = 1.0e-8) {

  p <- ncol(I)

  if (rcond(I) < 1e-12) {
    I <- I + lambda * diag(p)
  }

  solve(I, u)
}
