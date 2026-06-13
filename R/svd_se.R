#' Standard errors
#'
#' Robust standard errors obtained by svd of the Fisher information matrix
#'
#' @param I Fisher information matrix
#' @param J Jacobian matrix
#' @param tol tolerance
#' @return standard errors
#' @author David J. Hessen
#' @export
svd_se <- function(I, J, tol = NULL) {

  s <- svd(I)

  if (is.null(tol)) {
    tol <- max(dim(I)) *
      max(s$d) *
      .Machine$double.eps
  }

  d.inv <- ifelse(s$d > tol,
                  1 / s$d,
                  0)

  vcov <- s$v %*%
    diag(d.inv, length(d.inv)) %*%
    t(s$u)

  se <- sqrt(diag(J%*%vcov%*%t(J)))

  list(
    vcov = vcov,
    se = se,
    singular.values = s$d,
    tol = tol,
    retained = s$d > tol
  )
}
