#' Jacobian
#'
#' @param Phi matrix
#' @param Gamma matrix
#' @return Jacobian
#' @export
jacob_stn <- function(Phi, Gamma) {
  m <- ncol(Phi)
  q <- nrow(Phi)
  k <- nrow(Gamma)

  nrow_J <- m * (1 + k)
  ncol_J <- m + k*q + q*(m - q)

  J <- matrix(0, nrow_J, ncol_J)

  idx_alpha <- 1:m
  idx_Gamma <- (m + 1):(m + k*q)
  idx_Phi2  <- (m + k*q + 1):ncol_J

  row_ptr <- 1

  for (s in 1:m) {
    J[row_ptr, idx_alpha[s]] <- 1
    row_ptr <- row_ptr + 1

    J[row_ptr:(row_ptr + k - 1), idx_Gamma] <-
      kronecker(t(Phi[, s, drop = FALSE]), diag(k))

    if (s > q) {
      e <- matrix(0, 1, m - q)
      e[1, s - q] <- 1

      J[row_ptr:(row_ptr + k - 1), idx_Phi2] <-
        kronecker(e, Gamma)
    }

    row_ptr <- row_ptr + k
  }

  return(J)
}
