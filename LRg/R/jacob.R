#' Jacobian
#'
#' @param alpha vector of intercepts
#' @param gamma vector of regressor effects
#' @param Lambda matrix of regressor effects
#' @param Phi matrix of weights
#' @return Jacobian matrix
#' @export
jacob <- function(alpha, gamma, Lambda, Phi) {

  m <- length(alpha)
  k <- length(gamma)
  q <- ncol(Lambda)

  stopifnot(nrow(Lambda) == k)
  stopifnot(nrow(Phi) == q)
  stopifnot(ncol(Phi) == m)

  # =====================================================
  # aantallen vrije parameters
  # =====================================================

  n_alpha  <- m
  n_gamma  <- k
  n_lambda <- k * (q - 1)
  n_phi2   <- q * (m - q)

  p_theta <- n_alpha + n_gamma + n_lambda + n_phi2
  p_eta   <- m * (k + 1)

  J <- matrix(0, p_eta, p_theta)

  # =====================================================
  # kolomindices
  # =====================================================

  col_alpha <- seq_len(n_alpha)

  col_gamma <- max(col_alpha) + seq_len(n_gamma)

  if (n_lambda > 0) {
    col_lambda <- max(col_gamma) + seq_len(n_lambda)
  } else {
    col_lambda <- integer(0)
  }

  if (n_phi2 > 0) {
    start_phi <- if (n_lambda > 0) {
      max(col_lambda)
    } else {
      max(col_gamma)
    }

    col_phi <- start_phi + seq_len(n_phi2)
  } else {
    col_phi <- integer(0)
  }

  # =====================================================
  # Gamma = diag(gamma) Lambda
  # =====================================================

  Gamma <- diag(gamma) %*% Lambda

  # =====================================================
  # loop over categorieën
  # =====================================================

  for (s in seq_len(m)) {

    # -------------------------------------
    # alpha_s
    # -------------------------------------

    row_alpha <- (s - 1) * (k + 1) + 1

    J[row_alpha, col_alpha[s]] <- 1

    # -------------------------------------
    # beta_s
    # -------------------------------------

    rows_beta <- row_alpha + seq_len(k)

    phi_s <- Phi[, s, drop = FALSE]

    # =====================================
    # d beta_s / d gamma'
    # =====================================

    J_gamma <- diag(as.vector(Lambda %*% phi_s))

    J[rows_beta, col_gamma] <- J_gamma

    # =====================================
    # d beta_s / d vec(Lambda_*)'
    # =====================================

    if (q > 1) {

      phi_star <- Phi[-1, s, drop = FALSE]

      J_lambda <- kronecker(
        t(phi_star),
        diag(gamma)
      )

      J[rows_beta, col_lambda] <- J_lambda
    }

    # =====================================
    # d beta_s / d vec(Phi_2)'
    # =====================================

    if (s > q) {

      e <- numeric(m - q)
      e[s - q] <- 1

      J_phi <- kronecker(
        matrix(e, nrow = 1),
        Gamma
      )

      J[rows_beta, col_phi] <- J_phi
    }
  }

  J
}
