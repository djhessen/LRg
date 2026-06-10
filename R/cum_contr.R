#' Cumulative binary recoding
#'
#' @param x ordered factor
#' @return cumulative contrast
#' @export
cum_contr <- function(x) {

  if (!is.factor(x)) {
    x <- factor(x, ordered = TRUE)
  }

  K <- nlevels(x)

  # Contrast matrix
  C <- matrix(0, nrow = K, ncol = K - 1)

  for (i in 2:K) {
    C[i, 1:(i - 1)] <- 1
  }

  rownames(C) <- levels(x)
  colnames(C) <- paste0("c", 1:(K - 1))

  return(C)
}
