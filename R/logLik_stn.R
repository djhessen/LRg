#' The log-likelihood function
#'
#' @param x parameter vector
#' @param X desgin matrix
#' @param Z response dummies
#' @param q rank
#' @return log-likelihood value
#' @export
logLik_stn<-function(x, X, Z, q){
  m <- dim(Z)[2]
  k <- dim(X)[2]-1
  Gamma <- matrix(x[(m+1):(m+k*q)],nrow=k)
  U <- upper.tri(diag(q),diag=TRUE)*1
  if (q==m) {
    Phi <- U
  } else {
    Phi <- cbind(U,matrix(x[(m+k*q+1):(m+q*(k+m)-q^2)],nrow=q))
  }
  theta <- rbind(x[1:m],Gamma%*%Phi)
  eta <- X%*%theta
  sum(eta*Z)-sum(log(1+rowSums(exp(eta))))
}
