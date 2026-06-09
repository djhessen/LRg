#' The log-likelihood function
#'
#' @param x parameter vector
#' @param X desgin matrix
#' @param Z response dummies
#' @return log-likelihood value
#' @export
logLik<-function(x, X, Z){
#  sum(eta*Z)-sum(log(1+rowSums(exp(eta))))
  t(c(t(X)%*%Z))%*%x-sum(log(1+rowSums(exp(X%*%matrix(x,nrow=dim(X)[2])))))
}
