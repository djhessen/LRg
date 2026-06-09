#' Multinomial logistic regression
#'
#' @param formula model expression
#' @param data data frame
#' @param maxit maximum number of iterations
#' @return A list of
#' \item{coefficients}{Parameter estimates}
#' \item{fitted.values}{The model implied probabilities}
#' \item{residuals}{Dummy responses minus fitted values}
#' \item{pred.cat}{The model implied categories}
#' @export
mlr <- function(formula, data, maxit=100) {
  new_data <- auto_factor(formula, data)
  mf <- model.frame(formula, new_data)
  y <- as.integer(model.response(mf))
  X <- model.matrix(formula, mf)

  n <- dim(mf)[1]
  k <- dim(X)[2]-1
  m <- length(unique(y))-1

  Z <- matrix(,nrow=n,ncol=m)
  for (j in 1:m) {
    Z[,j] <- ifelse(y==j,1,0)
  }

  sv <- rep(0,m*(k+1))
  x <- sv

  it <- 0
  for (ii in 1:maxit) {
    it <- it+1
    z <- x
    theta <- matrix(x,nrow=k+1)
    eta <- X%*%theta
    eta_max <- apply(eta,1,max)
    eta_centered <- eta-eta_max
    exp_eta <- exp(eta_centered)
    dn <- exp(-eta_max)+rowSums(exp_eta)
    P <- exp_eta/dn

    u <- rep(0,m*(k+1))
    I <- matrix(0,m*(k+1),m*(k+1))

    for (i in 1:n) {
      D <- kronecker(diag(1,m),t(X[i,]))
      u <- u+t(D)%*%(Z[i,]-P[i,])
      W <- diag(P[i,])-P[i,]%*%t(P[i,])
      I <- I+t(D)%*%W%*%D
    }
    dx <- solve_safe(I,u)
    # clip step
    if (max(abs(dx)) > 1) {
      dx <- dx/max(abs(dx))
    }
    x <- btls(x, X, Z, dx, u)

    d <- sum(abs(z-x))
    if (d < 1.0e-6) {
      break
    }
    cat(paste('Iteration:',it),'\r')
  }
  if (!inherits(tryCatch(solve(I), error=function(e) e), "error")) {
    se <- sqrt(diag(solve(I)))
  } else {
    se <- sqrt(diag(solve(I+1.0e-8*diag(1,ncol(I)))))
  }
  z <- x/se
  C <- round(cbind(x,se,z,1-pnorm(abs(z))),3)
  rownames(C) <- paste(colnames(X),rep(1:m,rep(k+1,m)),sep=':')
  colnames(C) <- c('Estimate','Std. Error','z value','Pr(>|z|)')
  fitv <- cbind(1-rowSums(P),P)
  res <- cbind(1-rowSums(Z),Z)-fitv
  pcat <- apply(fitv,1,which.max)
  out <- list(coefficients=x,fitted.values=fitv,
              residuals=res,pred.cat=pcat,obs.cat=y)
  cat(paste('Number of iterations:',it),'\n')
  cat(paste('Difference:',round(d,8)),'\n\n')
  cat('Coefficients:\n\n')
  print(C)
  invisible(out)
}



