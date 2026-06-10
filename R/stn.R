#' Stereotype nominal regression
#'
#' This function can be used to fit a q-dimensional stereotype nominal regression model.
#'
#' @param formula model expression
#' @param q rank
#' @param data data frame
#' @param maxit maximum number of iterations
#' @param stv vector of starting values
#' @param csv common starting value
#' @param tol convergence tolerance
#' @return A list of
#' \item{alpha, Beta, Gamma, Phi}{Parameter estimates}
#' \item{loglik}{Value of the log-likelihood function}
#' \item{npars}{Number of estimated parameters}
#' \item{AIC}{Akaike information criterion}
#' \item{BIC}{Bayesian information criterion}
#' \item{fitted.values}{Model implied probabilities}
#' \item{residuals}{Dummy responses minus fitted values}
#' \item{pred.cat}{Model implied categories}
#' \item{obs.cat}{Observed categories}
#' @author David J. Hessen
#' @export
stn <- function(formula, q=NULL, data, stv=NULL, csv=0.001, tol=1.0e-8, maxit=100) {
  if (!is.null(q) && q==0) {
    formula <- stats::update(formula,.~1)
  }
  new_data <- auto_factor(formula, data)
  mf <- stats::model.frame(formula, new_data)
  n <- dim(mf)[1]
  X <- stats::model.matrix(formula, mf)
  k <- dim(X)[2]-1
  if (is.ordered(stats::model.response(mf))) {
    y <- as.integer(stats::model.response(mf))-1
  } else if (min(stats::model.response(mf))==1) {
    y <- as.integer(stats::model.response(mf))-1
  } else {
    y <- as.integer(stats::model.response(mf))
  }
  m <- length(unique(y))-1

  if (k>0) {
    preds <- colnames(X)[-1]
  } else if (k==0) {
    q <- 0
    preds <- NULL
  }

  if (is.null(q)) {
    q <- min(k,m)
  }

  Z <- matrix(0,nrow=n,ncol=m)
  for (j in 1:m) {
    Z[,j] <- ifelse(y==j,1,0)
  }
  if (q<0) {
    stop('The minimum rank is 0.')
  }
  if (q>min(k,m)) {
    stop('The maximum rank is ',min(k,m),'.')
  }
  freq <- table(y)
  sva <- log(freq[2:(m+1)]/freq[1])
  ll0 <- t(freq[2:(m+1)])%*%sva-n*log(1+sum(exp(sva)))
  if (q>0) {
    if (is.null(stv)) {
      sv <- c(sva,rep(0,k),rep(csv,q*(m+k)-k-q^2))
    }
    x <- sv

    it <- 0
    for (ii in 1:maxit) {
      it <- it+1
      z <- x
      Gamma <- matrix(x[(m+1):(m+k*q)],nrow=k)
      U <- upper.tri(diag(q),diag=TRUE)*1
      if (q==m) {
        Phi <- U
      } else {
        Phi <- cbind(U,matrix(x[(m+k*q+1):(m+q*(k+m)-q^2)],nrow=q))
      }

      theta <- rbind(x[1:m],Gamma%*%Phi)

      eta <- X%*%theta
      eta_max <- apply(eta,1,max)
      eta_centered <- eta-eta_max
      exp_eta <- exp(eta_centered)
      dn <- exp(-eta_max)+rowSums(exp_eta)
      P <- exp_eta/dn

      u <- rep(0,(m+q*(k+m)-q^2))
      I <- matrix(0,(m+q*(k+m)-q^2),(m+q*(k+m)-q^2))

      for (i in 1:n) {
        D2 <- t(Phi)%*%kronecker(diag(1,q),t(X[i,-1]))
        D3 <- kronecker(rbind(matrix(0,q,m-q),diag(1,m-q)),t(X[i,-1])%*%Gamma)
        D <- cbind(diag(1,m),D2,D3)

        u <- u+t(D)%*%(Z[i,]-P[i,])
        W <- diag(P[i,])-P[i,]%*%t(P[i,])
        I <- I+t(D)%*%W%*%D
      }
      dx <- solve_safe(I,u)
      # clip step
      #if (max(abs(dx)) > 1) {
      #  dx <- dx/max(abs(dx))
      #}
      x <- btls_stn(x, X, Z, q, dx, u)

      d <- sum(abs(z-x))
      if (d < tol) {
        break
      }
      cat(paste('Iteration:',it),'\r')
    }
    if (it==maxit && d>tol) {
      warning('Iteration limit reached without convergence')
    }
    J <- jacob_stn(Phi, Gamma)
    if (!inherits(tryCatch(solve(I,t(J)), error=function(e) e), "error")) {
      #se <- sqrt(diag(solve(I)))
      se <- sqrt(diag(J%*%solve(I,t(J))))
    } else {
      #se <- sqrt(diag(solve(I+1.0e-8*diag(1,ncol(I)))))
      se <- sqrt(diag(J%*%solve(I+1.0e-8*diag(1,ncol(I)),t(J))))
    }
    loglik <- logLik_stn(x, X, Z, q)
    x2 <- c(rbind(x[1:m],Gamma%*%Phi))
  } else if (q==0) {
    x <- x2 <- sva
    P <- exp(sva)/(1+sum(exp(sva)))
    I <- n*(diag(P)-P%*%t(P))
    se <- sqrt(diag(solve(I)))
    loglik <- ll0
  }
  z <- x2/se
  npars <- m+q*(k+m)-q^2
  aic <- 2*npars-2*loglik
  bic <- npars*log(n)-2*loglik
  #C <- round(cbind(x,se,z,1-stats::pnorm(abs(z))),3)
  beta <- round(cbind(x2,se,z,1-stats::pnorm(abs(z))),3)
  rownames(beta) <- paste(colnames(X),rep(1:m,rep(k+1,m)),sep=':')
  colnames(beta) <- c('Estimate','Std. Error','z value','Pr(>|z|)')
  if (q==0) {
    P <- matrix(rep(P,n),nrow=n,byrow=T)
  }
  fitv <- cbind(1-rowSums(P),P)
  res <- cbind(1-rowSums(Z),Z)-fitv
  pcat <- apply(fitv,1,which.max)
  if (q>0) {
    Beta <- Gamma%*%Phi
    CS <- 1-exp(2*(ll0-loglik)/n)
    NK <- CS/(1-exp(2*ll0/n))
    MF <- 1-(loglik/ll0)
    out <- list(alpha=x[1:m],Beta=Beta,Gamma=Gamma,Phi=Phi,loglik=loglik,
                npars=npars,AIC=aic,BIC=bic,fitted.values=fitv,residuals=res,pred.cat=pcat,
                obs.cat=y,CS=CS,NK=NK,MF=MF)
  } else if (q==0) {
    out <- list(alpha=x[1:m],loglik=loglik,npars=npars,AIC=aic,fitted.values=fitv,
                residuals=res,pred.cat=pcat,obs.cat=y)
  }
  if (!is.null(preds)) {
    cat(paste('Number of iterations:',it),'\n')
    cat(paste('Difference:',round(d,8)),'\n\n')
    cat(paste('Rank',q,'stereotype unordered'),'\n\n')
  } else {
    cat(paste('Null model'),'\n\n')
  }
  cat(paste('Number of observations:',n,'\n\n'))
  cat('Coefficients:\n\n')
  print(beta)
  invisible(out)
}



