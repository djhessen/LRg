#' Stereotype ordered regression
#'
#' This function can be used to fit a q-dimensional stereotype ordered regression model.
#'
#' @param formula model expression
#' @param q rank of the matrix of regressor effects. If not supplied, then it is set to its maximum.
#' @param data data frame
#' @param maxit maximum number of iterations
#' @param tol convergence tolerance
#' @return A list of the following components
#' \item{alpha}{vector of estimated intercepts}
#' \item{Beta}{matrix of estimated regressor effects}
#' \item{Gamma, Phi}{matrices of parameter estimates}
#' \item{loglik}{Value of the log-likelihood function}
#' \item{npars}{Number of estimated parameters}
#' \item{AIC}{Akaike information criterion}
#' \item{BIC}{Bayesian information criterion}
#' \item{CS, NK, MF}{Cox and Snell, Nagelkerke, and McFadden pseudo-R-squared values}
#' @author David J. Hessen
#' @export
sto <- function(formula, q=NULL, data, tol=1.0e-8, maxit=100) {
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

  if (k==1) {stop('Program is not yet ready for one regressor variable.')}

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

  if (q==1) {
    sv <- c(sva,rep(0,k),seq(1.1,(1.1+(m-2)*0.1),0.1))
  } else if (q>1 && q<m) {
    svP <- matrix(rep(seq(1.1,(1.1+(m-q-1)*0.1),0.1),q),nrow=q,byrow=TRUE)
    sv <- c(sva,rep(0.01,k),rep(0.01,k*(q-1)),c(svP))
  } else if (q==m) {
    sv <- c(sva,rep(0.01,q*(k+m)-q^2))
  }

  ll0 <- t(freq[2:(m+1)])%*%sva-n*log(1+sum(exp(sva)))
  if (q>0) {

    x <- sv
    grad <- function(x) {
      if (k>1) {
        gamma <- diag(x[(m+1):(m+k)])
      } else if (k==1) {
        gamma <- x[m+1]
      }

      if (q==1) {
        Lambda <- rep(1,k)
      } else if (q>1) {
        Lambda <- matrix(c(rep(1,k),x[(m+k+1):(m+k*q)]),nrow=k)
      }
      U <- upper.tri(diag(q),diag=TRUE)*1
      if (q==m) {
        Phi <- U
      }
      if (q<m && k>1) {
        Phi <- cbind(U,matrix(x[(m+k*q+1):(m+q*(k+m)-q^2)],nrow=q))
      } else if (q==1 && k==1) {
        Phi <- matrix(c(U,x[(m+2):(2*m)]),nrow=1)
      }

      if (k==1) {
        theta <- matrix(c(x[1:m],c(gamma)*Phi),nrow=2,byrow=TRUE)
      } else if (k>1) {
        theta <- rbind(x[1:m],gamma%*%Lambda%*%Phi)
      }
      eta <- X%*%theta
      eta_max <- apply(eta,1,max)
      eta_centered <- eta-eta_max
      exp_eta <- exp(eta_centered)
      dn <- exp(-eta_max)+rowSums(exp_eta)
      P <- exp_eta/dn

      u <- rep(0,(m+q*(k+m)-q^2))

      for (i in 1:n) {
        x_i <- X[i,-1]
        D2 <- t(Phi)%*%t(Lambda)%*%diag(x_i)
        Phi_min <- Phi[-1,,drop=FALSE]
        xgamma <- matrix(t(x_i)%*%gamma,nrow=1)
        D3 <- kronecker(t(Phi_min),xgamma)
        xgamlam <- xgamma%*%Lambda
        S_phi <- rbind(matrix(0,q,m-q),diag(1,(m-q)))
        D4 <- kronecker(S_phi,xgamlam)
        D <- cbind(diag(1,m),D2,D3,D4)
        u <- u+t(D)%*%(Z[i,]-P[i,])
      }
      return(c(u))
    }

    F <- function(x) {
      if (k>1) {
        gamma <- diag(x[(m+1):(m+k)])
      } else if (k==1) {
        gamma <- x[m+1]
      }

      if (q==1) {
        Lambda <- rep(1,k)
      } else if (q>1) {
        Lambda <- matrix(c(rep(1,k),x[(m+k+1):(m+k*q)]),nrow=k)
      }
      U <- upper.tri(diag(q),diag=TRUE)*1
      if (q==m) {
        Phi <- U
      }
      if (q<m && k>1) {
        Phi <- cbind(U,matrix(x[(m+k*q+1):(m+q*(k+m)-q^2)],nrow=q))
      } else if (q==1 && k==1) {
        Phi <- matrix(c(U,x[(m+2):(2*m)]),nrow=1)
      }

      if (k==1) {
        theta <- matrix(c(x[1:m],c(gamma)*Phi),nrow=2,byrow=TRUE)
      } else if (k>1) {
        theta <- rbind(x[1:m],gamma%*%Lambda%*%Phi)
      }

      eta <- X%*%theta
      eta_max <- apply(eta,1,max)
      eta_centered <- eta-eta_max
      exp_eta <- exp(eta_centered)
      dn <- exp(-eta_max)+rowSums(exp_eta)
      P <- exp_eta/dn

      I <- matrix(0,(m+q*(k+m)-q^2),(m+q*(k+m)-q^2))

      for (i in 1:n) {
        x_i <- X[i,-1]
        D2 <- t(Phi)%*%t(Lambda)%*%diag(x_i)
        Phi_min <- Phi[-1,,drop=FALSE]
        xgamma <- matrix(t(x_i)%*%gamma,nrow=1)
        D3 <- kronecker(t(Phi_min),xgamma)
        xgamlam <- xgamma%*%Lambda
        S_phi <- rbind(matrix(0,q,m-q),diag(1,(m-q)))
        D4 <- kronecker(S_phi,xgamlam)
        D <- cbind(diag(1,m),D2,D3,D4)
        W <- diag(P[i,])-P[i,]%*%t(P[i,])
        I <- I+t(D)%*%W%*%D
      }
      return(I)
    }

    logLik<-function(x){
      if (k>1) {
        gamma <- diag(x[(m+1):(m+k)])
      } else if (k==1) {
        gamma <- x[m+1]
      }
      if (q==1) {
        Lambda <- rep(1,k)
      } else if (q>1) {
        Lambda <- matrix(c(rep(1,k),x[(m+k+1):(m+k*q)]),nrow=k)
      }
      U <- upper.tri(diag(q),diag=TRUE)*1
      if (q==m) {
        Phi <- U
      } else if (q<m && k>1) {
        Phi <- cbind(U,matrix(x[(m+k*q+1):(m+q*(k+m)-q^2)],nrow=q))
      } else if (q==1 && k==1) {
        Phi <- matrix(c(U,x[(m+2):(2*m)]),nrow=1)
      }
      if (k==1) {
        theta <- matrix(c(x[1:m],c(gamma)*Phi),nrow=2,byrow=TRUE)
      } else if (k>1) {
        theta <- rbind(x[1:m],gamma%*%Lambda%*%Phi)
      }
      eta <- X%*%theta
      sum(eta*Z)-sum(log(1+rowSums(exp(eta))))
    }

    if (q==1) {
      A<-matrix(0,nrow=m-1,ncol=(2*m+k-1))
      A[1,(m+k+1)]<-1
      for (i in 2:(m-1)) {
        A[i,(m+k+i-1)]<--1
        A[i,(m+k+i)]<-1
      }
      B<-c(-1,rep(0,m-2))
    } else if (q>1 && q<m) {
      A<-diag(m+q*(k+m)-q^2)[(m+k+1):(m+q*(k+m)-q^2),]
      #A<-matrix(0,nrow=(k*(q-1)+q*(m-q)),ncol=(m+q*(k+m)-q^2))
      #for (i in 1:(k*(q-1)+q)) {
      #  A[i,m+k+i]<-1
      #}
      if ((m-q)>1) {
        for (i in (k*(q-1)+q+1):(k*(q-1)+q*(m-q))) {
          A[i,m+k+i-q]<--1
          #A[i,m+k+i]<-1
        }
      }
      #A[k*(q-1)+1,m+k*q+1]<-1
      B<-c(rep(0,k*(q-1)),rep(-1,q),rep(0,q*(m-q)-q))
    } else if (q==m) {
      A<-diag(m+q*(k+m)-q^2)[(m+k+1):(m+k*q),]
      B<-rep(0,k*(q-1))
    }

    #fit<-maxLik::maxLik(logLik,grad=grad,hess=H,start=sv,constraints=list(ineqA=A,ineqB=B))
    fit<-stats::constrOptim(sv,logLik,grad,A,B,control=list(fnscale=-1,reltol=1e-14),outer.eps=1e-14)

    alpha<-fit$par[1:m]
    gamma <- fit$par[(m+1):(m+k)]
    if (k>1) {
      Gamma<-diag(gamma)
    } else if (k==1) {
      Gamma<-gamma
    }
    if (q==1) {
      Phi<-matrix(c(1,fit$par[(m+k+1):(2*m+k-1)]),nrow=1)
      Lambda<-matrix(rep(1,k),nrow=k)
    }
    if (q>1) {
      Lambda<-matrix(c(rep(1,k),fit$par[(m+k+1):(m+q*k)]),nrow=k)
      Phi1<-upper.tri(diag(q),diag=TRUE)*1
    }
    if (q>1 && q<m) {
      Phi<-cbind(Phi1,matrix(fit$par[(m+q*k+1):(m+q*(k+m)-q^2)],nrow=q))
    }
    if (q==m) {
      Phi<-Phi1
    }
    Beta<-Gamma%*%Lambda%*%Phi
    JC<-jacob(alpha,gamma,Lambda,Phi)
    x<-c(rbind(alpha,Gamma%*%Lambda%*%Phi))
    #se<-sqrt(diag(JC%*%solve(F(fit$par))%*%t(JC)))
    se<-svd_se(F(fit$par),JC)$se
    loglik <- logLik(fit$par)
  } else if (q==0) {
    x <- sva
    P <- exp(sva)/(1+sum(exp(sva)))
    I <- n*(diag(P)-P%*%t(P))
    se <- sqrt(diag(solve(I)))
    loglik <- ll0
  }
  z <- x/se
  npars <- m+q*(k+m)-q^2
  aic <- 2*npars-2*loglik
  bic <- npars*log(n)-2*loglik
  #C <- round(cbind(x,se,z,1-stats::pnorm(abs(z))),3)
  res <- round(cbind(x,se,z,1-stats::pnorm(abs(z))),3)
  rownames(res) <- paste(colnames(X),rep(1:m,rep(k+1,m)),sep=':')
  colnames(res) <- c('Estimate','Std. Error','z value','Pr(>|z|)')
  #if (q==0) {
  #  P <- matrix(rep(P,n),nrow=n,byrow=T)
  #}
  #fitv <- cbind(1-rowSums(P),P)
  #res <- cbind(1-rowSums(Z),Z)-fitv
  #pcat <- apply(fitv,1,which.max)
  if (q>0) {
    CS <- 1-exp(2*(ll0-loglik)/n)
    NK <- CS/(1-exp(2*ll0/n))
    MF <- 1-(loglik/ll0)
    out <- list(pars=fit$par,npars=npars,A=A,B=B,
                alpha=alpha,Beta=Beta,Gamma=Gamma,Phi=Phi,loglik=loglik,
                npars=npars,AIC=aic,BiC=bic,
                #fitted.values=fitv,residuals=res,pred.cat=pcat,
                obs.cat=y,CS=CS,NK=NK,MF=MF)
  } else if (q==0) {
    out <- list(pars=sva,npars=npars,loglik=loglik,npars=npars,AIC=aic,BIC=bic)
     #fitted.values=fitv,residuals=res,pred.cat=pcat,obs.cat=y)
  }

  if (!is.null(preds)) {
  #  cat(paste('Number of iterations:',it),'\n')
  #  cat(paste('Difference:',round(d,8)),'\n\n')
  #  if (!is.null(ord)) {
    cat(paste('Rank',q,'stereotype ordered'),'\n\n')
  #  }
  } else {
    cat(paste('Null model'),'\n\n')
  }
  cat(paste('Number of observations:',n,'\n\n'))
  cat('Coefficients:\n\n')
  print(res)
  invisible(out)
}



