library(MASS)
library(numDeriv)

dat<-read.table('eyedisease.dat',header=TRUE)
n<-dim(dat)[1]
X<-as.matrix(dat[,1:4])
k<-dim(X)[2]
m<-length(unique(dat$erl))-1

U<-cbind(rep(1,n),X)

Z<-matrix(,nrow=n,ncol=m)
Z[,1]<-ifelse(dat$erl==1,1,0)
Z[,2]<-ifelse(dat$erl==2,1,0)
Z[,3]<-ifelse(dat$erl==3,1,0)

if (m==1){ns<-sum(Z)} else if (m>1){ns<-colSums(Z)}

########################################################################

logl_fixed<-function(y){
  if (m==1){
    Dt<-rep(1,k)
    Th<-Dt
  } else if (m>1){
    Dt<-matrix(c(rep(1,k),y[(m+k+1):(m+m*k)]^2),nrow=k)
    Th<-t(apply(Dt,1,cumsum))
  }
  linpred<-U%*%rbind(y[1:m],y[(m+1):(m+k)]*Th)
  if (is.null(dim(linpred))) linpred<-matrix(linpred,ncol=1)
  row_m<-apply(linpred,1,max)
  row_m<-pmax(0,row_m)
  s<-rowSums(exp(linpred-row_m))
  log_dn<-row_m+log(exp(-row_m)+s)
  if (m==1){
    out<--ns*y[1:m]-t(y[(m+1):(m+k)])%*%rowSums(Th*t(X)%*%Z)+sum(log_dn)
  } else {
    out<--ns%*%y[1:m]-t(y[(m+1):(m+k)])%*%rowSums(Th*t(X)%*%Z)+sum(log_dn)
  }
  as.numeric(out)
}

dl<-function(y){if (m==1){dt<-rep(1,k)} else if (m>1){dt<-matrix(c(rep(1,k),y[(m+k+1):(m+m*k)]),nrow=k)}
if (m==1){Th<-dt} else if (m>1){Th<-t(apply(dt^2,1,cumsum))}
E<-U%*%rbind(y[1:m],y[(m+1):(m+k)]*Th)
R<-t(apply(E,1,function(row){
  M<-max(c(0,row))
  num<-exp(row-M)
  denom<-exp(-M)+sum(num)
  num/denom
}))
da<-colSums(R)-ns
V<-t(X)%*%(R-Z)
dg<-rowSums(V*Th)
if (m==1){return(c(da,dg))} else
if (m>1){dd<-2*y[(m+1):(m+k)]*dt*t(apply(t(apply(t(apply(V,1,rev)),1,cumsum)),1,rev))
return(c(da,dg,dd[,2:m]))}}

sv<-rnorm(m*k+m) # starting values

fit<-optim(sv,logl_fixed,dl,method='BFGS',lower=-Inf,upper=Inf,
           control=list(trace=TRUE,reltol=1e-100,abstol=1e-100,maxit=50000),hessian=F)
dl(fit$par)

HES<-hessian(logl_fixed,fit$par) # package numDeriv
se<-sqrt(diag(solve(HES)))
round(cbind(fit$par,se,1-pnorm(abs(fit$par/se))),3)

if (m==1){dt<-rep(1,k)} else if (m>1){dt<-matrix(c(rep(1,k),fit$par[(m+k+1):(m+m*k)]),nrow=k)}

if (m==1){round(rbind(fit$par[1:m],fit$par[(m+1):(m+k)]*dt),3)} else
if (m>1){round(rbind(fit$par[1:m],fit$par[(m+1):(m+k)]*t(apply(dt^2,1,cumsum))),3)}

# Jacobian Matrix
D<-t(apply(dt^2,1,cumsum))
A<-rbind(diag(1,m),matrix(0,nrow=m*k,ncol=m))
E<-matrix(0,nrow=m,ncol=k)
for (i in 1:m){E<-rbind(E,diag(D[,i]))}
J<-cbind(A,E)
for (s in 2:m){
C<-matrix(0,nrow=(m+(s-1)*k),ncol=k)
for (i in s:m){C<-rbind(C,2*diag(dt[,s]*fit$par[(m+1):(m+k)]))}
J<-cbind(J,C)}

std<-sqrt(diag(J%*%solve(HES)%*%t(J)))
Std<-matrix(rbind(std[1:m],matrix(std[(m+1):(m+m*k)],nrow=k)))

par<-matrix(rbind(fit$par[1:m],fit$par[(m+1):(m+k)]*t(apply(dt^2,1,cumsum))),nrow=m+m*k)

out<-round(cbind(par,Std,1-pnorm(abs(par/Std))),3)
rownames(out)<-paste(c('Intercept',paste('x',1:k,sep='')),rep(1:m,rep(k+1,m)),sep=':')
colnames(out)<-c('Est.','s.e.','p')
out

fit$val

library(nnet)

mnlg<-multinom(erl~1,data=dat) # null model

# Cox and Snell:
CS<-1-exp(2*(fit$val-915.354542)/n)
CS

# Nagelkerke:
NK<-CS/(1-exp(-2*915.354542/n))
NK

# McFadden
MF<-1-(fit$val/915.354542)
MF

