#' Caculate the estimator of beta on the COR
#'
#' @param K is the number of subsets
#' @param nk is the length of subsets
#' @param alpha is the significance level
#' @param X is the observation matrix
#' @param y is the response vector
#' 
#' @return betaC
#' @export

#' @examples 
#'  p=6;n=1000;K=2;nk=200;alpha=0.05;sigma=1
#'  e=rnorm(n,0,sigma); beta=c(sort(c(runif(p,0,1)))); 
#'  data=c(rnorm(n*p,5,10));X=matrix(data, ncol=p);
#'  y=X%*%beta+e;
#'  beta_cor(K=K,nk=nk,alpha=alpha,X=X,y=y)

beta_cor=function(K=K,nk=nk,alpha=alpha,X=X,y=y){
  n=nrow(X);p=ncol(X)
  M=N=W=c(rep(1,K))
  R=matrix(rep(0,n*nk), ncol=n); Io=matrix(rep(0,nk*K), ncol=nk);
  mr=matrix(rep(0,K*nk),ncol=nk)
  for (i in 1:K){ 
    mr[i,]=sample(1:n,nk,replace=T);
    r=matrix(c(1:nk,mr[i,]),ncol=nk,byrow=T);
    R[t(r)]=1
    Io[i,]=r[2,]
    X1=R%*%X;y1=R%*%y;
    ux=solve(crossprod(X1)) 
    W[i]= sum(diag(t(ux)%*% ux))
    M[i]=  det(X1%*%t(X1))
    N[i]=t(y1)%*% y1
  }
  int=intersect(intersect(Io[which.min(W),],Io[which.min(M),]),Io[which.min(N),])  
  lWMN=length(intersect(intersect(Io[which.min(W),],Io[which.max(M),]),Io[which.min(N),]))
  Xcor= X[intersect(intersect(Io[which.min(W),],Io[which.max(M),]),Io[which.min(N),]),];
  ycor= y[intersect(intersect(Io[which.min(W),],Io[which.max(M),]),Io[which.min(N),])]  
  betaC=solve(crossprod(Xcor))%*%t(Xcor)%*% ycor
  return(list(betaC=betaC))
}
