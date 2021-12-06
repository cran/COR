#' Caculate the optimal subset lengths on the COR
#'
#' @param K is the number of subsets
#' @param nk is the length of subsets
#' @param alpha is the significance level
#' @param X is the observation matrix
#' @param y is the response vector
#' 
#' @return  seqL, seqN,lWMN
#' @export

#' @examples 
#'  p=6;n=1000;K=2;nk=200;alpha=0.05;sigma=1
#'  e=rnorm(n,0,sigma); beta=c(sort(c(runif(p,0,1)))); 
#'  data=c(rnorm(n*p,5,10));X=matrix(data, ncol=p);
#'  y=X%*%beta+e;
#'  COR(K=K,nk=nk,alpha=alpha,X=X,y=y)

COR=function(K=K,nk=nk,alpha=alpha,X=X,y=y){
  n=nrow(X);p=ncol(X)
  L=M=N=E=W=c(rep(1,K));I=diag(rep(1,nk));betam=matrix(rep(0,p*K), ncol=K)
  R=matrix(rep(0,n*nk), ncol=n); Io=matrix(rep(0,nk*K), ncol=nk);
  mr=matrix(rep(0,K*nk),ncol=nk)
  for (i in 1:K){ 
    mr[i,]=sample(1:n,nk,replace=T);
    r=matrix(c(1:nk,mr[i,]),ncol=nk,byrow=T);
    R[t(r)]=1
    Io[i,]=r[2,]
    X1=R%*%X;y1=R%*%y;
    ux=solve(crossprod(X1)) 
    sy=sqrt((t(y1)%*%(I-X1%*%solve(crossprod(X1))%*%t(X1))%*%y1)/(length(y1)-p))
    L[i]= sy*sum(sqrt(diag(ux)))*(qt(1-alpha/2, length(y1)-p)-qt(alpha/2, length(y1)-p)) 
    W[i]= sum(diag(t(ux)%*% ux))
    M[i]=  det(X1%*%t(X1))
    N[i]=t(y1)%*% y1
    E[i]=t(y1)%*%(I-X1%*%solve(crossprod(X1)) %*%t(X1))%*% y1/(length(y1)-p)
    betam[,i]=solve(t(X1)%*%X1)%*%t(X1)%*%y1
  }
  seqL=which.min(L);seqN=which.min(N)
  int=intersect(intersect(Io[which.min(W),],Io[which.min(M),]),Io[which.min(N),])
  Xc=X[int,];yc= y[int]; I=diag(rep(1,length(int))) 
  t(yc)%*%(I-Xc%*%solve(crossprod(Xc))%*%t(Xc))%*% yc/(n-p)
  minL=L[which.min(L)]
  minM=M[which.min(M)]  
  minN=N[which.min(N)]
  minE=E[which.min(E)]
  lW=length(Io[which.min(W),])
  lWM=length(intersect(Io[which.min(W),],Io[which.max(M),]))
  lWMN=length(intersect(intersect(Io[which.min(W),],Io[which.max(M),]),Io[which.min(N),]))
  return(list(seqL=seqL,seqN=seqN,lWMN=lWMN))
}
