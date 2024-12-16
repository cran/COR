#' Caculate the estimators of beta on the LEV-opt#'
#' @param K is the number of subsets
#' @param nk is the length of subsets
#' @param X is the observation matrix
#' @param Y is the response vector
#'
#' @return A list containing:
#' \item{betalev}{The estimator of beta on the LEV-opt subset.}
#' \item{betam}{The mean of the beta estimators across all K subsets.}
#' \item{AMSE}{The Average Mean Squared Error (AMSE) for the estimator.}
#' \item{WMSE}{The Weighted Mean Squared Error (WMSE) for the estimator.}
#' \item{MSElevb}{The Mean Squared Error (MSE) of the LEV-opt estimator compared to the true beta.}
#' \item{MSEb}{The Mean Squared Error (MSE) of the mean estimator (betam) compared to the true beta.}
#' \item{MSEyleva}{The Mean Squared Error (MSE) of the LEV-opt estimator on the subset with the maximum hat value (Xleva).}
#' \item{MSEyleviy}{The Mean Squared Error (MSE) of the LEV-opt estimator on the subset with the minimum hat value (Xlevi).}
#' \item{MSEW}{The Mean Squared Error (MSE) of the weighted estimator (Wbeta) compared to the true beta.}
#' \item{MSEw}{The Mean Squared Error (MSE) of the weighted estimator (wbeta) compared to the true beta.}
#' @export
#' @references
#' Guo, G., Song, H. & Zhu, L. The COR criterion for optimal subset selection in distributed estimation. \emph{Statistics and Computing}, 34, 163 (2024). \doi{10.1007/s11222-024-10471-z}
#' @importFrom stats var
beta_LW=function(X,Y,K,nk){
  Y=Y;n=nrow(X);p=ncol(X)
  I=diag(c(rep(1,nk)))
  mh=c(rep(1,K))
  R=matrix(rep(0,n*nk),ncol=n)
  mr=matrix(rep(0,nk*K),ncol=nk)
  Io=matrix(rep(0,nk*K), ncol=nk)
  betal=matrix(rep(0,p*K), ncol=K)
  be=bw=bW=matrix(rep(0,K*p),ncol=ncol(X));v=w=rep(0,K)
  for (i in 1:K){
    mr[i,]=sample(1:n,nk,replace=T);
    r=matrix(c(1:nk,mr[i,]),ncol=nk,byrow=T);
    R[t(r)]=1
    Io[i,]=r[2,]
    Xk=R%*%X; yk=R%*%Y
    mh[i]=sum(diag(Xk%*%solve(t(Xk)%*%Xk)%*%t(Xk)))
    betal[,i]=solve(t(Xk)%*%Xk)%*%t(Xk)%*%yk
    ###################
    #be[i,]=solve(crossprod(Xk))%*%t(Xk)%*%yk
    v[i]=1/var(betal[,i])
    #################
    ykhat=Xk%*%solve(t(Xk)%*%Xk)%*%t(Xk)%*%yk
    sigma2hat=mse=sum(yk-ykhat)^2/(nk-ncol(X)-1)#
    sigmahat=sqrt(sum(yk-ykhat)^2/(nk-ncol(X)-1))#sqrt(mse)
    Dbatahat=diag(solve(t(Xk)%*%Xk))*sigmahat#
    sdbatahat=sqrt(diag(solve(t(Xk)%*%Xk))*sigmahat)#
    wsd=(1/sdbatahat)/sum(1/sdbatahat)
    WDb=(1/Dbatahat)/sum(1/Dbatahat)
    bw[i,]=(solve(t(Xk)%*%Xk)%*%t(Xk)%*%yk)*wsd
    bW[i,]=(solve(t(Xk)%*%Xk)%*%t(Xk)%*%yk)*WDb
  }
  for (i in 1:K){
    w[i]=(v[i])/(sum(v))}
  Xleva= X[Io[which.max(mh),],];yleva= Y[Io[which.max(mh),]]
  Xlevi= X[Io[which.min(mh),],];ylevi= Y[Io[which.min(mh),]]
  betalev=solve(t(Xleva)%*%Xleva)%*%t(Xleva)%*% yleva
  betam=rowMeans(betal)
  wbeta=Wbeta=rep(0,ncol(X))
  betaA=colMeans(be)
  betaW=colMeans(be*w)
  for (j in 1:ncol(X)) {
    wbeta[j]=sum(bw[,j])/K
    Wbeta[j]=sum(bW[,j])/K
  }
  MSEyleva=t(yleva)%*%(I-Xleva%*%solve(crossprod(Xleva))%*%t(Xleva))%*%yleva/(length(yleva)-ncol(X));
  MSEyleviy=t(ylevi)%*%(I-Xlevi%*%solve(crossprod(Xlevi))%*%t(Xlevi))%*%ylevi/(length(ylevi)-ncol(X));
  MSElevb=sum((beta-betalev)^2)/abs(nk-ncol(X))
  MSEb=sum((beta-betam)^2)/abs(nk-ncol(X))
  MSEw=sum((beta-wbeta)^2)/abs(nk-ncol(X));
  MSEW=sum((beta-Wbeta)^2)/abs(nk-ncol(X))
  AMSE=sum((beta-betaA)^2)/abs(nk-ncol(X))
  WMSE=sum((beta-betaW)^2)/abs(nk-ncol(X))
  return(list(betalev=c(betalev),betam=rowMeans(betal),AMSE=AMSE,WMSE=WMSE,MSElevb=MSElevb, MSEb=MSEb,MSEyleva=MSEyleva,MSEyleviy=MSEyleviy,MSEW=MSEW,MSEw=MSEw))
}
