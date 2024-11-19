#' Calculate MSE values for different beta estimation methods
#'
#' @param X The design matrix (observations).
#' @param Y The response vector.
#' @param alpha The significance level.
#' @param K The number of subsets.
#' @param nk The length of subsets (number of observations in each subset).
#'
#' @return A list containing:
#' \item{MSECOR}{The MSE of the COR beta estimator.}
#' \item{MSEAopt}{The MSE of the A-optimal beta estimator.}
#' \item{MSEDopt}{The MSE of the D-optimal beta estimator.}
#' \item{MSElic}{The MSE of the LIC beta estimator.}
#' @export
#'
MSEbeta=function (X, Y, alpha, K, nk)
{
  betaAD_result=beta_AD(K=K,nk=nk,alpha=alpha,X=X,y=Y)
  betaA=betaAD_result$betaA
  betaD=betaAD_result$betaD
  betacor_result=beta_cor(K=K,nk=nk,alpha=alpha,X=X,y=Y)
  betaCOR=betacor_result$betaC
  MSEDopt=sum((betaD-beta)^2)/abs(nk-ncol(X))
  MSEAopt=sum((betaA-beta)^2)/abs(nk-ncol(X))
  MSECOR=sum((betaCOR-beta)^2)/abs(nk-ncol(X))
  MSElic=LICbeta(X=X, Y=Y, alpha=alpha, K=K,nk=nk)
  return(list(MSECOR = MSECOR,MSEAopt = MSEAopt, MSEDopt = MSEDopt,MSElic= MSElic))
}
