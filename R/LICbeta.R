#' Calculate the LIC estimator for linear regression
#'
#' This function estimates the coefficients of a linear regression model using a design matrix `X` and a response vector `Y`. It implements an A-optimal and D-optimal design criteria to choose optimal subsets of observations.
#'
#' @param X The observation matrix (n x p)
#' @param Y The response vector (n x 1)
#' @param alpha The significance level for computing confidence intervals
#' @param K The number of subsets
#' @param nk The number of observations per subset
#'
#' @return A list containing:
#' \item{E5}{The LIC estimator for linear regression.}
#' @export
#' @importFrom stats qt
LICbeta=function (X, Y, alpha, K, nk)
{
  n = nrow(X)
  p = ncol(X)
  N = L1 = c(1:K)
  Rm = matrix(rep(0, nk * K), ncol = K)
  mr = matrix(rep(0, K * nk), ncol = nk)
  for (i in 1:K) {
    mr[i, ] = sample(1:n, nk, replace = T)
    r = matrix(c(1:nk, mr[i, ]), ncol = nk, byrow = T)
    Rm[, i] = r[2, ]
    R = matrix(rep(0, nk * n), ncol = n)
    R[t(r)] = 1
    X1 = R %*% X
    Y1 = R %*% Y
    Hr = X1 %*% solve(crossprod(X1)) %*% t(X1)
    I1 = diag(rep(1, nk))
    SX = (t(Y1) %*% (I1 - Hr) %*% Y1)/(nk - p)
    SY = sqrt(t(Y1) %*% (I1 - Hr) %*% Y1)/(nk - p)
    C1 = sum(diag(X1 %*% solve(crossprod(X1)) %*% t(X1)))/nk
    L1[i] = 2 * SY * C1 * qt(1 - alpha/2, nk - p)
    N[i] = det(t(X1) %*% X1)
  }
  opt1 = Rm[, which.min(L1)]
  opt2 = Rm[, which.max(N)]
  opt = intersect(opt1, opt2)
  Yopt = Y[opt]
  Xopt = X[opt, ]
  I = diag(rep(1, length(Yopt)))
  betalic= solve(t(Xopt) %*% Xopt) %*% t(Xopt) %*% Yopt
  beta
  E5=sum((betalic-beta)^2)/abs(nk-p)
  return(E5)
}
