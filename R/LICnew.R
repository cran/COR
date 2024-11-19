#' Calculate the LIC estimator based on A-optimal and D-optimal criterion
#'
#' @param X A matrix of observations (design matrix) with size n x p
#' @param Y A vector of responses with length n
#' @param alpha The significance level for confidence intervals
#' @param K The number of subsets to consider
#' @param nk The size of each subset
#'
#' @return A list containing:
#' \item{E5}{The LIC estimator based on A-optimal and D-optimal criterion.}
#' @export
#'
#' @examples
#' p = 6; n = 1000; K = 2; nk = 200; alpha = 0.05; sigma = 1
#' e = rnorm(n, 0, sigma); beta = c(sort(c(runif(p, 0, 1))));
#' data = c(rnorm(n * p, 5, 10)); X = matrix(data, ncol = p);
#' Y = X %*% beta + e;
#' LICnew(X = X, Y = Y, alpha = alpha, K = K, nk = nk)
#' @importFrom stats qt
LICnew=function (X, Y, alpha, K, nk)
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
  E5=t(Yopt) %*% (I - Xopt %*% solve(crossprod(Xopt)) %*% t(Xopt)) %*% Yopt/(length(Yopt) - p)
  return(E5)
}
