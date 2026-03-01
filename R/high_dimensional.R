
#' @title Spectral_filter subroutine
#
#' @description Canconne et al (2023)
#'
#' @param Y dataset as matrix of n rows and p columns
#' @param gamma2 Tuning parameter
#'
#' @export
spectral_filter <- function(Y, gamma2) {

  n <- nrow(Y)
  p <- ncol(Y)

  w <- rep(1, n)
  t <- 1

  repeat {

    # M(w,S) = Y^T diag(w) Y
    M <- t(Y) %*% (c(w) * Y)

    A <- M - n * diag(p)

    # largest singular value
    sv <- svd(A, nu = 1, nv = 0)
    op_norm <- sv$d[1]

    if (op_norm < 5 * gamma2) break

    v <- sv$u[,1]

    # tau_i = <v, Y_i>^2 * 1[w_i > 0]
    projections <- Y %*% v
    tau <- (projections^2) * (w > 0)

    tau1 <- max(tau)
    if (tau1 == 0) break

    # Filter step
    w <- (1 - tau / tau1) * w

    t <- t + 1
  }

  return(w)
}

#' @title Rowsum_filter subroutine
#
#' @description Canconne et al (2023)
#'
#' @param Y dataset as matrix of n rows and p columns
#' @param w Processed weights
#' @param u contamination level
#'
#' @export
rowsum_filter <- function(Y, w, u) {

  n <- nrow(Y)
  p <- ncol(Y)

  sqrt_w <- sqrt(w)

  # Sum(w,S) = sum sqrt(w_i) Y_i
  S_vec <- colSums(Y*c(sqrt_w))

  # tau_i definition
  inner_terms <- (Y %*% S_vec) * sqrt_w
  tau <- abs(inner_terms - w * p) * (w > 0)

  # Remove top un indices
  num_remove <- floor(u * n)
  ord <- order(tau, decreasing = TRUE)

  w_new <- w
  if (num_remove > 0) {
    w_new[ord[1:num_remove]] <- 0
  }

  return(w_new)
}

#' @title Robust Mean Testing routine
#
#' @description Canconne et al (2023)
#'
#' @param Y dataset as matrix of n rows and p columns
#' @param kappa0 Norm of mean difference to be detected
#' @param delta confidence level
#' @param epsilon contamination level
#' @param C_gamma tuning parameter
#'
#' @export
RobustMeanTest <- function(Y, kappa0,
                           delta,
                           epsilon,
                           C_gamma = 1) {

  # Y: n x p matrix
  Y <- as.matrix(Y)
  n <- nrow(Y)
  p <- ncol(Y)

  if (n < p) {
    stop("Require n >= p")
  }
  # contamination fraction u
  u <- epsilon + 1/n+sqrt(log(1/delta)/(2*n))
  if (u>0.1){
    stop("n is too small. u needs to be less than 0.1")
  }

  # Compute gamma_2
  gamma2 <- C_gamma * (
    u * n * p * log(1/u) +
      sqrt(n * p) * log(2 * p / delta) +
      kappa0 * n
  )

  w <- spectral_filter(Y, gamma2)

  w_prime <- rowsum_filter(Y, w, u)

  sqrt_w <- sqrt(pmax(w_prime, 0))
  Sum_wS <- colSums(c(sqrt_w) * Y)

  test_stat <- abs(sum(Sum_wS^2) - p * sum(w_prime))

  if (test_stat >= 0.1 * kappa0^2 * n^2) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' @title High dimensional t-distribution contaminated with normal samples
#'
#' @description Generate samples from a t-distribution contaminated with normal
#' distribution
#'
#' @param n The number of samples to generate. A positive integer.
#' @param p Dimension
#' @param df The degrees of freedom of the t-distribution. A positive number, defaults to 3.
#' @param epsilon Contamination level for the Normal(0,?) distribution. A number in (0,1).
#' @param mu Mean of the inlier distribution.
#'
#' @importFrom MASS mvrnorm
#' @returns Samples from a t-distribution contaminated with normal
#'
#' @export
rt_hd<-function(n,p=1,df=3,epsilon=0,mu=rep(0,p)){
  all_indices<-1:n
  mean_matrix<-t(matrix(mu, nrow=p, ncol=n))
  final_sample<-matrix(rt(n*p,df), nrow=n, ncol=p)+mean_matrix
  contaminated_indices<-all_indices[rbinom(n,1,epsilon)==1]
  for (i in contaminated_indices){
    final_sample[i]<-mvrnorm(n=1, mu=rep(0,p), Sigma=diag(p))
  }
  return(final_sample)
}


#' @title High dimensional Laplace distribution
#'
#' @description Generate samples from a t-distribution contaminated with normal
#' distribution
#'
#' @param n The number of samples to generate. A positive integer.
#' @param p Dimension
#' @param s Scale parameter of the Laplace distribution. Defaults to 1.
#' @param epsilon Contamination level for the Normal(0,?) distribution. A number in (0,1).
#' @param mu Mean of the inlier distribution.
#'
#' @importFrom MASS mvrnorm
#' @importFrom rmutil rlaplace
#' @returns Samples from a t-distribution contaminated with normal
#'
#' @export
rlaplace_hd<-function(n,p=1,s=1,epsilon=0,mu=rep(0,p)){
  all_indices<-1:n
  mean_matrix<-t(matrix(mu, nrow=p, ncol=n))
  final_sample<-matrix(rlaplace(n*p,m=0,s), nrow=n, ncol=p)+mean_matrix
  contaminated_indices<-all_indices[rbinom(n,1,epsilon)==1]
  for (i in contaminated_indices){
    final_sample[i]<-mvrnorm(n=1, mu=rep(0,p), Sigma=diag(p))
  }
  return(final_sample)
}
