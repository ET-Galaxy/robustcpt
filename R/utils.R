# Created by: Edwin Tang
# Created on: 12/10/2025

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

#' @title RUME estimator
#
#' @description Robust Univariate Mean Estimation, suggested by Prasad et al. (2021)
#'
#' @param X Univariate dataset
#' @param epsilon Contamination level, a number in (0,1) 0
#' @param delta Confidence level, a number in (0,1)
#' @returns The robust estimate of mean
#'
#' @export
rume<-function(X,epsilon=0,delta=0.05){
  if (length(X)%%2!=0){
    stop("X must be of even length.")
  }
  n<-length(X)/2
  # Sample splitting
  X1<-sort(X[1:n])
  X2<-X[(n+1):(2*n)]

  # Find shortest interval in X1
  vareps<-max(epsilon, log(1/delta)/n)
  datapts<-floor(n*(1-2*vareps-2*sqrt(vareps*log(1/delta)/n)-log(1/delta)/n))
  if (2*vareps+2*sqrt(vareps*log(1/delta)/n)+log(1/delta)/n>=1/2){
    stop("Not enough sample size to achieve required confidence level.")
  }
  candidates<-c()
  for (i in 1:(n-datapts)){
    candidates[i]<-X1[i+datapts]-X1[i]
  }
  min_index<-which.min(candidates)
  interval<-c(X1[min_index], X1[min_index+datapts])

  # Find mean in X2
  filtered_X2<-X2[X2<=interval[2] & X2>=interval[1]]  # Could be empty?
  return(mean(filtered_X2))
}

#' @title t-distribution contaminated with normal samples
#'
#' @description Generate samples from a t-distribution contaminated with normal
#' distribution with large variance
#'
#' @param n The number of samples to generate. A positive integer.
#' @param df The degrees of freedom of the t-distribution. A positive number, defaults to 3.
#' @param epsilon Contamination level for the Normal(0,25) distribution. A number in (0,1)
#' @param mu Mean of the inlier distribution.
#' @returns Samples from a t-distribution contaminated with normal
#'
#' @export
contaminated_sample_t<-function(n,df=3,epsilon=0,mu=0){
  all_indices<-1:n
  final_sample<-rt(n,df)+mu
  contaminated_indices<-all_indices[rbinom(n,1,epsilon)==1]
  for (i in contaminated_indices){
    final_sample[i]<-rnorm(1,0,100)
  }
  return(final_sample)
}

#' @title Generate samples from a change point model
#'
#' @description Generate samples from a change point model with user-specified
#' contamination level and underlying data generating mechanism.
#'
#' @param n The number of samples to generate. A positive integer.
#' @param mechanism A function for the sample generating mechanism with two inputs, n and epsilon.
#' Defaults to \code{contaminated_sample_t}.
#' @param cpt Location of change point. A positive integer less than n or \code{NA}, which
#' means no change point. Defaults to NA
#' @param kappa Change in mean at the change point.
#' @returns Samples from a t-distribution contaminated with normal distributed data.
#'
#' @export
change_point_model<-function(n, mechanism=contaminated_sample_t,cpt=NA,kappa=1){
  if (is.na(cpt)){
    # No change point: generate all n points using mechanism
    final_sample<-mechanism(n=n)
  } else {
    # Generate first segment (before change point)
    final_sample<-mechanism(n=cpt)

    # Generate second segment (after change point) with shift mu = kappa
    final_sample[(cpt+1):n]<-mechanism(n=n-cpt, mu=kappa)
  }
  return(final_sample)
}


#' @title rumedian changepoint algorithm
#'
#' @description Runs the rumedian algorithm on \code{online_data}.
#'
#' @param online_data The complete dataset to analyse. A vector of numbers.
#' @param epsilon Contamination level. A number in (0,1).
#' @param alpha Type I error. A number in (0,1).
#' @returns Location where change point is detected
#'
#' @export
rumedian<-function(online_data, epsilon=0, alpha=0.1){
  if (epsilon>0.1){
    stop("epsilon>0.1 is not supported.")
  }
  for (t in 2:n){
    delta_t=(8*alpha)/(3*t^3-3*t)
    h_t<-ceiling(20*log(1/delta_t))
    for (s in 1:floor(t/2)){
      if (s%%2==0 && s>=h_t){
        # Use RUME
        vareps<-max(epsilon, 2*log(1/delta_t)/s)
        diff_rume<-abs(rume(online_data[(t-s+1):t], epsilon=epsilon)-rume(online_data[1:s], epsilon=epsilon))
        zeta<-2*sigma*1.59*max(vareps^(1-1/v), sqrt(2/s*log(1/delta_t)))
        if (diff_rume>zeta){
          return(list(method="RUME",subsample=s,location=t))
        }
      } else {
        # Use median
        diff_median<-abs(median(online_data[(t-s+1):t])-median(online_data[1:s]))
        chi<-0.13*2*2.25*(0.5*exp(-1)*(delta_t/2)^(2/s)-epsilon)^(-1/v)
        if (diff_median>chi){
          return(list(method="median",subsample=s,location=t))
        }
      }
    }
  }
  return(list(method="no changepoint",subsample=s+1,location=t+1))
}
