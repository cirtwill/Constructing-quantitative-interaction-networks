# priordata is a list of interaction probabilies used to build the prior
# Returns prior distribution when n=0 and k=0. 
# Otherwise, returns posterior distribution
calculate_parameters<-function(priordata,n,k){
  library(MASS)
  # Calculate prior parameters
  pars=fitdistr(x=priordata,"beta",start=list(shape1=1,shape2=1),lower=c(0,0))$estimate
  # The lower=c(0,0) prevents R from fitting invalid (negative) parameters
  alpha=pars[[1]]
  beta=pars[[2]]
  # Update parameters with data. If n=0 and k=0, no change.
  alpha_prime=alpha+k
  beta_prime=beta+n-k
  pars2=c(alpha_prime,beta_prime)
  return(pars2)
}


# Returns maximum likelihood estimator for interaction probability
# Based on parameters returned by calculate_parameters
calculate_mean_MLE<-function(priordata,n,k){
  library(MASS)
  pars=fitdistr(x=priordata,"beta",start=list(shape1=1,shape2=1),lower=c(0,0))$estimate
  # The lower=c(0,0) prevents R from fitting invalid (negative) parameters
  alpha=pars[[1]]
  beta=pars[[2]]

  numerator=alpha+k
  denominator=alpha+beta+n
  MLE=numerator/denominator
  return(MLE)
}

# Returns mean and variance of an estimated interaction probability
# Based on parameters returned by calculate_parameters
calculate_distribution<-function(pars){
  alpha=pars[[1]]
  beta=pars[[2]]
  # Mean
  mu_num=alpha
  mu_den=alpha+beta
  mu=mu_num/mu_den
  # Variance
  sig_num=alpha*beta
  den1=alpha+beta
  den2=den1**2
  sig_den=den1*den2
  sigma2=sig_num/sig_den

  return(c(mu,sigma2))
}


# Returns a credible interval between two probability bounds 
# (e.g., 0.025, 0.975 for a 95% CI)
# Based on parameters returned by calculate_parameters
credible_interval<-function(pars,p_lower,p_upper){
  alpha=pars[[1]]
  beta=pars[[2]]
  lowCI=qbeta(p=p_lower,shape1=alpha,shape2=beta)
  highCI=qbeta(p=p_upper,shape1=alpha,shape2=beta)
  return(c(lowCI,highCI))
}


# Returns the number of samples required to reach a given level
# of confidence that the probability of interaction is below a
# given threshold.
# Based on parameters returned by calculate_parameters
samples_for_threshold<-function(threshold,confidence,pars){
  alpha=pars[[1]]
  beta=pars[[2]]
  n=seq(0,1000,1)
  k=0
  cdf=pbeta(threshold,shape1=alpha,shape2=beta+n)
  samples=length(which(cdf<confidence))
  return(samples)
}


