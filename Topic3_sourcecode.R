rm(list = ls())
#-------------------------------------------------------------------------------
#' Two-piece (TP)  asymmetric lognormal  hazard function
#' https://github.com/Ewnetu-github/TPAN
#-------------------------------------------------------------------------------
#' @param eta  : location parameter
#' @param phi  : scale parameter
#' @param alpha: index parameter
#' @param f0:  : pdf for the reference distribution e.g. f0=dnorm()
#' @param F0:  : CDF for the reference distribution e.g. F0 = pnorm()
#' @param QF:  : quantile function for the reference distribution. e.g. qnorm() 
#' @param g:   : link function e.g. g=log(t)
#' @param log: : log scale (TRUE or FALSE)
#' @t,q        : these are each a vector of quantile.
#' @tau : vector of probabilities 
#' @return the value of the TP asymmetric lognormal  hazard function
#' @export
#------------------------------------------------------------------------------
#' Basic properties for standard normal distribution 
#------------------------------------------------------------------------------
f0_N <- function(s){dnorm(s, mean = 0,sd = 1)} # density 
F0_N <- function(s){pnorm(s, mean = 0,sd = 1)} # CDF
QF_N <- function(tau){qnorm(tau, mean = 0, sd = 1, 
                            lower.tail = TRUE, log.p = FALSE)} 
#-------------------------------------------------------------------------------
#' Basic properties for standard Laplace distribution 
#-------------------------------------------------------------------------------
f0_La <- function(s){0.5*exp(-abs(s))} # density 
F0_La <- function(s){0.5 + 0.5*sign(s)*(1-exp(-abs(s)))} # CDF
QF_La <- function(beta){-sign(beta-0.5)*log(1-2*abs(beta-0.5))}# quantile 
#-------------------------------------------------------------------------------
#' Basic properties for standard logistic distribution 
#-------------------------------------------------------------------------------
f0_Lo <- function(s){ exp(s)*(1+exp(s))^-2} # density 
F0_Lo <- function(s){1/(1+ exp(-s))} # CDF
QF_Lo <- function(tau){log(tau/(1-tau))}# quantile  
#-------------------------------------------------------------------------------

#Link function
glog <- function(s){log(s)}
glog.inv <- function(s){exp(s)} 
glogit    <-function(s){log(exp(s)-1)}
glogit.inv<-function(s){log(exp(s)+1)}

# How to install old version of QBAsyDist package?
# library(devtools)
# install_version("QBAsyDist", version = "0.1.2", repos = "http://cran.us.r-project.org")
library(QBAsyDist)
library(nlme)# for hessian matrix 

# density function for TPA  family 
dtpa <- function(t, eta, phi, alpha, f0, g, log = FALSE){
  g.prime <- Deriv::Deriv(g)
  com <- 2*alpha*(1-alpha)*g.prime(t)/phi
  den <-   ifelse(t < eta, (com*f0((1-alpha)*(g(eta)-g(t))/phi)),
                           (com*f0(alpha*(g(t)-g(eta))/phi)))
  if(log) return(log(den)) else return(den)
}
# cumulative distribution  function for TPA  family 
ptpa <- function(q,eta, phi, alpha, F0, g, lower.tail = TRUE){
  p <- ifelse(q < eta, (2*alpha*F0((1 - alpha)*(g(q) - g(eta))/phi)),
                       (2*alpha - 1 + 2*(1 - alpha)*F0(alpha*(g(q) - g(eta))/phi)))
  if(lower.tail) return(p) else return(1-p)
}
# hazard function for TPA  family
htpa <- function(t, eta, phi, alpha, f0, F0, g, log = FALSE){
  dens <- dtpa(t, eta, phi, alpha, f0, g, log = FALSE)
  surv <- ptpa(t, eta, phi, alpha, F0, g, lower.tail = FALSE)
  h <- dens/surv
  if(log) return(log(h)) else return(h)
}
# cumulative hazard function for TPA family 
Htpa <- function(t,eta, phi, alpha, F0, g, log = FALSE){
  cmh <- -log(ptpa(t, eta, phi, alpha, F0, g, lower.tail = FALSE))
  if(log) return(log(cmh)) else return(cmh)
}
# quantile function for TPA family with log-link function 
qtpa <- function(tau, eta, phi, alpha, F0,g, g.inv = NULL, QF =NULL){
   if(is.null(QF)){
     QF <- GoFKernel::inverse(F0, lower = -Inf, upper = Inf)
   }
if(is.null(g.inv)){
  g.inv <- GoFKernel::inverse(g, lower = 0, upper = Inf)
}
    q <- vector()
  for (i in 1:length(tau)) {
  Ca <- ifelse(tau[i] < alpha, (1/(1-alpha))*QF(tau[i]/(2*alpha)),
                                (1/alpha)*QF((1+tau[i]-2*alpha)/(2*(1-alpha))))
   q[i] <-  g.inv(g(eta) + phi*Ca)
  }                       
  return(q)
}
# Random number generation from TPA family with log-link function
rtpa <- function(n, eta, phi, alpha, F0, g, QF){
  u <- runif(n, min = 0, max = 1)
  r <- vector()
  for (i in 1:n) {
    r[i] <-  qtpa(u[i],eta, phi, alpha,F0, g, QF)
  }
  return(r)
}

#----------------------------------------------------------------------------------------
#' lognormal (lnorm) hazard function.
#----------------------------------------------------------------------------------------

hlnorm <- function(t,mu,sigma, log = FALSE){
  lpdf0 <-  dlnorm(t,mu,sigma, log = T)
  ls0 <- plnorm(t,mu,sigma, lower.tail = FALSE, log.p = T)
  val <- lpdf0 - ls0
  if(log) return(val) else return(exp(val))
}
# lognormal (lnorm) cumulative hazard 
chlnorm <- function(t,mu,sigma){
  H0 <- -plnorm(t,mu,sigma, lower.tail = FALSE, log.p = TRUE)
  return(H0)
}
#----------------------------------------------------------------------------------------
#' loglogistic (llogis) hazard function.
#----------------------------------------------------------------------------------------
hllogis <- function(t,mu,sigma, log = FALSE){
  lpdf0 <-  dlogis(log(t),mu,sigma, log = T) - log(t)
  ls0 <- plogis(log(t),mu,sigma, lower.tail = FALSE, log.p = T)
  val <- lpdf0 - ls0
  if(log) return(val) else return(exp(val))
}
#loglogistic (llogis) cumulative hazard 
chllogis <- function(t,mu,sigma){
  H0 <- -plogis(log(t),mu,sigma, lower.tail = FALSE, log.p = TRUE)
  return(H0)
}
#--------------------------------------------------------------------------------------------------------------------------
#' Exponentiated Weibull (EW) hazard function.
#--------------------------------------------------------------------------------------------------------------------------
#' @param rho     : scale parameter
#' @param kappa      : shape parameter
#' @param gamma   : shape parameter
#' @param t       : positive argument
#' @param log: log scale (TRUE or FALSE)
#' @return the value of the EW hazard function
#' @export
# Probability density function
dexpweibull<- function(t, rho, kappa, gamma,log=FALSE){
  log.pdf <-  log(gamma) + (gamma-1)*pweibull(t,scale=rho,shape=kappa,log=TRUE) + 
                            dweibull(t,scale=rho,shape=kappa,log=TRUE)
  ifelse(log, return(log.pdf), return(exp(log.pdf)))
}
# Cumulative distribution function
pexpweibull<- function(t,rho,kappa,gamma,log.p=FALSE){
  log.cdf <- gamma*pweibull(t,scale=rho,shape=kappa,log.p=TRUE)
  ifelse(log.p, return(log.cdf), return(exp(log.cdf)))
}  
# hazard function for EW
hexpweibull <- function(t,rho,kappa,gamma,log=FALSE){
  log.pdf <-  log(gamma) + (gamma-1)*pweibull(t,scale=rho,shape=kappa,log.p=TRUE) + dweibull(t,scale=rho,shape=kappa,log=TRUE)
  cdf <- exp(gamma*pweibull(t,scale=rho,shape=kappa,log.p=TRUE) )
  log.h <- log.pdf - log(1-cdf)
  ifelse(log, return(log.h), return(exp(log.h)))
}
# cumulative hazard function for EW 
chexpweibull <- function(t,rho,kappa,gamma,log.p=FALSE){
  cdf <- exp(gamma*pweibull(t,scale=rho,shape=kappa,log.p=TRUE) )
  return(-log(1-cdf))
}
# Quantile function for EW
qexpweibull<- function(p,rho,kappa,gamma){
  quant <-  qweibull(p^(1/gamma),scale=rho,shape=kappa)
  return(quant)
} 
#-------------------------------------------------------------------------------
#' Log likelihood and MLE for the parametric hazard based  models.
#' Baseline hazards: 
#-------------------------------------------------------------------------------
#' @param init  : initial points for optimization
#' @param x1    : design matrix for  time scale  covariates (nxd1)
#' @param x2    : design matrix for   hazard  scale covariates (nxd2)
#' @param y     : observed survival times
#' @param delta : censoring indicator or event status 
#' @param dist  : baseline distribution 
#' @param method: "nlminb" or a method from "optim"
#' @param maxit: The maximum number of iterations. Defaults to 1000
#' @return a list containing the output of the optimization (OPT) and 
#'     the negative log-likelihood function (nloglik)
#' @export

#---------------------------------------------------------------------
# HH model: for parametric and nonparametric risk function, g(t)=glog
# semiparametric risk function is not yet included 
#---------------------------------------------------------------------
parHH_mle <- function(init, y, delta, x1, x2, dist, psi2.hat=NULL, 
                       method = "Nelder-Mead", maxit = 10000,...){
  # Required variables
  y <- as.vector(y)
  x1 <- as.matrix(x1)
  x2 <- as.matrix(x2)
  d1 <- dim(x1)[2]
  d2 <- dim(x2)[2]
  
  # TPA log-normal  HH model  
  if(dist == "tpalnorm"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      # baseline parameters 
      eta0 <- exp(par[1]); phi0 <- exp(par[2]);  alpha0 <- expit(par[3]) 
      beta1 <- par[4:(3+d1)] # vector of regression parameters 
      psi.x1 <- as.vector(x1%*%beta1)
      exp.psi.x1 <- exp(psi.x1)
      # profiling out  psi(x2) if it is defined nonparametrically 
      if(is.null(psi2.hat)){
        beta2 <- par[(4+d1):(3+d1+d2)]
        psi.x2 <- as.vector(x2%*%beta2)
      } else{psi.x2 <- psi2.hat }
      
      exp.psi.dif <- as.vector(exp(- psi.x1  + psi.x2 ))
      
      lhaz0 <- htpa(y*exp.psi.x1, eta0, phi0, alpha0, f0_N, F0_N,
                    glog, log = TRUE) + psi.x2
      chaz0 <- Htpa(y*exp.psi.x1, eta0, phi0, alpha0, F0_N, 
                    glog, log = FALSE)*exp.psi.dif
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0, na.rm = T) + sum(chaz0, na.rm = T)#  nloglik
      return(val) 
    }
  }
  # TPA log-logistic HH model  
  if(dist == "tpallogis"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      eta0 <- exp(par[1]); phi0 <- exp(par[2]);  alpha0 <- expit(par[3]) 
      beta1 <- par[4:(3+d1)] # vector of regression parameters 
      psi.x1 <- as.vector(x1%*%beta1)
      exp.psi.x1 <- exp(psi.x1)
      # profiling out  psi(x2) if it is defined nonparametrically 
      if(is.null(psi2.hat)){
        beta2 <- par[(4+d1):(3+d1+d2)]
        psi.x2 <- as.vector(x2%*%beta2)
      } else{psi.x2 <- psi2.hat }
      
      exp.psi.dif <- as.vector(exp(- psi.x1  + psi.x2 ))
      
      lhaz0 <- htpa(y*exp.psi.x1, eta0, phi0, alpha0, f0_Lo, F0_Lo, 
                    glog, log = TRUE) + psi.x2
      chaz0 <- Htpa(y*exp.psi.x1, eta0, phi0, alpha0, F0_Lo, 
                    glog, log = FALSE)*exp.psi.dif
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0, na.rm = T) + sum(chaz0, na.rm = T)#  nloglik
      return(val) 
    }
  }
  
  # TPA log-Laplace HH model 
  if(dist == "tpalLap"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      eta0 <- exp(par[1]); phi0 <- exp(par[2]);  alpha0 <- expit(par[3]) 
      beta1 <- par[4:(3+d1)] # vector of regression parameters 
      
      psi.x1 <- as.vector(x1%*%beta1)
      exp.psi.x1 <- exp(psi.x1)
      # profiling out  psi(x2) if it is defined nonparametrically 
      if(is.null(psi2.hat)){
        beta2 <- par[(4+d1):(3+d1+d2)]
        psi.x2 <- as.vector(x2%*%beta2)
      } else{psi.x2 <- psi2.hat }
      
      exp.psi.dif <- as.vector(exp(- psi.x1  + psi.x2 ))
      
      lhaz0 <- htpa(y*exp.psi.x1, eta0, phi0, alpha0, f0_La, F0_La, 
                    glog, log = TRUE) + psi.x2
      chaz0 <- Htpa(y*exp.psi.x1, eta0, phi0, alpha0, F0_La, 
                    glog, log = FALSE)*exp.psi.dif
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0, na.rm = T) + sum(chaz0, na.rm = T)#  nloglik
      return(val) 
    }
  }
  # EW  HH Model
  if(dist == "EW"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      rh0 <- exp(par[1]); kap0 <- exp(par[2]);  gam0 <- exp(par[3])
      beta1 <- par[4:(3+d1)] # vector of regression parameters 
      psi.x1 <- as.vector(x1%*%beta1)
      exp.psi.x1 <- exp(psi.x1)
      # profiling out  psi(x2) if it is defined nonparametrically 
      if(is.null(psi2.hat)){
        beta2 <- par[(4+d1):(3+d1+d2)]
        psi.x2 <- as.vector(x2%*%beta2)
      } else{psi.x2 <- psi2.hat }
      
      exp.psi.dif <- as.vector(exp(- psi.x1  + psi.x2 ))
      
      lhaz0 <- hexpweibull(y*exp.psi.x1, rh0, kap0, gam0, log = TRUE) + psi.x2
      chaz0 <- chexpweibull(y*exp.psi.x1,rh0, kap0, gam0)*exp.psi.dif
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0, na.rm = T) + sum(chaz0, na.rm = T)#  nloglik
      return(val) 
    } 
  }
  if(method != "nlminb") OPT <- optim(init,nlog.lik,control=list(maxit=maxit),
                                      method = method)
  if(method == "nlminb") OPT <- nlminb(init,nlog.lik,control=list(iter.max=maxit))
  OUT <- list(OPT = OPT, nloglik = nlog.lik)
  return(OUT)
}

#------------------------------------------------------------------------------
 # GH Model:  with parametric risk function 
#------------------------------------------------------------------------------
parGH_mle <- function(init, y, delta, x, dist, method = "Nelder-Mead",
                      maxit = 10000,...){
  # Required variables
  y <- as.vector(y)
  delta <- as.vector(delta)
  x <- as.matrix(x)
  d <- dim(x)[2]
  
  # TPA log-normal  HH model  
  if(dist == "tpalnorm"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      eta0 <- exp(par[1]); phi0 <- exp(par[2]);  alpha0 <- expit(par[3]) 
      beta1 <- par[4:(3+d)] # vector of regression parameters 
      beta2 <- par[(4+d):(3+2*d)]
      
      psi.x1 <- as.vector(x%*%beta1)
      psi.x2 <- as.vector(x%*%beta2)
      
      exp.psi.x1 <- exp(psi.x1)
      exp.psi.dif <- as.vector(exp(-psi.x1  + psi.x2 ))
      
      lhaz0 <- htpa(y*exp.psi.x1, eta0, phi0, alpha0, f0_N, F0_N, glog, log = TRUE) + psi.x2
      chaz0 <- Htpa(y*exp.psi.x1, eta0, phi0, alpha0, F0_N, glog, log = FALSE)*exp.psi.dif
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0, na.rm = T) + sum(chaz0, na.rm = T)#  nloglik
      return(val) 
    }
  }
  # TPA log-logistic GH model  
  if(dist == "tpallogis"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      eta0 <- exp(par[1]); phi0 <- exp(par[2]);  alpha0 <- expit(par[3]) # par. for baseline distribution
      beta1 <- par[4:(3+d)] # vector of regression parameters 
      beta2 <- par[(4+d):(3+2*d)]
      psi.x1 <- as.vector(x%*%beta1)
      psi.x2 <- as.vector(x%*%beta2)
      
      exp.psi.x1 <- exp(psi.x1)
      exp.psi.dif <- as.vector(exp(-psi.x1  + psi.x2 ))
      
      lhaz0 <- htpa(y*exp.psi.x1, eta0, phi0, alpha0, f0_Lo, F0_Lo, glog, log = TRUE) + psi.x2
      chaz0 <- Htpa(y*exp.psi.x1, eta0, phi0, alpha0, F0_Lo, glog, log = FALSE)*exp.psi.dif
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0, na.rm = T) + sum(chaz0, na.rm = T)#  nloglik
      return(val) 
    }
  }
  
  # TPA log-Laplace GH model 
  if(dist == "tpalLap"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      eta0 <- exp(par[1]); phi0 <- exp(par[2]);  alpha0 <- expit(par[3]) # par. for baseline distribution
      beta1 <- par[4:(3+d)] # vector of regression parameters 
      beta2 <- par[(4+d):(3+2*d)]
      
      psi.x1 <- as.vector(x%*%beta1)
      psi.x2 <- as.vector(x%*%beta2)
      
      exp.psi.x1 <- exp(psi.x1)
      exp.psi.dif <- as.vector(exp(-psi.x1  + psi.x2 ))
      
      lhaz0 <- htpa(y*exp.psi.x1, eta0, phi0, alpha0, f0_La, F0_La, glog, log = TRUE) + psi.x2
      chaz0 <- Htpa(y*exp.psi.x1, eta0, phi0, alpha0, F0_La, glog, log = FALSE)*exp.psi.dif
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0, na.rm = T) + sum(chaz0, na.rm = T)#  nloglik
      return(val) 
    }
  }
  # EW  GH Model
  if(dist == "EW"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      rh0 <- exp(par[1]); kap0 <- exp(par[2]);  gam0 <- exp(par[3])
      beta1 <- par[4:(3+d)] # vector of regression parameters 
      beta2 <- par[(4+d):(3+2*d)]
      
      psi.x1 <- as.vector(x%*%beta1)
      psi.x2 <- as.vector(x%*%beta2)
      
      exp.psi.x1 <- exp(psi.x1)
      exp.psi.dif <- as.vector(exp(-psi.x1  + psi.x2 ))
      
      lhaz0 <- hexpweibull(y*exp.psi.x1, rh0, kap0, gam0, log = TRUE) + psi.x2
      chaz0 <- chexpweibull(y*exp.psi.x1, rh0, kap0, gam0)*exp.psi.dif
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0, na.rm = T) + sum(chaz0, na.rm = T)#  nloglik
      return(val) 
    } 
  }
  if(method != "nlminb") OPT <- optim(init,nlog.lik,control=list(maxit=maxit),
                                      method = method)
  if(method == "nlminb") OPT <- nlminb(init,nlog.lik,control=list(iter.max=maxit))
  OUT <- list(OPT = OPT, nloglik = nlog.lik)
  return(OUT)
}

#-------------------------------------------------------------------------------
# AFT hazard model: parametric risk function  
#-------------------------------------------------------------------------------
#' @param init  : initial points for optimization
#' @param x     : design matrix for covariates  (n x d)
#' @param delta     : censoring indicator or event status 
#' @param y     : observed survival times
#' @param dist  : baseline distribution 
#' @param method :  methods from "optim" or "nlminb"
#' @param maxit : the maximum number of iterations. Defaults to 10000
#' @return a list containing the output of the optimization (OPT) and the nloglik
#' @export

parAFT_mle <- function(init, y, delta, x,  dist, method = "Nelder-Mead", 
                       psihat=NULL, maxit = 10000,...){
  # Required variables
  y <- as.vector(y)
  x <- as.matrix(x)
  d <- dim(x)[2]
  if(is.null(psihat)){
    psi.x2  <- as.matrix(rep(0, length(y)))
  } else {psi.x2  <- psihat} 
  # TPA log-normal AFT model
  if(dist == "tpalnorm"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      eta0 <- exp(par[1]); phi0 <- exp(par[2]);  alpha0 <- expit(par[3]) 
      beta <- par[4:(3+d)] # vector of regression parameters 
      psi.x1 <- as.vector(x%*%as.vector(beta))
      psi.x <- psi.x1 + psi.x2 
      
      exp.psi.x <- exp(psi.x)
      
      lhaz0 <- htpa(y*exp.psi.x, eta0, phi0, alpha0, f0_N, F0_N, glog, log = TRUE) + psi.x
      chaz0 <- Htpa(y*exp.psi.x, eta0, phi0, alpha0, F0_N, glog, log = FALSE)
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0, na.rm = T) + sum(chaz0, na.rm = T)#  nloglik
      return(val) 
    }
  }
  # TPA log-logistic AFT model
  if(dist == "tpallogis"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      eta0 <- exp(par[1]); phi0 <- exp(par[2]);  alpha0 <- expit(par[3]) 
      beta <- par[4:(3+d)] # vector of regression parameters 
      psi.x1 <- as.vector(x%*%as.vector(beta))
      psi.x <- psi.x1 + psi.x2 
      exp.psi.x <- exp(psi.x)
      
      lhaz0 <- htpa(y*exp.psi.x, eta0, phi0, alpha0, f0_Lo, F0_Lo, glog,log = TRUE) + psi.x
      chaz0 <- Htpa(y*exp.psi.x, eta0, phi0, alpha0, F0_Lo, glog, log = FALSE)
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0, na.rm = T) + sum(chaz0, na.rm = T)#  nloglik
      return(val) 
    }
  }
  # TPA log-Laplace AFT model
  if(dist == "tpalLap"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      eta0 <- exp(par[1]); phi0 <- exp(par[2]);  alpha0 <- expit(par[3]) # par. for baseline distribution
      beta <- par[4:(3+d)] # vector of regression parameters 
      psi.x1 <- as.vector(x%*%as.vector(beta))
      psi.x <- psi.x1 + psi.x2 
      exp.psi.x <- exp(psi.x)
      
      lhaz0 <- htpa(y*exp.psi.x, eta0, phi0, alpha0, f0_La,F0_La, glog,log = TRUE) + psi.x
      chaz0 <- Htpa(y*exp.psi.x, eta0, phi0, alpha0, F0_La, glog, log = FALSE)
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0, na.rm = T) + sum(chaz0, na.rm = T)#  nloglik
      return(val) 
    }
  }
  # EW  AFT Model
  if(dist == "EW"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      rh0 <- exp(par[1]); kap0 <- exp(par[2]);  gam0 <- exp(par[3])
      beta <- par[4:(3+d)] # vector of regression parameters 
      psi.x <- as.vector(x%*%beta)
      exp.psi.x <- exp(psi.x)
      
      lhaz0 <- hexpweibull(y*exp.psi.x, rh0, kap0, gam0, log = TRUE) + psi.x
      chaz0 <- chexpweibull(y*exp.psi.x, rh0, kap0, gam0)
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  -sum(delta*lhaz0, na.rm = T) + sum(chaz0, na.rm = T)#  nloglik
      return(val) 
    } 
  }
  if(method != "nlminb") OPT <- optim(init,nlog.lik,control=list(maxit=maxit), method = method)
  if(method == "nlminb") OPT <- nlminb(init,nlog.lik,control=list(iter.max=maxit))
  OUT <- list(OPT = OPT, nloglik = nlog.lik)
  return(OUT)
}

#-------------------------------------------------------------------------------
# PH hazard model: can be used for parametric, semiparametric  and nonparametric psi2(x) 
#-------------------------------------------------------------------------------
#' @param init  : initial points for optimization
#' @param x     : design matrix for covariate effects (n x d)
#' @param delta : censoring indicator or event status 
#' @param y     : observed survival times
#' @param dist  : baseline distribution 
#' @param method :  methods from "optim" or "nlminb"
#' @param maxit : the maximum number of iterations. Defaults to 10000
#' @return a list containing the output of the optimization (OPT) and (nloglik)
#' @export

parPH_mle <- function(init, y, delta, x,  dist, psihat=NULL, 
                       method = "Nelder-Mead", maxit = 10000,...){
  
  # Required variables
  y <- as.vector(y)
  x <- as.matrix(x)
  d <- dim(x)[2]
  if(is.null(psihat)){
        psi.x2  <- as.matrix(rep(0, length(y)))
  } else {psi.x2  <- psihat} 
  # TPA log-normal 
  if(dist == "tpalnorm"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      eta0 <- exp(par[1]); phi0 <- exp(par[2]);  alpha0 <- expit(par[3]) 
        beta <- par[4:(3+d)] # vector of regression parameters 
        psi.x1 <- as.vector(x%*%as.vector(beta))
        psi.x <- psi.x1 + psi.x2 
        exp.psi.x <- exp(psi.x)
      
      lhaz0 <- htpa(y, eta0, phi0, alpha0, f0_N, F0_N, glog, log = TRUE) + psi.x
      chaz0 <- Htpa(y, eta0, phi0, alpha0, F0_N, glog, log = FALSE)*exp.psi.x
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0, na.rm = T) + sum(chaz0, na.rm = T)#  nloglik
      return(val) 
    }
  }
  # TPA log-logistic  
  if(dist == "tpallogis"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      eta0 <- exp(par[1]); phi0 <- exp(par[2]);  alpha0 <- expit(par[3]) 
      if(is.null(psihat)){
        beta <- par[4:(3+d)] # vector of regression parameters 
        psi.x <- as.vector(x%*%beta)
        exp.psi.x <- exp(psi.x)
      } else{psi.x <- psihat ; exp.psi.x <- exp(psihat)}
      
      lhaz0 <- htpa(y, eta0, phi0, alpha0, f0_Lo, F0_Lo,glog, log = TRUE) + psi.x
      chaz0 <- Htpa(y, eta0, phi0, alpha0, F0_Lo, glog, log = FALSE)*exp.psi.x
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0, na.rm = T) + sum(chaz0, na.rm = T)#  nloglik
      return(val) 
    }
  }
  # TPA log-Laplace 
  if(dist == "tpalLap"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      eta0 <- exp(par[1]); phi0 <- exp(par[2]);  alpha0 <- expit(par[3]) # par. for baseline distribution
      if(is.null(psihat)){
        beta <- par[4:(3+d)] # vector of regression parameters 
        psi.x <- as.vector(x%*%beta)
        exp.psi.x <- exp(psi.x)
      } else{psi.x <- psihat ; exp.psi.x <- exp(psihat)}
      
      lhaz0 <- htpa(y, eta0, phi0, alpha0, f0_La, F0_La, glog, log = TRUE) + psi.x
      chaz0 <- Htpa(y, eta0, phi0, alpha0, F0_La, glog, log = FALSE)*exp.psi.x
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0, na.rm = T) + sum(chaz0, na.rm = T)#  nloglik
      return(val) 
    }
  }
  # EW
  if(dist == "EW"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      rh0 <- exp(par[1]); kap0 <- exp(par[2]);  gam0 <- exp(par[3])
      if(is.null(psihat)){
        beta <- par[4:(3+d)] # vector of regression parameters 
        psi.x <- as.vector(x%*%beta)
        exp.psi.x <- exp(psi.x)
      } else{psi.x <- psihat ; exp.psi.x <- exp(psihat)}
      
      lhaz0 <- hexpweibull(y, rh0, kap0, gam0, log = TRUE) + psi.x
      chaz0 <- chexpweibull(y,rh0, kap0, gam0)*exp.psi.x
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0, na.rm = T) + sum(chaz0, na.rm = T)#  nloglik
      return(val) 
    } 
  }
  if(method != "nlminb") OPT <- optim(init,nlog.lik,control=list(maxit=maxit), method = method)
  if(method == "nlminb") OPT <- nlminb(init,nlog.lik,control=list(iter.max=maxit))
  OUT <- list(OPT = OPT, nloglik = nlog.lik)
  return(OUT)
}


#-------------------------------------------------------------------------------
# AH   model: parametric risk function  
#-------------------------------------------------------------------------------
#' @param init  : initial points for optimization
#' @param x     : design matrix for covariate effects (n x p)
#' @param delta     : censoring indicator or event status 
#' @param y     : observed survival times
#' @param dist  : baseline distribution (TPA lognormal,TPA loglogistic,TPA logLaplace, EW, lognormal, loglogistc and logLaplace)
#' @param method :  methods from "optim" or "nlminb"
#' @param maxit : the maximum number of iterations. Defaults to 10000
#' @return a list containing the output of the optimization (OPT) and the negative log-likelihood function (nloglik)
#' @export

parAH_mle <- function(init, y, delta, x, dist, psihat=NULL, 
                      method = "Nelder-Mead", maxit = 10000,...){
  # Required variables
  y <- as.vector(y)
  x <- as.matrix(x)
  d <- dim(x)[2]
  if(is.null(psihat)){
    psi.x2  <- as.matrix(rep(0, length(y)))
  } else {psi.x2  <- psihat} 
  if(dist == "tpalnorm"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      eta0 <- exp(par[1]); phi0 <- exp(par[2]);  alpha0 <- expit(par[3]) # par. for baseline distribution
      beta <- par[4:(3+d)] # vector of regression parameters 
      psi.x1 <- as.vector(x%*%beta) 
      psi.x <-  psi.x1 + psi.x2
      exp.psi.x <- exp(psi.x)
      
      lhaz0 <- htpa(y*exp.psi.x, eta0, phi0, alpha0, f0_N, F0_N, glog,  log = TRUE)
      chaz0 <- Htpa(y*exp.psi.x, eta0, phi0, alpha0, F0_N, glog, log = FALSE)*exp(-psi.x)
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0, na.rm = T) + sum(chaz0, na.rm = T)#  nloglik
      return(val) 
    }
  }
  if(dist == "tpallogis"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      eta0 <- exp(par[1]); phi0 <- exp(par[2]);  alpha0 <- expit(par[3]) # par. for baseline distribution
      beta <- par[4:(3+d)] # vector of regression parameters 
      psi.x1 <- as.vector(x%*%beta) 
      psi.x <-  psi.x1 + psi.x2
      exp.psi.x <- exp(psi.x)
      
      lhaz0 <- htpa(y*exp.psi.x, eta0, phi0, alpha0, f0_Lo, F0_Lo, glog,  log = TRUE)
      chaz0 <- Htpa(y*exp.psi.x, eta0, phi0, alpha0, F0_Lo, glog, log = FALSE)*exp(-psi.x)
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0, na.rm = T) + sum(chaz0, na.rm = T)#  nloglik
      return(val) 
    }
  }
  if(dist == "tpalLap"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      eta0 <- exp(par[1]); phi0 <- exp(par[2]);  alpha0 <- expit(par[3]) # par. for baseline distribution
      beta <- par[4:(3+d)] # vector of regression parameters 
      psi.x1 <- as.vector(x%*%beta) 
      psi.x <-  psi.x1 + psi.x2
      exp.psi.x <- exp(psi.x)
      
      lhaz0 <- htpa(y*exp.psi.x, eta0, phi0, alpha0, f0_La, F0_La, glog,  log = TRUE)
      chaz0 <- Htpa(y*exp.psi.x, eta0, phi0, alpha0, F0_La, glog, log = FALSE)*exp(-psi.x)
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0, na.rm = T) + sum(chaz0, na.rm = T)#  nloglik
      return(val) 
    }
  }
  # EW - PH Model
  if(dist == "EW"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      rh0 <- exp(par[1]); kap0 <- exp(par[2]);  gam0 <- exp(par[3])
      beta <- par[4:(3+d)] # vector of regression parameters 
      psi.x1 <- as.vector(x%*%beta) 
      psi.x <-  psi.x1 + psi.x2
      exp.psi.x <- exp(psi.x)
      
      lhaz0 <- hexpweibull(y*exp.psi.x, rh0, kap0, gam0, log = TRUE)
      chaz0 <- chexpweibull(y*exp.psi.x,rh0, kap0, gam0)*exp(-psi.x)
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0, na.rm = T) + sum(chaz0, na.rm = T)#  nloglik
      return(val) 
    } 
  }
  if(method != "nlminb") OPT <- optim(init,nlog.lik,control=list(maxit=maxit), method = method)
  if(method == "nlminb") OPT <- nlminb(init,nlog.lik,control=list(iter.max=maxit))
  OUT <- list(OPT = OPT, nloglik = nlog.lik)
  return(OUT)
}


#--------------------------------------------------------------------------
# profile local log-likelihood and MLE for the Semiparametric Hazard models .
# only HH and PH can be accommodate unspecified covariate effects on the relative hazard
# only AFT and AH can be accommodate unspecified psi1 and psi2 at the same time 
# HH cannot be suited to include  psi1 and psi2 at the same time : one of the two 
# functions should be parametric due to identifiable  issue 
#--------------------------------------------------------------------------
#' @param init  : initial points for optimization
#' @param x1    : design matrix for covariate effects on time scale  (nxp1)
#' @param x2    : design matrix for  nonlinear covariate  effects on relative hazard  (nxp2)
#' @param d     : censoring indicator or event status 
#' @param y     : observed survival times
#' @param x0 : local point (grid value) that the parameter needs to be estimated
#' @param h : bandwidth parameter 
#' @param ktype : type of kernel function. Defaults Epanechnikov (Epa) 
#' @param dist  : baseline distribution (TPA lognormal,TPA loglogistic,TPA logLaplace, EW, lognormal, loglogistc and logLaplace)
#' @param method: methods from "optim" or  "nlminb" 
#' @param maxit: The maximum number of iterations. Defaults to 10000
#' @return a list containing the output of the optimization (OPT) and the negative log-likelihood function (nloglik)
#' @export
#-----------------------------------
# Kernel function 
#-----------------------------------
ker <- function(ktype, x) {
  if(ktype=="uniform"){
    result <- dunif(x, -1, 1)
  }
  if (ktype == "Gaussian") {
    result <- dnorm(x)
    
  }
  else if (ktype == "Epan") {
    result <- ifelse(abs(x) <= 1 , 3/4*(1 - x^2), 0 )
    
  }
  else if (ktype == "biweight") {
    result <- ifelse(abs(x) <= 1 , 15/16*(1 - x^2)^2, 0 )
    
  }
  else if (ktype == "triweight") {
    result <- ifelse(abs(x) <= 1 , 35/32*(1 - x^2)^3, 0 )
  } else if(ktype =="cosine"){
    result <- ifelse(abs(x) <= 1 , pi/4*cos(pi/2*x), 0 )
  }
  return(result)
}
#-----------------------------------------------------------
# HH model: to estimate the nonparametric part at the hazard level only
#-----------------------------------------------------------
#' init: initial values to start the optimization
#' y: vector of observed time 
#' delta: vector of censoring indicator
#' x1: matrix  of time scale covariates 
#' x2: matrix  of hazard scale covariates (specified  + unspecified)
#' x21.index : index of  x21 
#' theta.hat: estimates of the transformed baseline parameters 
#' beta1.hat: estimates of the time scale parametric component 
#' beta2.hat: estimates of the hazard ratio parametric component 
#' h: bandwidth value
#' x0: local point 
#' ktype: kernel function 


nonparHH_mle <- function(init=NULL, y, delta, x1, x2, x21.index  = NULL, 
                         theta.hat, beta1.hat, beta2.hat = NULL, dist, 
                         x0, h, ktype = "Epan", method = "Nelder-Mead", maxit = 10000,...){
  if (is.na(h) || h <= 0L)
    stop("invalid 'bandwidth value.' ")
  if(is.null(colnames(x2)))
    stop("give a column names for the 'covariates'")
  # Required variables
  y <- as.vector(y)
  delta <- as.vector(delta)
  
  x1 <- as.matrix(x1)
  d1 <- dim(x1)[2]
  n <- length(y)

  # profiling out the parametric coefficients 
  if(!is.null(x21.index) | !is.null(beta2.hat)){
    x21 <- as.matrix(x2[, x21.index])
    d21 <- ncol(x21)
    x22 <- as.matrix(x2[,-c(1:d21)])
    d22 <- ncol(x22)
    psi.x21 <- as.vector(x21%*%beta2.hat) 
  } else {
    psi.x21 <- rep(0, n)
    x22 <- as.matrix(x2)
    d22 <- dim(x22)[2]
  }
  # local points 
  newx0 = x22 - t(matrix(rep(x0,d22*n),d22,n))
  ker.h <- as.vector(apply(newx0/h, MARGIN = 2, FUN = ker, ktype = ktype)/h)
  newx <- cbind(1, newx0, newx0^2)# upto second order derivative 
  
  if(is.null(init)){
    fit.poly <- lm(y[delta==1]~poly(newx0[delta==1], degree = 2))
    init <- as.numeric(fit.poly$coefficients)
  }
  # TPA lognormal  PH model  
  if(dist == "tpalnorm"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      eta0 <- exp(theta.hat[1]); phi0 <- exp(theta.hat[2]);  alpha0 <- expit(theta.hat[3]) 
      b2 <- par # vector of local regression coefficients
      psi.x22 <- as.vector(newx%*%b2)
      
      psi.x2 <-  psi.x21 + psi.x22
      psi.x1 <- as.vector(x1%*%beta1.hat)
      exp.psi.x1 <- exp(psi.x1)
      exp.psi.dif <- exp(-psi.x1 + psi.x2)

      lhaz0 <- htpa(y*exp.psi.x1, eta0, phi0, alpha0, f0_N, F0_N, glog, log = TRUE) + psi.x2
      chaz0 <- Htpa(y*exp.psi.x1, eta0, phi0, alpha0, F0_N, glog, log = FALSE)*exp.psi.dif
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0*ker.h, na.rm = T) + sum(chaz0*ker.h, na.rm = T)#  nloglik
      return(val) 
    }
  }
  # TPA loglogistic PH model  
  if(dist == "tpallogis"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      eta0 <- exp(theta.hat[1]); phi0 <- exp(theta.hat[2]);  alpha0 <- expit(theta.hat[3]) 
      b2 <- par # vector of local regression coefficients
      psi.x22 <- as.vector(newx%*%b2)
      
      psi.x2 <-  psi.x21 + psi.x22
      psi.x1 <- as.vector(x1%*%beta1.hat)
      exp.psi.x1 <- exp(psi.x1)
      exp.psi.dif <- exp(-psi.x1 + psi.x2)
      
      lhaz0 <- htpa(y*exp.psi.x1, eta0, phi0, alpha0, f0_Lo, F0_Lo, glog, log = TRUE) + psi.x2
      chaz0 <- Htpa(y*exp.psi.x1, eta0, phi0, alpha0, F0_Lo, glog, log = FALSE)*exp.psi.dif
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0*ker.h, na.rm = T) + sum(chaz0*ker.h, na.rm = T)#  nloglik
      return(val) 
    }
  }
  # EW - Model
  if(dist == "EW"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      par0 <- exp(theta.hat)
      b2 <- par # vector of local regression coefficients  
      psi.x22 <- as.vector(newx%*%b2)
      
      psi.x2 <-  psi.x21 + psi.x22
      psi.x1 <- as.vector(x1%*%beta1.hat)
      exp.psi.x1 <- exp(psi.x1)
      exp.psi.dif <- exp(-psi.x1 + psi.x2)
      
      lhaz0 <- hexpweibull(y*exp.psi.x1, par0[1], par0[2], par0[3], log = TRUE) + psi.x2
      chaz0 <- chexpweibull(y*exp.psi.x1, par0[1], par0[2], par0[3])*exp.psi.dif
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0*ker.h, na.rm = T) + sum(chaz0*ker.h, na.rm = T)#  nloglik
      return(val) 
    } 
  }
  if(method != "nlminb") OPT <- optim(init,nlog.lik,control=list(maxit=maxit), method = method)
  if(method == "nlminb") OPT <- nlminb(init,nlog.lik,control=list(iter.max=maxit))
  OUT <- OPT$par
  return(OUT)
}
#--------------------------------------------------------------------
# PH model: to estimate the nonparametric part only 
#--------------------------------------------------------------------
nonparPH_mle <- function(init = NULL, y, delta, x2, x21.index = NULL, 
                         theta.hat, beta2.hat = NULL, dist, 
                          x0, h, ktype = "Epan", method = "Nelder-Mead", maxit = 10000, ...){
  if (is.na(h) || h <= 0L)
    stop("invalid 'bandwidth value.' ")
  if(is.null(colnames(x2)))
    stop("give a column names for the 'covariates'")
  # Required variables
  y <- as.vector(y)
  delta <- as.vector(delta) 
  n <- length(y)
  # profiling out the parametric coefficients 
  if(!is.null(x21.index) | !is.null(beta2.hat)){
    x21 <- as.matrix(x2[, x21.index])
    d21 <- ncol(x21)
    x22 <- as.matrix(x2[,-c(1:d21)])
    d22 <- ncol(x22)
    psi.x21 <- as.vector(x21%*%beta2.hat) 
  } else {psi.x21 <- rep(0, n)
          x22 <- as.matrix(x2)
          d22 <- dim(x22)[2]
    }
  # local points 
  newx0 = x22 - t(matrix(rep(x0,d22*n),d22,n))
  ker.h <- as.vector(apply(newx0/h, MARGIN = 2, FUN = ker, ktype = ktype)/h)
  newx <- cbind(1, newx0, newx0^2)# second order derivative 
  if(is.null(init)){
    fit.poly <- lm(y[delta==1]~poly(newx0[delta==1], degree = 2))
    init <- as.numeric(fit.poly$coefficients)
  }
  # TPA lognormal  PH model  
  if(dist == "tpalnorm"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      eta0 <- exp(theta.hat[1]); phi0 <- exp(theta.hat[2]);  alpha0 <- expit(theta.hat[3]) # par. for baseline distribution
      b2 <- par # vector of local regression coefficients 
      psi.x22 <- as.vector(newx%*%b2)
      psi.x2 <- psi.x21 + psi.x22
      
      lhaz0 <- htpa(y, eta0, phi0, alpha0, f0_N, F0_N, glog, log = TRUE) + psi.x2
      chaz0 <- Htpa(y, eta0, phi0, alpha0, F0_N, glog,  log = FALSE)*exp(psi.x2)
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0*ker.h, na.rm = T) + sum(chaz0*ker.h, na.rm = T)#  nloglik
      return(val) 
    }
  }
  # TPA loglogistic PH model  
  if(dist == "tpallogis"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      eta0 <- exp(theta.hat[1]); phi0 <- exp(theta.hat[2]);  alpha0 <- expit(theta.hat[3]) # par. for baseline distribution
      b2 <- par # vector of local regression coefficients  
      psi.x22 <- as.vector(newx%*%b2)
      psi.x2 <- psi.x21 + psi.x22
      
      lhaz0 <- htpa(y, eta0, phi0, alpha0, f0_Lo, F0_Lo, glog,  log = TRUE) + psi.x2
      chaz0 <- Htpa(y, eta0, phi0, alpha0, F0_Lo, glog, log = FALSE)*exp(psi.x2)
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0*ker.h, na.rm = T) + sum(chaz0*ker.h, na.rm = T)#  nloglik
      return(val) 
    }
  }
  
  # TPA logLaplace  model 
  if(dist == "tpalLap"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      eta0 <- exp(theta.hat[1]); phi0 <- exp(theta.hat[2]);  alpha0 <- expit(theta.hat[3]) # par. for baseline distribution
      b2 <- par # vector of local regression coefficients  
      psi.x22 <- as.vector(newx%*%b2)
      psi.x2 <- psi.x21 + psi.x22
      
      lhaz0 <- htpa(y, eta0, phi0, alpha0, f0_La, F0_La, glog,  log = TRUE) + psi.x2
      chaz0 <- Htpa(y, eta0, phi0, alpha0, F0_La, glog, log = FALSE)*exp(psi.x2)
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0*ker.h, na.rm = T) + sum(chaz0*ker.h, na.rm = T)#  nloglik
      return(val) 
    }
  }
  # EW - Model
  if(dist == "EW"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      par0 <- exp(theta.hat)
      b2 <- par # vector of local regression coefficients  
      psi.x22 <- as.vector(newx%*%b2)
      psi.x2 <- psi.x21 + psi.x22
      
      lhaz0 <- hexpweibull(y, par0[1], par0[2], par0[3], log = TRUE) + psi.x2
      chaz0 <- chexpweibull(y, par0[1], par0[2], par0[3])*exp(psi.x2)
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0*ker.h, na.rm = T) + sum(chaz0*ker.h, na.rm = T)#  nloglik
      return(val) 
    } 
  }
  if(method != "nlminb") OPT <- optim(init,nlog.lik,control=list(maxit=maxit), method = method)
  if(method == "nlminb") OPT <- nlminb(init,nlog.lik,control=list(iter.max=maxit))
  OUT <- OPT$par
  return(OUT)

}


#--------------------------------------------------------------------
# AFT model: to estimate the nonparametric part only 
#--------------------------------------------------------------------
nonparAH_mle <- function(init = NULL, y, delta, x, x11.index = NULL, 
                         theta.hat, beta1.hat = NULL, dist, 
                         x0, h, ktype = "Epan", method = "Nelder-Mead",
                         maxit = 10000, ...){
  if (is.na(h) || h <= 0L)
    stop("invalid 'bandwidth value.' ")
  if(is.null(colnames(x)))
    stop("give a column names for the 'covariates'")
  # Required variables
  y <- as.vector(y)
  delta <- as.vector(delta) 
  n <- length(y)
  # profiling out the parametric coefficients 
  if(!is.null(x11.index) | !is.null(beta1.hat)){
    x11 <- as.matrix(x[, x11.index])
    d11 <- ncol(x11)
    x12 <- as.matrix(x[,-c(1:d11)])
    d12 <- ncol(x12)
    psi.x11 <- as.vector(x11%*%beta1.hat) 
  } 
  else {psi.x11 <- rep(0, n)
         x12 <- as.matrix(x)
         d12 <- dim(x12)[2]
  }
  # local points 
  newx0 = x12 - t(matrix(rep(x0,d12*n),d12,n))
  ker.h <- as.vector(apply(newx0/h, MARGIN = 2, FUN = ker, ktype = ktype)/h)
  newx <- cbind(1, newx0, newx0^2)# second order derivative 
  
  if(is.null(init)){
    fit.poly <- lm(y[delta==1]~poly(newx0[delta==1], degree = 2))
    init <- as.numeric(fit.poly$coefficients)
  }
  # TPA lognormal  PH model  
  if(dist == "tpalnorm"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      eta0 <- exp(theta.hat[1]); phi0 <- exp(theta.hat[2]);  alpha0 <- expit(theta.hat[3]) # par. for baseline distribution
      b1 <- par # vector of local regression coefficients 
      psi.x12 <- as.vector(newx%*%b1)
      psi.x <- psi.x11 + psi.x12
      
      lhaz0 <- htpa(y*exp(psi.x), eta0, phi0, alpha0, f0_N, F0_N, glog, log = TRUE) 
      chaz0 <- Htpa(y*exp(psi.x), eta0, phi0, alpha0, F0_N, glog,  log = FALSE)*exp(-psi.x)
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0*ker.h, na.rm = T) + sum(chaz0*ker.h, na.rm = T)#  nloglik
      return(val) 
    }
  }
  # TPA loglogistic PH model  
  if(dist == "tpallogis"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      eta0 <- exp(theta.hat[1]); phi0 <- exp(theta.hat[2]);  alpha0 <- expit(theta.hat[3]) # par. for baseline distribution
      b1 <- par # vector of local regression coefficients 
      psi.x12 <- as.vector(newx%*%b1)
      psi.x <- psi.x11 + psi.x12
      
      lhaz0 <- htpa(y*exp(psi.x), eta0, phi0, alpha0, f0_Lo, F0_Lo, glog, log = TRUE) 
      chaz0 <- Htpa(y*exp(psi.x), eta0, phi0, alpha0, F0_Lo, glog,  log = FALSE)*exp(-psi.x)
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0*ker.h, na.rm = T) + sum(chaz0*ker.h, na.rm = T)#  nloglik
      return(val) 
    }
  }
  
  # TPA logLaplace  model 
  if(dist == "tpalLap"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      eta0 <- exp(theta.hat[1]); phi0 <- exp(theta.hat[2]);  alpha0 <- expit(theta.hat[3]) # par. for baseline distribution
      b1 <- par # vector of local regression coefficients 
      psi.x12 <- as.vector(newx%*%b1)
      psi.x <- psi.x11 + psi.x12
      
      lhaz0 <- htpa(y*exp(psi.x), eta0, phi0, alpha0, f0_La, F0_La, glog, log = TRUE) 
      chaz0 <- Htpa(y*exp(psi.x), eta0, phi0, alpha0, F0_La, glog,  log = FALSE)*exp(-psi.x)
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0*ker.h, na.rm = T) + sum(chaz0*ker.h, na.rm = T)#  nloglik
      return(val) 
    }
  }
  # EW - Model
  if(dist == "EW"){
    nlog.lik <- function(par){
      if(anyNA(par)) return(Inf)
      par0 <- exp(theta.hat)
      b1 <- par # vector of local regression coefficients 
      psi.x12 <- as.vector(newx%*%b1)
      psi.x <- psi.x11 + psi.x12
      
      lhaz0 <- hexpweibull(y*exp(psi.x), par0[1], par0[2], par0[3], log = TRUE) 
      chaz0 <- chexpweibull(y*exp(psi.x), par0[1], par0[2], par0[3])*exp(-psi.x)
      lhaz0[is.infinite(lhaz0)] <- NA
      chaz0[is.infinite(chaz0)] <- NA
      val <-  - sum(delta*lhaz0*ker.h, na.rm = T) + sum(chaz0*ker.h, na.rm = T)#  nloglik
      return(val) 
    } 
  }
  if(method != "nlminb") OPT <- optim(init,nlog.lik,control=list(maxit=maxit), method = method)
  if(method == "nlminb") OPT <- nlminb(init,nlog.lik,control=list(iter.max=maxit))
  OUT <- OPT$par
  return(OUT)
  
}
#===================================================================
# distance measures 
#===================================================================

# ASE for the survival  function 
surv_ase <- function(true_sp, pred_sp){
  if(length(true_sp)!=length(pred_sp))
    stop(" Error: surv_sp and pred_sp must have equal length")
  d <- (true_sp - pred_sp)^2
  return(ase = mean(d, na.rm=TRUE))
}
# absolute error 

surv_abs <- function(true_sp, pred_sp){
  if(length(true_sp)!=length(pred_sp))
    stop(" Error: surv_sp and pred_sp must have equal length")
  d <- abs(true_sp - pred_sp)
  return(ase = mean(d, na.rm=TRUE))
}
#====================================================================
# Predictions: survival, hazard, cumulative hazard and quantile
#====================================================================
#' theta: transformed baseline parameters 
#' tau: order of quantile (must be b/n 0 and 1)
#' dist: baseline distribution 
#' time: a single time 
#' in computation psihat.1 and psihat.2 are a  vector of zeros for AH models. 
#' whereas for AFT model psihat.1 = psihat.2 

predict_fun <- function(tau=NULL, theta,  time, psihat.1, psihat.2,
                        dist, type ="survival",  ...){
  
  time <- sort(time)
  
  psihat.dif <- psihat.2 - psihat.1
  
  if(dist=="EW") par <- exp(theta) else par <- c(exp(theta[1:2]), expit(theta[3]))
  
  if(type=="survival"){
    if(dist =="tpalnorm"){
      value <- ptpa(time*exp(psihat.1), par[1], par[2], par[3], 
                    F0_N, glog, lower.tail = FALSE)^exp(psihat.dif)
    }
    else if(dist =="tpallogis"){
      value <- ptpa(time*exp(psihat.1), par[1], par[2], par[3], 
                    F0_Lo, glog, lower.tail = FALSE)^exp(psihat.dif)
    }
    else if(dist =="EW"){
      value <- 1 - pexpweibull(time*exp(psihat.1), par[1], par[2], par[3])^exp(psihat.dif)
    } else stop("Baseline distribution is not provided")
  }
  else if(type=="hazard"){
    if(dist =="tpalnorm"){
      value <- htpa(time*exp(psihat.1), par[1], par[2], par[3], 
                    f0_N, F0_N, glog)*exp(psihat.2)
    }
    else if(dist =="tpallogis"){
      value <- htpa(time*exp(psihat.1), par[1], par[2], par[3], 
                    f0_Lo, F0_Lo, glog)*exp(psihat.2)
    }
    else if(dist =="EW"){
      value <- hexpweibull(time*exp(psihat.1), par[1], par[2], par[3])*exp(psihat.2)
    } else stop("Baseline distribution is not provided")
  }
  else if(type=="cmhazard"){
     if(dist =="tpalnorm"){
      value <- Htpa(time*exp(psihat.1), par[1], par[2], par[3], F0_N, glog)*exp(psihat.dif)
    }
   else  if(dist =="tpallogis"){
      value <- Htpa(time*exp(psihat.1), par[1], par[2], par[3], F0_Lo, glog)*exp(psihat.dif)
    }
    else if(dist =="EW"){
      value <- chexpweibull(time*exp(psihat.1), par[1], par[2], par[3])*exp(psihat.dif)
    } else stop("Baseline distribution is not provided")
  }
  else if(type=="quantile"){
     if(dist =="tpalnorm"){
      value <- exp(-psihat.1)*qtpa( 1- exp(log(1-tau)/exp(psihat.dif)), par[1], par[2], par[3], 
                                   F0_N, glog, glog.inv, QF_N)
    }
    else if(dist =="tpallogis"){
      value <- exp(-psihat.1)*qtpa(1- exp(log(1-tau)/exp(psihat.dif)), par[1], par[2], par[3], 
                                   F0_Lo, glog, glog.inv, QF_Lo)
    }
    else if(dist =="EW"){
      value <- exp(-psihat.1)*qexpweibull(1- exp(log(1-tau)/exp(psihat.dif)), par[1], par[2], par[3])
    } else stop("Baseline distribution is not provided")
  } 
  return(value)
}
#predict_fun <- Vectorize(predict_fun, vectorize.args = "time")
# ================================================================
# Predicted plot 
#==================================================================
predict_plot <- function(object, time,  xlab="Time", ylab = "Survival probability"){
  ord <- order(time)
  
  plot(time[ord], object[ord], xlab=xlab, ylab=ylab)
}
#====================================================================
# trapezoidal rule to estimate psi
#====================================================================
#' x0: vector of grid values 
#' dpsi.hat: estimate of the derivative of psi



# using  trapz function in pracma package 

trap_fun <- function(x0, dpsi.hat){
  if(is.unsorted(x0)){
    ord <- order(x0)
    x0 <- x0[ord]
    dpsi.hat <- dpsi.hat[ord]
  }

  m <- length(x0)
  psi.hat <- as.null(x0)
  psi.hat[1] <- (x0[1]/2)*dpsi.hat[1]
  for (k in 2:m) {
    s.x <- sapply(2:k, function(j) trapz(x= x0[(j-1):j], y=dpsi.hat[(j-1):j]))
    psi.hat[k] <- sum(s.x)
  }
  return(psi.hat = as.numeric(psi.hat))
}
# or using my own function 
trap_phihat <- function(x0, dpsi.hat){
  if(!is.sorted(x0))
    stop(" 'x0' should be sorted")
  m <- length(x0)
  psi.hat <- as.null(x0)
  psi.hat[1] <- (x0[1]/2)*dpsi.hat[1]
  for (k in 2:m) {
    s.x <- sapply(2:k, function(j) (x0[j] - x0[j-1])*(dpsi.hat[j] + dpsi.hat[j-1])/2)
    psi.hat[k] <- sum(s.x)
  }
  return(psi.hat = as.numeric(psi.hat))
}

#===============================================================
# functions to compute confidence intervals 
#===============================================================
# based on hessian() function from  "numDeriv" library

# nlogLik: - log-likelihood value
# mle: MLE estimates 

Conf.Int.hess <- function(nlogLik,mle,level=0.95){
  sd.int <- abs(qnorm(0.5*(1-level)))
  HESS <- hessian(func = nlogLik,x=mle)
  Fisher.Info <- solve(HESS)
  se <- sqrt(diag(Fisher.Info))
  U <- mle + sd.int*se
  L <- mle - sd.int*se
  C.I <- cbind(mle, se, L, U)
  names.row <- paste0("par", seq_along(1:length(mle)))
  rownames(C.I)<- names.row
  colnames(C.I)<- c("Estimates", "Std. Error","Lower","Upper")
  return(C.I)
}

#  based on optimHess() function with stats library
Conf.Int.optim <- function(nlogLik,mle,level=0.95){
  sd.int <- abs(qnorm(0.5*(1-level)))
  HESS <- optimHess(par=mle, fn=nlogLik, gr=NULL)
  Fisher.Info <- solve(HESS)
  se <- sqrt(diag(Fisher.Info))
  U <- mle + sd.int*se
  L <- mle - sd.int*se
  C.I <- cbind(mle, se, L, U)
  names.row <- paste0("par", seq_along(1:length(mle)))
  rownames(C.I)<- names.row
  colnames(C.I)<- c("Estimates", "Std. Error","Lower","Upper")
  return(C.I)
}

#  based on a finite difference numerical approximation for 
# the hessian matrix using fdHess()  function from nlme library
Conf.Int.FD <- function(nlogLik,mle,level=0.95,...){
      sd.int <- abs(qnorm(0.5*(1-level)))
      HESS <- fdHess(pars=mle, fun = nlogLik)
      Fisher.Info <- solve(HESS$Hessian)
      se <- sqrt(diag(Fisher.Info))
      U <- mle + sd.int*se
      L <- mle - sd.int*se
      C.I <- cbind(mle, se, L, U)
      names.row <- paste0("par", seq_along(1:length(mle)))
      rownames(C.I)<- names.row
      colnames(C.I)<- c("Estimates", "Std. Error","Lower","Upper")
      return(C.I)
}

 # CI.hessian <- Conf.Int.hess(nlogLik= norm.parPH$nloglik, mle=norm.parPH$OPT$par,level=0.95)
 # CI.optim <- Conf.Int.optim(nlogLik= norm.parPH$nloglik, mle=norm.parPH$OPT$par,level=0.95)
# library(knitr)
# print(kable(CI.hessian,digits=4))

#---------------------------------------------------------------
# Confidence intervals for the survival based on Bootstrap
#---------------------------------------------------------------

#=====================================================
# Generate initial value for EW [Section 3.1.1]
# This function is used in fit.model() [see below]
# theta0: initial values for log(tau) and beta0;
# c(1,1) is used in fit.model()
# st = lifetime
# status = censoring indicator
# xx = design matrix of the form (1,x_1,...,x_p)
# reference: Exponentiated Weibull Regression for Time-to-Event Data
#=====================================================
init.EW <- function(theta0, st, status, xx) {
      pllik <- function(theta0, st, status) {
            y <- log(st)
            r <- sum(status)
            tau <- exp(theta0[1])
            beta0 <- theta0[2]
            w <- (y - beta0) / tau
            a <- -expm1(-exp(w))
            gam <- -r / sum(log(a[status == 1]))
            lf0 <- status * ((gam - 1) * log(a) + w - exp(w)) + (1 - status) * log1p(-a^gam)
            lf <- r * log(gam) - r * log(tau) + sum(lf0)
            return(-as.numeric(lf))
      }
      method <- "nlminb"
      options(warn = -1)
      fit <- nlminb(start = theta0, objective = pllik, st = st,
                    status = status)
      conv <- 0
      if (fit$message != "relative convergence (4)") {
            fit <- optim(par = theta0, fn = pllik, st = st,
                         status = status)# method = "BFGS", 
            method = "optim-BFGS"
            conv <- fit$convergence
      }
      options(warn = 0)
      
      init.logtau <- fit$par[1]
      init.beta0 <- fit$par[2]
      w <- (log(st) - init.beta0) / exp(init.logtau)
      a <- -expm1(-exp(w))
      init.loggam <- as.numeric(log(-sum(status) / sum(log(a[status == 1]))))
      if (ncol(xx) > 1) {
            xx0 <- data.frame(xx[, -1])
            fit.w <- survreg(Surv(st, status)~., data = xx0, dist = "weibull")
            init.beta <- fit.w$coef[-1]
            init <- c(init.logtau, init.loggam, init.beta0, init.beta)
            names(init) <- c("logtau", "loggam", "intercept", names(xx)[-1])
      } else{
            init <- c(init.logtau, init.loggam, init.beta0)
            names(init) <- c("logtau", "loggam", "intercept")
      }
      return(list(init = init, convergence = conv, method = method, message = fit$message))
}

