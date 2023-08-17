if(!dir.exists('ols')) dir.create('ols')
setwd('ols')

id = as.numeric(commandArgs(trailingOnly=T))
set.seed(id)

library(magrittr)
library(geoR)
library(mvtnorm)
library(dplyr)
z = qnorm(.975)
n = 1000



for(type in c('xcor', 'xuncor')){
  S = runif(n, 0, 10)
  Sig = exp(-as.matrix(dist(S))) + diag(n)
  X = rnorm(n)
  if(type == 'xcor') X=X+ as.vector(rmvnorm(1, sigma = Sig ))
  
  eps = as.vector(rmvnorm(1, sigma = Sig))
  
  Y = X + eps
  
  ols.mod = lm(Y~X)
  ols.est = coef(ols.mod)[2]
  ols.ci = confint(ols.mod)[2,]
  ols.cov = between(1, ols.ci[1], ols.ci[2])
  
  gls.mod = likfit(coords=cbind(S, 0),data=Y,ini.cov.pars=c(1,1),fix.kappa=T,kappa=.5,cov.model='matern',trend=trend.spatial(~X),messages=F,lik.method='REML')
  gls.est = gls.mod$beta[2]
  gls.se = sqrt(gls.mod$beta.var[2,2])
  gls.ci = c(gls.est - z*gls.se, gls.est + z*gls.se) 
  gls.cov = between(1, gls.ci[1], gls.ci[2])
  
  saveRDS(list(ols.est=ols.est, ols.ci=ols.ci, ols.cov=ols.cov, gls.est=gls.est, gls.ci=gls.ci, gls.cov=gls.cov), paste0(type, '_', id, '.rds'))
}

  
