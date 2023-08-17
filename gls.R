
library(geoR)

id = as.numeric(commandArgs(trailingOnly=T))

n = 1000
r = .05

s1 = runif(n, -1, 1)
s2 = runif(n, -1, 1)

U = s1^3+s2+sin(3*s2)

X = U + rnorm(n)

Y = X + U + rnorm(n)

sub = function(i, n, r, Y, XX, s1, s2){
  library(geoR)
  ind = sample.int(n, n*r)
  Y = Y[ind]
  X = XX[ind]
  s1 = s1[ind]
  s2 = s2[ind]
  mod = likfit(coords=cbind(s1,s2),data=Y,ini.cov.pars=c(1,1),fix.kappa=T,kappa=.5,cov.model='matern',trend=trend.spatial(~X),messages=F,lik.method='ML')
  mod$beta[2]
}

sub = function(i, n, r, Y, XX, s1, s2){
  library(geoR)
  ind = sample.int(n, n*r)
  Y = Y[ind]
  X = XX[ind]
  s1 = s1[ind]
  s2 = s2[ind]
  mod = likfit(coords=cbind(s1,s2),data=Y,ini.cov.pars=c(1,1),fix.kappa=T,kappa=.5,cov.model='matern',trend=trend.spatial(~X),messages=F,lik.method='ML')
  mod$beta[2]
}

jack = function(i, n, r, Y, XX, s1, s2){
  ind = setdiff(1:n, i)
  Y = Y[ind]
  X = XX[ind]
  s1 = s1[ind]
  s2 = s2[ind]
  mod = likfit(coords=cbind(s1,s2),data=Y,ini.cov.pars=c(1,1),fix.kappa=T,kappa=.5,cov.model='matern',trend=trend.spatial(~X),messages=F,lik.method='ML')
  mod$beta[2]
}

sub.samp = sapply( X = 1:120, FUN = sub, n =n, r = r, Y=Y, XX=X, s1=s1, s2=s2)
se = sd(sub.samp)*sqrt(r)

#jack.samp = sapply( X = 1:n, FUN = jack, n =n, r = r, Y=Y, XX=X, s1=s1, s2=s2)
#se = sd(jack.samp)

mod = likfit(coords=cbind(s1,s2),data=Y,ini.cov.pars=c(1,1),fix.kappa=T,kappa=.5,cov.model='matern',trend=trend.spatial(~X),messages=F,lik.method='ML')
est =  mod$beta[2]

res = data.frame(est = est, se = se)

saveRDS(res, paste0('gls_', id, '.rds'))

