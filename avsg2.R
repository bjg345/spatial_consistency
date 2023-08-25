library(lme4) 
library(magrittr)
library(geoR)


if(!dir.exists('area')) dir.create('area')
if(!dir.exists('geo')) dir.create('geo')

id = as.numeric(commandArgs(trailingOnly=T))
set.seed(id)

#if(file.exists( file.path('geo', paste0('res_', id, '.rds')))) quit()

#areal
N <- 100 
k <- 10 
grp <- rep(1:N, each = k) 
Z <- rep(rnorm(N, 1:N)/10, each = k) 
X <- rnorm(N*k, Z) 
Y <- X + Z + rnorm(N*k) 
gls.mod <- lmer(Y ~ X + (1| factor(grp))) 
ols.mod <- lm(Y~X)
est=gls.mod@beta[2]
r = 20 
ests = sapply(1:120, function(i){
    grps = sample.int(N, N/r)
    ind = which(grp %in% grps)
    Xsub=X[ind]
    Ysub=Y[ind]
    grpsub=grp[ind]
    df = data.frame(Xsub=Xsub, Ysub=Ysub, grpsub=factor(grpsub))
    mod = lmer(Ysub ~ Xsub + (1| grpsub)) 
    mod@beta[2]
  }
) 
se = sd(ests)/sqrt(r)
saveRDS(list(ols.mod=ols.mod, gls.mod=list(est = est, se = se)), file.path('area', paste0('res_', id, '.rds')))
quit()
#geostatistical
N <- 1000 
Z <- runif(N, 0, 10) 
X <- rnorm(N, Z) 
Y <- X + Z + rnorm(N) 
s = as.matrix(dist(Z)) 
ols.mod = lm(Y~X)
gls.mod = likfit(coords=cbind(Z,0),data=Y,ini.cov.pars=c(10,10),fix.kappa=T,
              kappa=.5,cov.model='matern',trend=trend.spatial(~X),messages=F,lik.method='REML') 
est=gls.mod$beta[2] 
r = 20 
ests = sapply(1:120, function(i){
    ind = sample.int(N, N/r)
    Xsub=X[ind]
    Zsub=Z[ind]
    Ysub=Y[ind]
    mod = likfit(coords=cbind(Zsub,0),data=Ysub,ini.cov.pars=c(10,10),fix.kappa=T,kappa=.5,
                 cov.model='matern',trend=trend.spatial(~Xsub),messages=F,lik.method='REML')
    mod$beta[2]
  }
) 
se = sd(ests)/sqrt(r)

saveRDS(list(ols.mod=ols.mod, gls.mod=list(est = est, se = se)), file.path('geo', paste0('res_', id, '.rds')))
