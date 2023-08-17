library(lme4) 
library(magrittr)
library(BRISC)
library(lme4)



if(!dir.exists('avsg')) dir.create('avsg')

id = as.numeric(commandArgs(trailingOnly=T))
set.seed(id)

z <- qnorm(.975)
N <- 300 
k <- 10  
Z <- rep(1:N/10, each = k) 
X <- rnorm(N*k, Z) 
Y <- X + Z + rnorm(N*k) 
r <- 20

ols.mod = lm(Y~X)
gls.mod.geo = BRISC_estimation(coords=cbind(Z, 0), y=Y, x=cbind(1,X))
gls.mod.area = lmer(Y ~ X + (1| factor(Z))) 

geo.est = gls.mod.geo$Beta[2]

gls.boot = BRISC_bootstrap(gls.mod.geo, n_boot = 120)
geo.ci = gls.boot$confidence.interval[,2]

area.est = coef(summary(gls.mod.area))[2,1]

area.ests = sapply(1:120, function(i){
    grps = sample.int(N, N/r)
    ind = which(as.integer(10*Z) %in% grps)
    Xsub=X[ind]
    Ysub=Y[ind]
    Zsub=Z[ind]
    df = data.frame(Xsub=Xsub, Ysub=Ysub, Zsub=factor(Zsub))
    mod = lmer(Ysub ~ Xsub + (1| Zsub)) 
    mod@beta[2]
  }
) 
area.se = sd(area.ests)/sqrt(r)
area.ci = c(area.est-z*area.se, area.est+z*area.se)

saveRDS(list(ols.mod=ols.mod, geo.mod=list(est = geo.est, ci = geo.ci), area.mod=list(est = area.est, ci = area.ci)),
 file.path('avsg', paste0('res_', id, '.rds')))
