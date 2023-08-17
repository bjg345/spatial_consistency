library(mgcv)
library(geoR)
library(mvtnorm)
library(BRISC)

source('functions.R')

out.dir = 'fixed'
if(!dir.exists(out.dir)) dir.create(out.dir)

id = as.numeric(commandArgs(trailingOnly=T))
set.seed(id+1000)

r = 20

confounder = readRDS('confounder.rds')

coords = confounder$coords
n = nrow(coords)
g = gam(confounder$output.data[,1]~s(coords[,1], coords[,2], k=200))$fitted.values
x = g + rnorm(n, sd=1)

y = x + g + rnorm(n,sd=1/2)

ols.out = mod.ols(x,y,coords)

gam.base = mod.gam(x,y,coords)
gam.sub = sapply(1:120, function(i){
  set.seed(i)
  sub = sample.int(n, replace=F, size = n/r)
  mod.gam(x[sub], y[sub], coords[sub,])$est
})
gam.out = list(est=gam.base$est, se.mod=gam.base$se.mod, se.sub = sd(gam.sub)/sqrt(r))

gam.fx.base = mod.gam.fx(x,y,coords)
gam.fx.sub = sapply(1:120, function(i){
  set.seed(i)
  sub = sample.int(n, replace=F, size = n/r)
  mod.gam.fx(x[sub], y[sub], coords[sub,])$est
})
gam.fx.out = list(est=gam.fx.base$est, se.mod=gam.fx.base$se.mod, se.sub = sd(gam.fx.sub)/sqrt(r))

gls.base = mod.gls(x,y,coords)
gls.sub = sapply(1:120, function(i){
  set.seed(i)
  sub = sample.int(n, replace=F, size = n/r)
  mod.gls(x[sub], y[sub], coords[sub,])$est
})
gls.boot = BRISC_bootstrap(gls.base$mod, n_boot = 120)
ci = gls.boot$confidence.interval[,2]
gls.out = list(est=gls.base$est, se.sub = sd(gls.sub)/sqrt(r), ci.mod = ci)

plus.base = mod.plus(x,y,coords)
plus.sub = sapply(1:120, function(i){
  sub = sample.int(n, replace=F, size = n/r)
  mod.plus(x[sub], y[sub], coords[sub,])$est
})
plus.out = list(est=plus.base$est, se.mod=plus.base$se.mod, se.sub = sd(plus.sub)/sqrt(r))

saveRDS(list(ols=ols.out, gam=gam.out, gam.fx=gam.fx.out, gls=gls.out, plus=plus.out), file.path(out.dir, paste0('res_', id, '.rds')))

