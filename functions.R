mod.ols = function(x, y, coords){
  mod = lm(y~x)
  est = coef(mod)[2]
  se.mod = sqrt(vcov(mod)[2,2])
  list(est = est, se.mod = se.mod)
}


mod.gam = function(x, y, coords){
  mod = gam(y~x+s(coords[,1], coords[,2], k=200))
  est = coef(mod)[2]
  se.mod = sqrt(vcov(mod)[2,2])
  list(est = est, se.mod = se.mod)
}

mod.gam.fx = function(x, y, coords){
  mod = gam(y~x+s(coords[,1], coords[,2], fx=T, k=200))
  est = coef(mod)[2]
  se.mod = sqrt(vcov(mod)[2,2])
  list(est = est, se.mod = se.mod)
}

mod.gls = function(x, y, coords){
  mod = BRISC_estimation(coords, y, cbind(1,x))
  est = mod$Beta[2]
  list(est = est, mod = mod)
}

mod.plus = function(x, y, coords){
  mod.init = gam(x~s(coords[,1], coords[,2], k=200))
  mod = gam(y~mod.init$resid + s(coords[,1], coords[,2], k=200))
  est = coef(mod)[2]
  se.mod = sqrt(vcov(mod)[2,2])
  list(est = est, se.mod = se.mod)
}


