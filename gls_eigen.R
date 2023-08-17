library(pracma)
library(EQL)

if(!dir.exists('gls_eigen')) dir.create('gls_eigen')

id = as.numeric(commandArgs(trailingOnly=T))
set.seed(id)
out.file = file.path('gls_eigen', paste0('run_', id, '.rds'))

# Squared Exponential Kernel
se_kernel <- function(s, l) {
  dist_matrix <- as.matrix(dist(s))
  kernel_matrix <-  exp(-dist_matrix^2 / (2 * l^2))
  return(kernel_matrix)
}


# Parameters
num_points <- 3000
l <- 1
mean <- 0
sd <- 1
tau = 1
tau.ratio = 2
kappa = .25
h.coef = g.coef = rep(0, num_points)
k.max = 5
h.coef[1:k.max] = rnorm(k.max)
g.coef[1:k.max] = rnorm(k.max)
beta = 1
n.trials = 500

# Hermite Eigenfunction
hermite_eigenfunction <- function(x, k, sd, l) {
  a_inv <- 4 * sd^2
  b_inv <- 2 * l^2
  c <- sqrt(a_inv^(-2) + 2*a_inv^(-1)*b_inv^(-1))
  A <- a_inv^(-1) + b_inv^(-1) + c
  B <- b_inv^(-1) / A
  return(exp(- (c - a_inv^(-1)) * x^2) * hermite(sqrt(2*c)*x, k, prob = F) )
}

bias.exact = function(X, g, W){
  num = t(X)%*%W%*%g/num_points
  denom = t(X)%*%W%*%X/num_points
  num/denom
}

bias.pred = function(h.coef, g.coef, sd, l, kappa){
  denom = t(X) %*% W %*% X/num_points
  a_inv <- 4 * sd^2
  b_inv <- 2 * l^2
  c <- sqrt(a_inv^(-2) + 2*a_inv^(-1)*b_inv^(-1))
  A <- a_inv^(-1) + b_inv^(-1) + c
  B <- b_inv^(-1) / A
  
  k.max = max( max(which(h.coef!=0), max(which(g.coef!=0))) )
  
  e.vals = sapply(0:(k.max-1), function(k){
    sqrt(2/A/a_inv)*B^k
  })
  
  num = sum(g.coef[1:k.max]*h.coef[1:k.max]/(num_points*e.vals+tau^2))
  return(num/denom)
}

be = vector('numeric', length=n.trials)
bp = vector('numeric', length=n.trials)

  # Generate points according to Gaussian density
  s <- rnorm(num_points, mean, sd)
  h = Reduce("+", lapply(0:(k.max-1), function(k){
    v = hermite_eigenfunction(s, k, sd=sd, l=l)
    h.coef[k+1]*v/(sqrt(sum(v^2)))*sqrt(num_points)
  }) )
  g = Reduce("+", lapply(0:(k.max-1), function(k){
    v = hermite_eigenfunction(s, k, sd=sd, l=l)
    g.coef[k+1]*v/(sqrt(sum(v^2)))*sqrt(num_points)
  }) )
  
  eta = rnorm(num_points, sd=kappa)
  X =  h + eta
  eps = rnorm(num_points, sd=tau)
  Y = beta*X + g+ eps
  
  # Calculate the kernel matrix
  kernel_matrix <- se_kernel(s, l)
  cov_matrix <- kernel_matrix + diag(tau^2, num_points)*tau.ratio
  W <- chol2inv(chol(cov_matrix))
  
  be = bias.exact(X,g,W)
  bp = bias.pred(h.coef,g.coef,sd,l,kappa)


saveRDS(list(be=be,bp=bp), file=out.file)

quit()
gls_estimator <- function(X, Y, W) {
  num = X%*%W%*%Y
  denom = X%*%W%*%X
  return(num/denom)
}

ests = sapply(1:200, function(i){
  Y = beta*X + g+ rnorm(num_points, sd=tau)
  gls_estimator(X,Y,W)
})
