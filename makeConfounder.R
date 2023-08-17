library(BRISC) 
n = 10000 
m = 1000 
set.seed(124) 
seeds = sample.int(1e6, replace=F, size = m) 
coords = 10*cbind(runif(n), runif(n)) 
g = BRISC_simulation(coords, phi = .25, 
cov.model='spherical', sim_number = m,
                     seeds = seeds)
saveRDS(g, 'confounder.rds')
