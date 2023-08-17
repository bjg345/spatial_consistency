library(ggplot2)
library(dplyr)
library(xtable)
library(gridExtra)
library(magrittr)
library(cowplot)
library(tidyr)
z = qnorm(.975)

setwd('gls_eigen')
data_list <- lapply(1:500, function(i) {
  file_name <- paste0("run_", i, ".rds")
  data <- readRDS(file_name)
  data_frame <- data.frame(bp = data$bp, be = data$be)
  return(data_frame)
})

data_all <- do.call(rbind, data_list)

# Create ggplot with 45-degree line
p <- ggplot(data_all, aes(x = bp, y = be)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Eigenfunction analysis of GLS bias",
       y = "Predicted Bias",
       x = "Exact Bias") 
ggsave('../gls_eigen.png', p, width = 7, height = 5, units = 'in')

setwd('..')
conf = readRDS('confounder.rds')
values = conf$output.data[,1]
coords = conf$coords
p = qplot(x=coords[,1], y= coords[,2]) + geom_point(aes(color=values)) + xlab('X-coordinate') + ylab('Y-coordinate') +
  ggtitle('Values of confounding variable')+theme(text = element_text(size = 16))
ggsave('confounder.png', p, width = 7, height = 5, units = 'in')
 

setwd('fixed')
f = list.files()
n = length(f)
res = lapply(f, readRDS)
ols.ests = sapply(res, function(x) x$ols$est)
gls.ests = sapply(res, function(x) x$gls$est)
gam.ests = sapply(res, function(x) x$gam$est)
gam.fx.ests = sapply(res, function(x) x$gam.fx$est)
plus.ests = sapply(res, function(x) x$plus$est)
# Create dataframe
df = data.frame('OLS'=ols.ests-1,'GLS'=gls.ests-1, 'GAM'=gam.ests-1, 'GAM.fx'=gam.fx.ests-1, 'Spatial.plus'=plus.ests-1) %>%
  pivot_longer(cols = OLS:Spatial.plus, names_to = 'Method', values_to= 'Error')

# Create the first plot and modify the order of the x-axis
p1 = ggplot(df, aes(x=Method, y=Error, fill=Method)) +
  geom_boxplot() + theme(legend.position='none', text = element_text(size = 8)) +
  scale_x_discrete(limits = c("OLS", "GLS", "GAM", "GAM.fx", "Spatial.plus"))

# Create the second plot without the OLS method
p2 = ggplot(df %>% filter(Method != 'OLS'), aes(x=Method, y=Error, fill=Method))+
  geom_boxplot() + theme(legend.position='none', text = element_text(size = 8))+
  scale_x_discrete(limits = c("GLS", "GAM", "GAM.fx", "Spatial.plus"))

# Create a title spanning both plots
title <- ggdraw() + draw_label("Estimator error distributions\nfor the fixed confounder", fontface='bold', size=12)

# Combine the two plots side-by-side with the title
combined_plots = plot_grid(title, plot_grid(p1, p2, ncol=2), nrow=2, rel_heights=c(0.1, 1))

# Save the combined plots as a PNG file
save_plot('../fixed.png', combined_plots, bg='white')
 

cov.mod = function(method){
  if(method == 'gls'){
    out = sapply(1:n, function(i){
      ci = res[[i]][[method]]$ci.mod
      between(1, ci[1], ci[2])
    }) 
    return(mean(out))
  } else
  out = sapply(1:n, function(i) {
    est = res[[i]][[method]]$est
    se = res[[i]][[method]]$se.mod
    if(is.null(se)) return(NA)
    between(1, est-z*se, est+z*se)
    })
    mean(out)
}
cov.sub = function(method){
  out = sapply(1:n, function(i) {
    est = res[[i]][[method]]$est
    se = res[[i]][[method]]$se.sub
    if(is.null(se)) return(NA)
    between(1, est-z*se, est+z*se)
    })
    mean(out)
}
methods = c('ols', 'gls', 'gam', 'gam.fx', 'plus')
cov.mods = sapply(methods, cov.mod)
print(cov.mods)
cov.subs = sapply(methods, cov.sub)

setwd('../random')
f = list.files()
n = length(f)
res = lapply(f, readRDS)
ols.ests = sapply(res, function(x) x$ols$est)
gls.ests = sapply(res, function(x) x$gls$est)
gam.ests = sapply(res, function(x) x$gam$est)
gam.fx.ests = sapply(res, function(x) x$gam.fx$est)
plus.ests = sapply(res, function(x) x$plus$est)
# Create dataframe
df = data.frame('OLS'=ols.ests-1,'GLS'=gls.ests-1, 'GAM'=gam.ests-1, 'GAM.fx'=gam.fx.ests-1, 'Spatial.plus'=plus.ests-1) %>%
  pivot_longer(cols = OLS:Spatial.plus, names_to = 'Method', values_to= 'Error')

# Create the first plot and modify the order of the x-axis
p1 = ggplot(df, aes(x=Method, y=Error, fill=Method)) +
  geom_boxplot() + theme(legend.position='none', text = element_text(size = 8)) +
  scale_x_discrete(limits = c("OLS", "GLS", "GAM", "GAM.fx", "Spatial.plus"))

# Create the second plot without the OLS method
p2 = ggplot(df %>% filter(Method != 'OLS'), aes(x=Method, y=Error, fill=Method))+
  geom_boxplot() + theme(legend.position='none', text = element_text(size = 8))+
  scale_x_discrete(limits = c("GLS", "GAM", "GAM.fx", "Spatial.plus"))

# Create a title spanning both plots
title <- ggdraw() + draw_label("Estimator error distributions\nfor the random confounder", fontface='bold', size=12)

# Combine the two plots side-by-side with the title
combined_plots = plot_grid(title, plot_grid(p1, p2, ncol=2), nrow=2, rel_heights=c(0.1, 1))

# Save the combined plots as a PNG file
save_plot('../random.png', combined_plots, bg='white')
 

methods = c('ols', 'gls', 'gam', 'gam.fx', 'plus')
cov.mods = sapply(methods, cov.mod)
print(cov.mods)
cov.subs = sapply(methods, cov.sub)



setwd('../ols')
f = list.files(pattern='xcor')
n = length(f)
res = lapply(f, readRDS)
ols.ests = sapply(res, function(x) x$ols.est)
gls.ests = sapply(res, function(x) x$gls.est)
ols.cov = sapply(res, function(x) x$ols.cov)
ols.cov.cor=mean(ols.cov)
gls.cov = sapply(res, function(x) x$gls.cov)
gls.cov.cor=mean(gls.cov)

df.ols = data.frame(ests=ols.ests,cov=ols.cov, type='OLS')
df.gls = data.frame(ests=gls.ests,cov=gls.cov, type='GLS')
df = rbind(df.ols, df.gls)

p.cor=ggplot(df, aes(x=ests, fill=type))+geom_density(alpha=.5)+xlab('estimate')+
  ggtitle('Estimates for the correlated-X process') +
  geom_vline(xintercept=1)+theme(text = element_text(size = 16))

ggsave('../xcor.png', p.cor)
 

f = list.files(pattern='xuncor')
n = length(f)
res = lapply(f, readRDS)
ols.ests = sapply(res, function(x) x$ols.est)
gls.ests = sapply(res, function(x) x$gls.est)
ols.cov = sapply(res, function(x) x$ols.cov)
ols.cov.uncor=mean(ols.cov)
gls.cov = sapply(res, function(x) x$gls.cov)
gls.cov.uncor=mean(gls.cov)


df.ols = data.frame(ests=ols.ests,cov=ols.cov, type='OLS')
df.gls = data.frame(ests=gls.ests,cov=gls.cov, type='GLS')
df = rbind(df.ols, df.gls)

p.uncor=ggplot(df, aes(x=ests, fill=type))+geom_density(alpha=.5)+xlab('estimate')+
  ggtitle('Estimates for the uncorrelated-X process')  +
  geom_vline(xintercept=1)+theme(text = element_text(size = 16))


ggsave('../xuncor.png', p.uncor)
 

tab = data.frame(Simulation=c('Correlated X', 'Correlated X', 
                              'Uncorrelated X', 'Uncorrelated X'),
                 Estimator=c('OLS', 'GLS','OLS','GLS'),
                 Coverage=c(ols.cov.cor, gls.cov.cor, ols.cov.uncor, gls.cov.uncor))
print(xtable(tab, digits=3))

setwd('../avsg')

f = list.files()
n = length(f)
res = lapply(f, readRDS)
ols.ests = sapply(res, function(x) coef(x$ols.mod)[2])
geo.ests = sapply(res, function(x) x$geo.mod$est)
area.ests = sapply(res, function(x) x$area.mod$est)
df = data.frame('OLS'=ols.ests-1, 'Grouped'=area.ests-1, 'Spatial'=geo.ests-1) %>%
  pivot_longer(cols = OLS:Spatial, names_to = 'Method', values_to= 'Error')
p = ggplot(df, aes(x=Method, y=Error, fill=Method)) +
  geom_boxplot() + theme(legend.position='none') +
  ggtitle('Estimator error distributions\nfor the linear confounder') +theme(text = element_text(size = 16))

ggsave('../avsg.png', p)
 

cov = function(method){
  out = sapply(1:n, function(i){
      ci = res[[i]][[method]]$ci
      between(1, ci[1], ci[2])
    }) 
    return(mean(out))
}
methods = c('area.mod', 'geo.mod')
cov = sapply(methods, cov)
print(cov)

quit()
setwd('../area')
f = list.files()
n = length(f)
res = lapply(f, readRDS)
ols.ests = sapply(res, function(x) coef(x$ols.mod)[2])
gls.ests = sapply(res, function(x) x$gls.mod$est)
gls.se = sapply(res, function(x) x$gls.mod$se)
gls.cov = sapply(1:n, function(i){
  between(1, gls.ests[i]-z*gls.se[i], gls.ests[i]+z*gls.se[i])
}
)
gls.cov.area=mean(gls.cov)

df.ols = data.frame(ests=ols.ests, type='OLS')
df.gls = data.frame(ests=gls.ests, type='GLS')
df = rbind(df.ols, df.gls)

p.area=ggplot(df, aes(x=ests, fill=type))+geom_density(alpha=.5)+xlab('estimate')+
  ggtitle('Estimates for the\ngrouped process')  +
  geom_vline(xintercept=1)+theme(text = element_text(size = 16))

ggsave('../area.png', p.area)
 

setwd('../geo')
f = list.files()
n = length(f)
res = lapply(f, readRDS)
ols.ests = sapply(res, function(x) coef(x$ols.mod)[2])
gls.ests = sapply(res, function(x) x$gls.mod$est)
gls.se = sapply(res, function(x) x$gls.mod$se)
gls.cov = sapply(1:n, function(i){
  between(1, gls.ests[i]-z*gls.se[i], gls.ests[i]+z*gls.se[i])
}
)
gls.cov.geo=mean(gls.cov)

df.ols = data.frame(ests=ols.ests, type='OLS')
df.gls = data.frame(ests=gls.ests, type='GLS')
df = rbind(df.ols, df.gls)

p.geo=ggplot(df, aes(x=ests, fill=type))+geom_density(alpha=.5)+xlab('estimate')+
  ggtitle('Estimates for the\ncontinuous process')  +
  geom_vline(xintercept=1)
ggsave('../geo.png', p.geo)+theme(text = element_text(size = 16))
 

tab = data.frame(Simulation=c('Grouped', 'Continuous'),
                 GLS.Coverage=c(gls.cov.area, gls.cov.geo))
print(xtable(tab, digits=3))
