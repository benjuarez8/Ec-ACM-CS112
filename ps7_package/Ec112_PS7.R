## ----------------------------------------------------------------------------------
rm(list=ls())
set.seed(123)
library(Rcpp)
library(rstan)
library(bayesplot)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

data = read.csv("~/Desktop/MetabolicRate.csv")
Mrate = data$Mrate
BodySize = data$BodySize
logMrate = log(Mrate)
logBodySize = log(BodySize)
Instar = data$Instar
N = length(data[,1])


## ----------------------------------------------------------------------------------
num_iter = 10000
fit = stan("Ec122_PS7.stan", iter = num_iter, chains = 4, 
           data = list(log_Mrate=logMrate,
                       log_BodySize=logBodySize,
                       Instar=Instar,
                       N=N))
print(fit,probs = c(0.25,0.5,0.75))


## ----------------------------------------------------------------------------------
traceplot(fit, pars = c("beta0", "beta1", "beta2", "sigma"), 
          inc_warmup = FALSE, nrow = 4, window = c(num_iter-2000, num_iter))


## ----------------------------------------------------------------------------------
pairs(fit, pars = c("beta0", "beta1", "beta2", "sigma", "lp__"))


## ----------------------------------------------------------------------------------
yrep = as.matrix(fit, pars = "yPred")
indices = sample(nrow(yrep),100)
y = logMrate
ppc_dens_overlay(logMrate, yrep[indices,])


## ----------------------------------------------------------------------------------
library(gridExtra)
new_data = subset(data, Instar==1)
fit = stan("Ec122_PS7.stan", iter = 10000, chains = 4,
         data = list(log_Mrate=log(new_data$Mrate),
                     log_BodySize=log(new_data$BodySize),
                     Instar=new_data$Instar,
                     N=length(new_data[,1])))
yrep = as.matrix(fit, pars = "yPred")
indices = sample(nrow(yrep),100)
p1 = ppc_dens_overlay(log(new_data$Mrate), yrep[indices,])

new_data = subset(data, Instar==2)
fit = stan("Ec122_PS7.stan", iter = 10000, chains = 4,
         data = list(log_Mrate=log(new_data$Mrate),
                     log_BodySize=log(new_data$BodySize),
                     Instar=new_data$Instar,
                     N=length(new_data[,1])))
yrep = as.matrix(fit, pars = "yPred")
indices = sample(nrow(yrep),100)
p1 = ppc_dens_overlay(log(new_data$Mrate), yrep[indices,])

new_data = subset(data, Instar==2)
fit = stan("Ec122_PS7.stan", iter = 10000, chains = 4,
         data = list(log_Mrate=log(new_data$Mrate),
                     log_BodySize=log(new_data$BodySize),
                     Instar=new_data$Instar,
                     N=length(new_data[,1])))
yrep = as.matrix(fit, pars = "yPred")
indices = sample(nrow(yrep),100)
p2 = ppc_dens_overlay(log(new_data$Mrate), yrep[indices,])

new_data = subset(data, Instar==3)
fit = stan("Ec122_PS7.stan", iter = 10000, chains = 4,
         data = list(log_Mrate=log(new_data$Mrate),
                     log_BodySize=log(new_data$BodySize),
                     Instar=new_data$Instar,
                     N=length(new_data[,1])))
yrep = as.matrix(fit, pars = "yPred")
indices = sample(nrow(yrep),100)
p3 = ppc_dens_overlay(log(new_data$Mrate), yrep[indices,])

new_data = subset(data, Instar==4)
fit = stan("Ec122_PS7.stan", iter = 10000, chains = 4,
         data = list(log_Mrate=log(new_data$Mrate),
                     log_BodySize=log(new_data$BodySize),
                     Instar=new_data$Instar,
                     N=length(new_data[,1])))
yrep = as.matrix(fit, pars = "yPred")
indices = sample(nrow(yrep),100)
p4 = ppc_dens_overlay(log(new_data$Mrate), yrep[indices,])

new_data = subset(data, Instar==5)
fit = stan("Ec122_PS7.stan", iter = 10000, chains = 4,
         data = list(log_Mrate=log(new_data$Mrate),
                     log_BodySize=log(new_data$BodySize),
                     Instar=new_data$Instar,
                     N=length(new_data[,1])))
yrep = as.matrix(fit, pars = "yPred")
indices = sample(nrow(yrep),100)
p5 = ppc_dens_overlay(log(new_data$Mrate), yrep[indices,])

grid.arrange(p1, p2, p3, p4, p5, nrow = 3)

