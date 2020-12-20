#### Intro Plots ####

#################
#### PRELIMS ####
#################

## Define directories ##
main_wd = "/Users/ls616/Google Drive/MPE CDT/PhD/Year 1/Code/Current"
fig_wd = "/Users/ls616/Desktop"

## Set directory ##
setwd(main_wd)

## Load functions from C ##
dyn.load("FFBS_spectral.so")
dyn.load("RML_spectral.so")
dyn.load("propagate_spectral.so")
dyn.load("tkf_spectral.so")
dyn.load("tkf2_spectral.so")

## Load additional libraries ##
library(mvtnorm)
library(coda)
library(microbenchmark)
library(lattice)
library(dlm)
library(RColorBrewer)
library(gtools)

## Load spde2019 functions ##
source("spde_v3.R")

#################

###################
#### All Plots ####
###################

## Prelims ##
n = 100; t = 10;
seed = 1
K = 200;
SPDE_FT = spde_initialise(n,t,K)
wave = SPDE_FT$wave
cosine_indices = SPDE_FT$cosine_indices
n_cos = SPDE_FT$n_cos
n_cos_sin = length(cosine_indices)
K = SPDE_FT$K

## Observation parameters ##
m = 1; m_indices = list(coords_indices); tau2_vals = list(0.01)
param_obs = param_obs_func(m,m_indices,tau2_vals)
param_obs_t_dep =  param_t_dep_func(list(1),list(param_obs))

## Observation bias ##
m1 = 0; m1_indices = list(coords_indices); bias_vals = list(0)
param_bias = param_bias_func(m1,m1_indices,bias_vals)
param_bias_t_dep = param_t_dep_func(list(1),list(param_bias))

###################

#############################
#### Plot 1 (true field) ####
#############################

## Parameters ##
param_sig = c(rho0 = 0.5, sigma2 = 0.2, zeta = 0.5,rho1 = 0.05, gamma = 2, 
              alpha = pi/4, mu_x = 0.3, mu_y = -0.3)
param_sig_i_t_dep = lapply(1:length(param_sig), 
                           function(i) param_t_dep_func(1,param_sig[i]))
param_sig_t_dep = 
  param_sig_t_dep_func(param_sig_i_t_dep[[1]],param_sig_i_t_dep[[2]],
                       param_sig_i_t_dep[[3]],param_sig_i_t_dep[[4]],
                       param_sig_i_t_dep[[5]],param_sig_i_t_dep[[6]],
                       param_sig_i_t_dep[[7]],param_sig_i_t_dep[[8]])

## Simulate ##
spde_sim1 = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                            param_bias_t_dep,seed = seed)

## Plot ##
#setwd("/Users/ls616/Desktop")
#pdf("intro_plot_1.pdf",height=7,width=7)
#par(mar=c(3,3,1,2))
plot(spde_sim1,which_t=t,main="NULL")
#dev.off()

#############################

#########################################################################
#### Plot 2 (filter estimate, true parameters, optimal observations) ####
#########################################################################

## True Parameters ##
param_sig = c(rho0 = 0.5, sigma2 = 0.2, zeta = 0.5,rho1 = 0.05, gamma = 2, 
             alpha = pi/4, mu_x = 0.3, mu_y = -0.3)
param_sig_i_t_dep = lapply(1:length(param_sig), 
                           function(i) param_t_dep_func(1,param_sig[i]))
param_sig_t_dep = 
  param_sig_t_dep_func(param_sig_i_t_dep[[1]],param_sig_i_t_dep[[2]],
                       param_sig_i_t_dep[[3]],param_sig_i_t_dep[[4]],
                       param_sig_i_t_dep[[5]],param_sig_i_t_dep[[6]],
                       param_sig_i_t_dep[[7]],param_sig_i_t_dep[[8]])

## Optimal Observations ##
grid = as.matrix(expand.grid(10*seq(0.5,9.5,.5),10*seq(0.5,9.5,.5)))
coords = lapply(seq_len(nrow(grid)), function(i) grid[i,])
coords_indices = sapply(coords,function(x) x[1]+n*x[2]+1)
n_obs = length(coords)
obs_loc = lapply(1:n_obs,function(j) matrix(coords[[j]],nrow=t,ncol=2,byrow=T))

## Simulate ##
spde_sim = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                           param_bias_t_dep,seed = 1,n_obs = n_obs,
                           obs_loc = obs_loc)
y_sim = spde_sim$y

## Kalman Filter ##
mu0 = rep(0,K); Sigma0 = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE)
Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs,coords)
XBeta = XBeta_f(param_bias,n,t,n_obs)
R = R_f(param_obs,n,n_obs)
G = G_f(param_sig,wave,cosine_indices,dt=1,n_cos=n_cos)
Q = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE)
kf = KalmanFilter(y = y_sim, Phi = Phi, R = R, G = G, Q = Q, t = t, K = K, 
             log_lik = TRUE, m0 = mu0, R0 = Sigma0)

## Optimal State Estimate ##
alpha_hat = kf$m_i_i[t+1,]
x_hat = FFT_real_matrix(wave,cosine_indices,n_cos,n)%*%alpha_hat

## Plot ##
#setwd("/Users/ls616/Desktop")
#pdf("intro_plot_2.pdf",height=7,width=7)
#par(mar=c(3,3,1,2))
SPDE_plot(x_i = t(x_hat),dt=1,main="")
#dev.off()

#########################################################################

########################################################################
#### Plot 3 (filter estimate, bad parameters, optimal observations) ####
########################################################################

## Bad Parameters ##
param_sig = c(rho0 = 0.5, sigma2 = 0.2, zeta = 0.5,rho1 = 0.3, gamma = 2, 
              alpha = pi/4, mu_x = 0.3, mu_y = -0.3)
param_sig_i_t_dep = lapply(1:length(param_sig), 
                           function(i) param_t_dep_func(1,param_sig[i]))
param_sig_t_dep = 
  param_sig_t_dep_func(param_sig_i_t_dep[[1]],param_sig_i_t_dep[[2]],
                       param_sig_i_t_dep[[3]],param_sig_i_t_dep[[4]],
                       param_sig_i_t_dep[[5]],param_sig_i_t_dep[[6]],
                       param_sig_i_t_dep[[7]],param_sig_i_t_dep[[8]])

## Optimal Observations ##
grid = as.matrix(expand.grid(10*seq(0.5,9.5,.5),10*seq(0.5,9.5,.5)))
coords = lapply(seq_len(nrow(grid)), function(i) grid[i,])
coords_indices = sapply(coords,function(x) x[1]+n*x[2]+1)
n_obs = length(coords)
obs_loc = lapply(1:n_obs,function(j) matrix(coords[[j]],nrow=t,ncol=2,byrow=T))

## Simulate ##
spde_sim = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                           param_bias_t_dep,seed = 1,n_obs = n_obs,
                           obs_loc = obs_loc)
y_sim = spde_sim$y

## Kalman Filter ##
mu0 = rep(0,K); Sigma0 = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE)
Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs,coords)
XBeta = XBeta_f(param_bias,n,t,n_obs)
R = R_f(param_obs,n,n_obs)
G = G_f(param_sig,wave,cosine_indices,dt=1,n_cos=n_cos)
Q = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE)
kf = KalmanFilter(y = y_sim, Phi = Phi, R = R, G = G, Q = Q, t = t, K = K, 
                  log_lik = TRUE, m0 = mu0, R0 = Sigma0)

## Optimal State Estimate ##
alpha_hat = kf$m_i_i[t+1,]
x_hat = FFT_real_matrix(wave,cosine_indices,n_cos,n)%*%alpha_hat

## Plot ##
setwd("/Users/ls616/Desktop")
pdf("intro_plot_3.pdf",height=7,width=7)
par(mar=c(3,3,1,2))
SPDE_plot(x_i = t(x_hat),dt=1,main="")
dev.off()

#########################################################################

#########################################################################
#### Plot 4 (filter estimate, true parameters, bad observations) ########
#########################################################################

## True Parameters ##
param_sig = c(rho0 = 0.5, sigma2 = 0.2, zeta = 0.5,rho1 = 0.05, gamma = 2, 
              alpha = pi/4, mu_x = 0.3, mu_y = -0.3)
param_sig_i_t_dep = lapply(1:length(param_sig), 
                           function(i) param_t_dep_func(1,param_sig[i]))
param_sig_t_dep = 
  param_sig_t_dep_func(param_sig_i_t_dep[[1]],param_sig_i_t_dep[[2]],
                       param_sig_i_t_dep[[3]],param_sig_i_t_dep[[4]],
                       param_sig_i_t_dep[[5]],param_sig_i_t_dep[[6]],
                       param_sig_i_t_dep[[7]],param_sig_i_t_dep[[8]])

## Bad Observations ##
grid = as.matrix(expand.grid(10*seq(4.1,5.9,0.1),10*seq(4.1,5.9,0.1)))
coords = lapply(seq_len(nrow(grid)), function(i) grid[i,])
coords_indices = sapply(coords,function(x) x[1]+n*x[2]+1)
n_obs = length(coords)
obs_loc = lapply(1:n_obs,function(j) matrix(coords[[j]],nrow=t,ncol=2,byrow=T))

## Simulate ##
spde_sim = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                           param_bias_t_dep,seed = 1,n_obs = n_obs,
                           obs_loc = obs_loc)
y_sim = spde_sim$y

## Kalman Filter ##
mu0 = rep(0,K); Sigma0 = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE)
Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs,coords)
XBeta = XBeta_f(param_bias,n,t,n_obs)
R = R_f(param_obs,n,n_obs)
G = G_f(param_sig,wave,cosine_indices,dt=1,n_cos=n_cos)
Q = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE)
kf = KalmanFilter(y = y_sim, Phi = Phi, R = R, G = G, Q = Q, t = t, K = K, 
                  log_lik = TRUE, m0 = mu0, R0 = Sigma0)

## Optimal State Estimate ##
alpha_hat = kf$m_i_i[t+1,]
x_hat = FFT_real_matrix(wave,cosine_indices,n_cos,n)%*%alpha_hat

## Plot ##
setwd("/Users/ls616/Desktop")
pdf("intro_plot_4.pdf",height=7,width=7)
par(mar=c(3,3,1,2))
SPDE_plot(x_i = t(x_hat),dt=1,main="")
dev.off()

#########################################################################

########################################################################
#### Plot 5 (filter estimate, bad parameters, bad observations) ########
########################################################################

## Bad Parameters ##
param_sig = c(rho0 = 0.5, sigma2 = 0.2, zeta = 0.5,rho1 = 0.3, gamma = 2, 
              alpha = pi/4, mu_x = 0.3, mu_y = -0.3)
param_sig_i_t_dep = lapply(1:length(param_sig), 
                           function(i) param_t_dep_func(1,param_sig[i]))
param_sig_t_dep = 
  param_sig_t_dep_func(param_sig_i_t_dep[[1]],param_sig_i_t_dep[[2]],
                       param_sig_i_t_dep[[3]],param_sig_i_t_dep[[4]],
                       param_sig_i_t_dep[[5]],param_sig_i_t_dep[[6]],
                       param_sig_i_t_dep[[7]],param_sig_i_t_dep[[8]])

## Bad Observations ##
grid = as.matrix(expand.grid(10*seq(4.1,5.9,0.1),10*seq(4.1,5.9,0.1)))
coords = lapply(seq_len(nrow(grid)), function(i) grid[i,])
coords_indices = sapply(coords,function(x) x[1]+n*x[2]+1)
n_obs = length(coords)
obs_loc = lapply(1:n_obs,function(j) matrix(coords[[j]],nrow=t,ncol=2,byrow=T))

## Simulate ##
spde_sim = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                           param_bias_t_dep,seed = 1,n_obs = n_obs,
                           obs_loc = obs_loc)
y_sim = spde_sim$y

## Kalman Filter ##
mu0 = rep(0,K); Sigma0 = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE)
Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs,coords)
XBeta = XBeta_f(param_bias,n,t,n_obs)
R = R_f(param_obs,n,n_obs)
G = G_f(param_sig,wave,cosine_indices,dt=1,n_cos=n_cos)
Q = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE)
kf = KalmanFilter(y = y_sim, Phi = Phi, R = R, G = G, Q = Q, t = t, K = K, 
                  log_lik = TRUE, m0 = mu0, R0 = Sigma0)

## Optimal State Estimate ##
alpha_hat = kf$m_i_i[t+1,]
x_hat = FFT_real_matrix(wave,cosine_indices,n_cos,n)%*%alpha_hat

## Plot ##
setwd("/Users/ls616/Desktop")
pdf("intro_plot_5.pdf",height=7,width=7)
par(mar=c(3,3,1,2))
SPDE_plot(x_i = t(x_hat),dt=1,main="")
dev.off()

########################################################################