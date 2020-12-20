#### Section 2 Plots ####

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

###############
#### Setup ####
###############

## Prelims ##
n = 80; t = 2;
seed = 1
K = 50;
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

###############

################
#### Plot 1 ####
################

## Parameters ##
param_sig = c(rho0 = 0.5, sigma2 = 0.2, zeta = 0.5, rho1 = 0.05, gamma = 2, 
             alpha = pi/4, mu_x = 0.3, mu_y = -0.3)
param_sig_i_t_dep = lapply(1:length(param_sig), 
                           function(i) param_t_dep_func(1,param_sig[i]))
param_sig_t_dep = 
  param_sig_t_dep_func(param_sig_i_t_dep[[1]],param_sig_i_t_dep[[2]],
                       param_sig_i_t_dep[[3]],param_sig_i_t_dep[[4]],
                       param_sig_i_t_dep[[5]],param_sig_i_t_dep[[6]],
                       param_sig_i_t_dep[[7]],param_sig_i_t_dep[[8]])


## Weighting Matrix (Physical Space) ##
grid = expand.grid(x = 0:(n-1),y = 0:(n-1))/n
centre_x = 2/5; centre_y = 1/2; sigma_x = 0.2; sigma_y = 0.2; rho = 0
weights = 1+6*with(grid,sech2d(x = x,y = y, centre_x = centre_x,
                           centre_y = centre_y, sigma_x = sigma_x,
                           sigma_y=sigma_y, rho=rho))
grid$weights = weights
weights_mat = weights_mat_func(weights)

## Weighting Matrix (Spectral Space) ##
spec_weights_mat = spec_weights_mat_func(n = n, K = K, 
                                         weights_mat = weights_mat)

if(FALSE){
  #pdf("signal_weights.pdf",width=6,height=6)
  wireframe(weights~x*y, data = grid, cuts = 100, main="", 
            col.regions = heat.colors(1000)[1:length(heat.colors(1000))],
            drape = TRUE,zlab="z",xlab="x",ylab="y", light.source = c(10,10,10), 
            aspect = c(1,1), screen = list(z = 20, x = -60), col = F, 
            zoom = 1, scales = list(arrows = F, tck = c(0.8, 0.6, 0.4), 
                                    distance =c(.7, .7, .5)))
  #dev.off()
}


## Loop Over Different Observation Locations ##
x_locs = seq(1,n,1)
y_locs = seq(1,n,1)
n_x = length(x_locs); n_y = length(y_locs)
obj_func = matrix(rep(0,n_x*n_y),nrow=n_x,ncol=n_y)
loop_times = rep(0,n_x*n_y)

## Kalman Filter (independent of obs. locations) ##
mu0 = rep(0,K); Sigma0 = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,
                             spec_weights_mat = spec_weights_mat)
XBeta = XBeta_f(param_bias,n,t,n_obs)
R = R_f(param_obs,n,n_obs)
G = G_f(param_sig,wave,cosine_indices,dt=1,n_cos=n_cos)
Q = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,spec_weights_mat=spec_weights_mat)

for (i in 1:n_x){
  for (j in 1:n_y){
    
    current_iter = (i-1)*n_x + j
    remaining_iter = n_x*n_y - current_iter
    current_time1 = Sys.time()
    
    last_loc = c(x_locs[i],y_locs[j])
    grid = matrix(data = c(c(n/5,n/5,n/5,4*n/5,4*n/5,n/5,4*n/5,4*n/5),last_loc),nrow=5,byrow=T)
    coords = lapply(seq_len(nrow(grid)), function(i) grid[i,])
    coords_indices = sapply(coords,function(x) x[1]+n*x[2]+1)
    n_obs = length(coords)
    obs_loc = lapply(1:n_obs,function(j) matrix(coords[[j]],nrow=t,ncol=2,byrow=T))
    
    ## Simulate ##
    spde_sim = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                               param_bias_t_dep,seed = 1,n_obs = n_obs,
                               obs_loc = obs_loc, weights = weights)
    y_sim = spde_sim$y
    
    ## Kalman Filter ##
    Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs,coords)
    kf = KalmanFilter(y = y_sim, Phi = Phi, R = R, G = G, Q = Q, t = t, K = K, 
                      log_lik = TRUE, m0 = mu0, R0 = Sigma0)
    
    ## Optimal State Estimate ##
    alpha_hat = kf$m_i_i[t+1,]
    #x_hat = FFT_real_matrix(wave,cosine_indices,n_cos,n)%*%alpha_hat
    
    ## Objective Function ##
    sigma_hat = kf$R_i_i[t+1,,]
    obj_func[i,j] = sum(diag(sigma_hat))
    
    current_time2 = Sys.time()
    
    loop_times[current_iter] = as.vector(current_time2-current_time1)
    average_loop_time = mean(loop_times[1:current_iter])
    
    if(current_iter%%100==0){
      cat(remaining_iter, "iterations remaining = ~", round(average_loop_time*remaining_iter/60,3), "mins\n")
    }
    
  }
}

## Save ##
saveRDS(obj_func,"section2_plot1_object")
obj_func = readRDS("section2_plot1_object")

## Plot ##
library("plot3D")
setwd(fig_wd)
pdf("section2_plot1a.pdf",height=7,width=7)
par(mar=c(1,.5,.5,2),mfrow=c(1,1))
#persp3D(z = obj_func,theta = 30,phi=20)
grid = expand.grid(x = 0:(n-1),y = 0:(n-1))/n
grid$obj_func = as.vector(obj_func)
wireframe(obj_func~x*y, data = grid, cuts = 100, main="", 
          col.regions = jet.col(1000)[1:1000],
          drape = TRUE,zlab="z",xlab="x",ylab="y", light.source = c(20,10,10), 
          aspect = c(1,1), screen = list(z = -45, x = -55), col = F, 
          zoom = 1, scales = list(arrows = F, tck = c(0.8, 0.6, 0.4), 
                                  distance =c(.7, .7, .5)))
dev.off()

pdf("section2_plot1b.pdf",height=7,width=7)
par(mar=c(2.6,2.8,1.1,1.1),mfrow=c(1,1),family="LM Roman 10")
#par(mar=c(2.6,2.8,1.1,3),mfrow=c(1,1))
image2D(z = obj_func,main="",xlab="",ylab="",cex.axis=1.5,colkey=F)#colkey = list(cex.axis=2))
points(x = grid$x[which(obj_func==min(obj_func))],y = grid$y[which(obj_func==min(obj_func))],
       col="white",pch=19,cex=4)
points(x = grid$x[which(obj_func==min(obj_func))],y = grid$y[which(obj_func==min(obj_func))],
       col="black",pch=19,cex=2)
dev.off()

################

################
#### Plot 2 ####
################

## Parameters ##
param_sig = c(rho0 = 0.5, sigma2 = 0.2, zeta = 0.5, rho1 = 0.05, gamma = 2, 
              alpha = pi/4, mu_x = 0.3, mu_y = -0.3)
param_sig_i_t_dep = lapply(1:length(param_sig), 
                           function(i) param_t_dep_func(1,param_sig[i]))
param_sig_t_dep = 
  param_sig_t_dep_func(param_sig_i_t_dep[[1]],param_sig_i_t_dep[[2]],
                       param_sig_i_t_dep[[3]],param_sig_i_t_dep[[4]],
                       param_sig_i_t_dep[[5]],param_sig_i_t_dep[[6]],
                       param_sig_i_t_dep[[7]],param_sig_i_t_dep[[8]])


## Weighting Matrix (Physical Space) ##
grid = expand.grid(x = 0:(n-1),y = 0:(n-1))/n
centre_x = 1/5; centre_y = 1/5; sigma_x = 0.2; sigma_y = 0.2; rho = 0
weights = 1+6*with(grid,sech2d(x = x,y = y, centre_x = centre_x,
                             centre_y = centre_y, sigma_x = sigma_x,
                             sigma_y=sigma_y, rho=rho))
grid$weights = weights
weights_mat = weights_mat_func(weights)

## Weighting Matrix (Spectral Space) ##
spec_weights_mat = spec_weights_mat_func(n = n, K = K, 
                                         weights_mat = weights_mat)

if(FALSE){
  #pdf("signal_weights.pdf",width=6,height=6)
  wireframe(weights~x*y, data = grid, cuts = 100, main="", 
            col.regions = heat.colors(1000)[1:length(heat.colors(1000))],
            drape = TRUE,zlab="z",xlab="x",ylab="y", light.source = c(10,10,10), 
            aspect = c(1,1), screen = list(z = 20, x = -60), col = F, 
            zoom = 1, scales = list(arrows = F, tck = c(0.8, 0.6, 0.4), 
                                    distance =c(.7, .7, .5)))
  #dev.off()
}


## Loop Over Different Observation Locations ##
x_locs = seq(1,n,1)
y_locs = seq(1,n,1)
n_x = length(x_locs); n_y = length(y_locs)
obj_func = matrix(rep(0,n_x*n_y),nrow=n_x,ncol=n_y)
loop_times = rep(0,n_x*n_y)

## Kalman Filter (independent of obs. locations) ##
mu0 = rep(0,K); Sigma0 = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,
                             spec_weights_mat = spec_weights_mat)
XBeta = XBeta_f(param_bias,n,t,n_obs)
R = R_f(param_obs,n,n_obs)
G = G_f(param_sig,wave,cosine_indices,dt=1,n_cos=n_cos)
Q = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,spec_weights_mat=spec_weights_mat)

for (i in 1:n_x){
  for (j in 1:n_y){
    
    current_iter = (i-1)*n_x + j
    remaining_iter = n_x*n_y - current_iter
    current_time1 = Sys.time()
    
    last_loc = c(x_locs[i],y_locs[j])
    grid = matrix(data = c(c(n/5,n/5,n/5,4*n/5,4*n/5,n/5,4*n/5,4*n/5),last_loc),nrow=5,byrow=T)
    coords = lapply(seq_len(nrow(grid)), function(i) grid[i,])
    coords_indices = sapply(coords,function(x) x[1]+n*x[2]+1)
    n_obs = length(coords)
    obs_loc = lapply(1:n_obs,function(j) matrix(coords[[j]],nrow=t,ncol=2,byrow=T))
    
    ## Simulate ##
    spde_sim = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                               param_bias_t_dep,seed = 1,n_obs = n_obs,
                               obs_loc = obs_loc, weights = weights)
    y_sim = spde_sim$y
    
    ## Kalman Filter ##
    Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs,coords)
    kf = KalmanFilter(y = y_sim, Phi = Phi, R = R, G = G, Q = Q, t = t, K = K, 
                      log_lik = TRUE, m0 = mu0, R0 = Sigma0)
    
    ## Optimal State Estimate ##
    alpha_hat = kf$m_i_i[t+1,]
    #x_hat = FFT_real_matrix(wave,cosine_indices,n_cos,n)%*%alpha_hat
    
    ## Objective Function ##
    sigma_hat = kf$R_i_i[t+1,,]
    obj_func[i,j] = sum(diag(sigma_hat))
    
    current_time2 = Sys.time()
    
    loop_times[current_iter] = as.vector(current_time2-current_time1)
    average_loop_time = mean(loop_times[1:current_iter])
    
    if(current_iter%%100==0){
      cat(remaining_iter, "iterations remaining = ~", round(average_loop_time*remaining_iter/60,3), "mins\n")
    }
    
  }
}

## Save ##
saveRDS(obj_func,"section2_plot2_object")

## Plot ##
library("plot3D")
setwd(fig_wd)
pdf("section2_plot2a.pdf",height=7,width=7)
par(mar=c(1,.5,.5,2),mfrow=c(1,1))
#persp3D(z = obj_func,theta = 30,phi=20)
grid = expand.grid(x = 0:(n-1),y = 0:(n-1))/n
grid$obj_func = as.vector(obj_func)
wireframe(obj_func~x*y, data = grid, cuts = 100, main="", 
          col.regions = jet.col(1000)[1:1000],
          drape = TRUE,zlab="z",xlab="x",ylab="y", light.source = c(10,10,10), 
          aspect = c(1,1), screen = list(z = -45, x = -55), col = F, 
          zoom = 1, scales = list(arrows = F, tck = c(0.8, 0.6, 0.4), 
                                  distance =c(.7, .7, .5)))
dev.off()

pdf("section2_plot2b.pdf",height=7,width=7)
par(mar=c(2.6,2.8,1.1,1.1),mfrow=c(1,1),family="LM Roman 10")
#par(mar=c(2.6,2.8,1.1,3),mfrow=c(1,1))
image2D(z = obj_func,main="",xlab="",ylab="",cex.axis=1.5,colkey=F)#colkey = list(cex.axis=2))
points(x = grid$x[which(obj_func==min(obj_func))],y = grid$y[which(obj_func==min(obj_func))],
       col="white",pch=19,cex=4)
points(x = grid$x[which(obj_func==min(obj_func))],y = grid$y[which(obj_func==min(obj_func))],
       col="black",pch=19,cex=2)
dev.off()

################