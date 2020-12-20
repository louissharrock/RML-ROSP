#### Section 4: Numerical Experiment 2 ####

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
n = 50; t = 50;
seed = 1
K = 21;
SPDE_FT = spde_initialise(n,t,K)
wave = SPDE_FT$wave
cosine_indices = SPDE_FT$cosine_indices
n_cos = SPDE_FT$n_cos
n_cos_sin = length(cosine_indices)
K = SPDE_FT$K

###############


##############################
#### Plot 1 (rho0 = 0.01) ####
##############################

## Parameters ##
param_sig = c(rho0 = 0.01, sigma2 = 0.5, zeta = 0.5,rho1 = 0.1, gamma = 5, 
              alpha = pi/3, mu_x = 0.3, mu_y = -0.3)
param_sig_i_t_dep = lapply(1:length(param_sig), 
                           function(i) param_t_dep_func(1,param_sig[i]))
param_sig_t_dep = 
  param_sig_t_dep_func(param_sig_i_t_dep[[1]],param_sig_i_t_dep[[2]],
                       param_sig_i_t_dep[[3]],param_sig_i_t_dep[[4]],
                       param_sig_i_t_dep[[5]],param_sig_i_t_dep[[6]],
                       param_sig_i_t_dep[[7]],param_sig_i_t_dep[[8]])


## Alternative: No Weighting Matrix ##
weights = NULL; spec_weights_mat = NULL

## FFT real matrix ##
FFT_matrix = FFT_real_matrix(wave,cosine_indices,n_cos,n)

## Loop Over Different Observation Locations ##
x_locs = seq(0,n-1,1)
y_locs = seq(0,n-1,1)
n_x = length(x_locs); n_y = length(y_locs)
obj_func = matrix(rep(0,n_x*n_y),nrow=n_x,ncol=n_y)
loop_times = rep(0,n_x*n_y)

## Kalman Filter (independent of obs. locations) ##
mu0 = rep(0,K); Sigma0 = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,
                             spec_weights_mat = spec_weights_mat)
G = G_f(param_sig,wave,cosine_indices,dt=1,n_cos=n_cos)
Q = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,spec_weights_mat=spec_weights_mat)

for (i in 1:n_x){
  for (j in 1:n_y){
    
    current_iter = (i-1)*n_x + j
    remaining_iter = n_x*n_y - current_iter
    current_time1 = Sys.time()
    
    last_loc = c(x_locs[i],y_locs[j])
    grid = matrix(data = c(c((n-1)/5,(n-1)/5,(n-1)/5,4*(n-1)/5,4*(n-1)/5,(n-1)/5,4*(n-1)/5,4*(n-1)/5),last_loc),nrow=5,byrow=T)
    coords = lapply(seq_len(nrow(grid)), function(i) grid[i,])
    coords_indices = sapply(coords,function(x) x[1]+n*x[2]+1)
    n_obs = length(coords)
    obs_loc = lapply(1:n_obs,function(j) matrix(coords[[j]],nrow=t,ncol=2,byrow=T))
    
    ## Observation parameters ##
    m = 1; m_indices = list(coords_indices); tau2_vals = list(0.01)
    param_obs = param_obs_func(m,m_indices,tau2_vals)
    param_obs_t_dep =  param_t_dep_func(list(1),list(param_obs))
    
    ## Observation bias ##
    m1 = 0; m1_indices = list(coords_indices); bias_vals = list(0)
    param_bias = param_bias_func(m1,m1_indices,bias_vals)
    param_bias_t_dep = param_t_dep_func(list(1),list(param_bias))
    
    ## Simulate ##
    spde_sim = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                               param_bias_t_dep,seed = 1,n_obs = n_obs,
                               obs_loc = obs_loc, weights = weights)
    y_sim = spde_sim$y
    
    ## Kalman Filter ##
    Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs,coords)
    XBeta = XBeta_f(param_bias,n,t,n_obs)
    R = R_f(param_obs,n,n_obs)
    kf = KalmanFilter(y = y_sim, Phi = Phi, R = R, G = G, Q = Q, t = t, K = K, 
                      log_lik = TRUE, m0 = mu0, R0 = Sigma0)
    
    ## Optimal State Estimate ##
    alpha_hat = kf$m_i_i[t+1,]
    #x_hat = FFT_real_matrix(wave,cosine_indices,n_cos,n)%*%alpha_hat
    
    ## Objective Function ##
    sigma_hat = kf$R_i_i[t+1,,]
    #sigma_hat = objective_function(t,kf$R_i_i,diag(1,nrow=K))
    #sigma_hat_true = FFT_matrix%*%sigma_hat%*%t(FFT_matrix) [not needed]
    obj_func[i,j] = sum(diag(sigma_hat))
    
    current_time2 = Sys.time()
    
    loop_times[current_iter] = as.vector(current_time2-current_time1)
    average_loop_time = mean(loop_times[1:current_iter])
    
    if(current_iter%%50==0){
      cat(remaining_iter, "iterations remaining = ~", round(average_loop_time*remaining_iter/60,3), "mins\n")
    }
    
  }
}

## Save ##
setwd(fig_wd)
saveRDS(obj_func,"param_dep_simulation_1")

## Load ##
if(FALSE){
  obj_func = readRDS("param_dep_simulation_1")
}

## Match obj_func values with (x,y) coords ##
grid = expand.grid(x = 0:(n-1),y = 0:(n-1))/n
grid$obj_func = as.vector(obj_func)

## Plot ##
library("plot3D")
library("extrafont")

## Wireframe Plot ##
if(FALSE){
  setwd(fig_wd)
  #pdf("section4_sim2_plot1.pdf",height=7,width=7)
  par(mar=c(1,.5,.5,2),mfrow=c(1,1))
  #persp3D(z = obj_func,theta = 30,phi=20)
  wireframe(obj_func~x*y, data = grid, cuts = 100, main="", 
            col.regions = jet.col(1000)[1:1000],
            drape = TRUE,zlab="z",xlab="x",ylab="y", light.source = c(20,10,10), 
            aspect = c(1,1), screen = list(z = -45, x = -55), col = F, 
            zoom = 1, scales = list(arrows = F, tck = c(0.8, 0.6, 0.4), 
                                    distance =c(.7, .7, .5)))
  #dev.off()
}

## Image Plot ##
if(TRUE){
  pdf("section4_sim2_plot1.pdf",height=7,width=7)
  par(mar=c(0.5,0.5,0.5,0.5),mgp=c(1.5,1,0),mfrow=c(1,1),family="LM Roman 10")
  image2D(z = obj_func, x=seq(0,1-1/n,1/n), y = seq(0,1-1/n,1/n),main="",xlab="",ylab="",cex.axis=1.5,colkey=F,
          xaxs="i",yaxs="i",resfac=10,xaxt="n",yaxt="n")#colkey = list(cex.axis=2))
  points(x = grid$x[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         y = grid$y[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         col="white",pch=19,cex=4)
  points(x = grid$x[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         y = grid$y[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         col="black",pch=19,cex=2)
  dev.off()
}

###############################

##############################
#### Plot 2 (rho0 = 0.03) ####
##############################  

## Parameters ##
param_sig = c(rho0 = 0.03, sigma2 = 0.2, zeta = 0.5,rho1 = 0.1, gamma = 5, 
              alpha = pi/4, mu_x = 0.3, mu_y = -0.3)
param_sig_i_t_dep = lapply(1:length(param_sig), 
                           function(i) param_t_dep_func(1,param_sig[i]))
param_sig_t_dep = 
  param_sig_t_dep_func(param_sig_i_t_dep[[1]],param_sig_i_t_dep[[2]],
                       param_sig_i_t_dep[[3]],param_sig_i_t_dep[[4]],
                       param_sig_i_t_dep[[5]],param_sig_i_t_dep[[6]],
                       param_sig_i_t_dep[[7]],param_sig_i_t_dep[[8]])


## Alternative: No Weighting Matrix ##
weights = NULL; spec_weights_mat = NULL

## FFT real matrix ##
FFT_matrix = FFT_real_matrix(wave,cosine_indices,n_cos,n)

## Loop Over Different Observation Locations ##
x_locs = seq(0,n-1,1)
y_locs = seq(0,n-1,1)
n_x = length(x_locs); n_y = length(y_locs)
obj_func = matrix(rep(0,n_x*n_y),nrow=n_x,ncol=n_y)
loop_times = rep(0,n_x*n_y)

## Kalman Filter (independent of obs. locations) ##
mu0 = rep(0,K); Sigma0 = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,
                             spec_weights_mat = spec_weights_mat)
G = G_f(param_sig,wave,cosine_indices,dt=1,n_cos=n_cos)
Q = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,spec_weights_mat=spec_weights_mat)

for (i in 1:n_x){
  for (j in 1:n_y){
    
    current_iter = (i-1)*n_x + j
    remaining_iter = n_x*n_y - current_iter
    current_time1 = Sys.time()
    
    last_loc = c(x_locs[i],y_locs[j])
    grid = matrix(data = c(c((n-1)/5,(n-1)/5,(n-1)/5,4*(n-1)/5,4*(n-1)/5,(n-1)/5,4*(n-1)/5,4*(n-1)/5),last_loc),nrow=5,byrow=T)
    coords = lapply(seq_len(nrow(grid)), function(i) grid[i,])
    coords_indices = sapply(coords,function(x) x[1]+n*x[2]+1)
    n_obs = length(coords)
    obs_loc = lapply(1:n_obs,function(j) matrix(coords[[j]],nrow=t,ncol=2,byrow=T))
    
    ## Observation parameters ##
    m = 1; m_indices = list(coords_indices); tau2_vals = list(0.01)
    param_obs = param_obs_func(m,m_indices,tau2_vals)
    param_obs_t_dep =  param_t_dep_func(list(1),list(param_obs))
    
    ## Observation bias ##
    m1 = 0; m1_indices = list(coords_indices); bias_vals = list(0)
    param_bias = param_bias_func(m1,m1_indices,bias_vals)
    param_bias_t_dep = param_t_dep_func(list(1),list(param_bias))
    
    ## Simulate ##
    spde_sim = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                               param_bias_t_dep,seed = 1,n_obs = n_obs,
                               obs_loc = obs_loc, weights = weights)
    y_sim = spde_sim$y
    
    ## Kalman Filter ##
    Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs,coords)
    XBeta = XBeta_f(param_bias,n,t,n_obs)
    R = R_f(param_obs,n,n_obs)
    kf = KalmanFilter(y = y_sim, Phi = Phi, R = R, G = G, Q = Q, t = t, K = K, 
                      log_lik = TRUE, m0 = mu0, R0 = Sigma0)
    
    ## Optimal State Estimate ##
    alpha_hat = kf$m_i_i[t+1,]
    #x_hat = FFT_real_matrix(wave,cosine_indices,n_cos,n)%*%alpha_hat
    
    ## Objective Function ##
    sigma_hat = kf$R_i_i[t+1,,]
    #sigma_hat = objective_function(t,kf$R_i_i,diag(1,nrow=K))
    #sigma_hat_true = FFT_matrix%*%sigma_hat%*%t(FFT_matrix) [not needed]
    obj_func[i,j] = sum(diag(sigma_hat))
    
    current_time2 = Sys.time()
    
    loop_times[current_iter] = as.vector(current_time2-current_time1)
    average_loop_time = mean(loop_times[1:current_iter])
    
    if(current_iter%%50==0){
      cat(remaining_iter, "iterations remaining = ~", round(average_loop_time*remaining_iter/60,3), "mins\n")
    }
    
  }
}

## Save ##
setwd(fig_wd)
saveRDS(obj_func,"param_dep_simulation_2")

## Load ##
if(FALSE){
  obj_func = readRDS("param_dep_simulation_2")
}

## Match obj_func values with (x,y) coords ##
grid = expand.grid(x = 0:(n-1),y = 0:(n-1))/n
grid$obj_func = as.vector(obj_func)

## Plot ##
library("plot3D")
library("extrafont")

## Wireframe Plot ##
if(FALSE){
setwd(fig_wd)
#pdf("section4_sim2_plot2.pdf",height=7,width=7)
par(mar=c(1,.5,.5,2),mfrow=c(1,1))
#persp3D(z = obj_func,theta = 30,phi=20)
wireframe(obj_func~x*y, data = grid, cuts = 100, main="", 
          col.regions = jet.col(1000)[1:1000],
          drape = TRUE,zlab="z",xlab="x",ylab="y", light.source = c(20,10,10), 
          aspect = c(1,1), screen = list(z = -45, x = -55), col = F, 
          zoom = 1, scales = list(arrows = F, tck = c(0.8, 0.6, 0.4), 
                                  distance =c(.7, .7, .5)))
#dev.off()
}

## Image Plot ##
if(TRUE){
pdf("section4_sim2_plot2.pdf",height=7,width=7)
par(mar=c(0.5,0.5,0.5,0.5),mgp=c(1.5,1,0),mfrow=c(1,1),family="LM Roman 10")
image2D(z = obj_func, x=seq(0,1-1/n,1/n), y = seq(0,1-1/n,1/n),main="",xlab="",ylab="",cex.axis=1.5,colkey=F,
        xaxs="i",yaxs="i",resfac=10,xaxt="n",yaxt="n")#colkey = list(cex.axis=2))
points(x = grid$x[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
       y = grid$y[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
       col="white",pch=19,cex=4)
points(x = grid$x[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
       y = grid$y[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
       col="black",pch=19,cex=2)
dev.off()
}

##############################

if(FALSE){
##############################
#### Plot 3 (rho0 = 0.05) ####
##############################

## Parameters ##
param_sig = c(rho0 = 0.05, sigma2 = 0.2, zeta = 0.5,rho1 = 0.1, gamma = 5, 
              alpha = pi/4, mu_x = 0.3, mu_y = -0.3)
param_sig_i_t_dep = lapply(1:length(param_sig), 
                           function(i) param_t_dep_func(1,param_sig[i]))
param_sig_t_dep = 
  param_sig_t_dep_func(param_sig_i_t_dep[[1]],param_sig_i_t_dep[[2]],
                       param_sig_i_t_dep[[3]],param_sig_i_t_dep[[4]],
                       param_sig_i_t_dep[[5]],param_sig_i_t_dep[[6]],
                       param_sig_i_t_dep[[7]],param_sig_i_t_dep[[8]])


## Alternative: No Weighting Matrix ##
weights = NULL; spec_weights_mat = NULL

## FFT real matrix ##
FFT_matrix = FFT_real_matrix(wave,cosine_indices,n_cos,n)

## Loop Over Different Observation Locations ##
x_locs = seq(0,n-1,1)
y_locs = seq(0,n-1,1)
n_x = length(x_locs); n_y = length(y_locs)
obj_func = matrix(rep(0,n_x*n_y),nrow=n_x,ncol=n_y)
loop_times = rep(0,n_x*n_y)

## Kalman Filter (independent of obs. locations) ##
mu0 = rep(0,K); Sigma0 = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,
                             spec_weights_mat = spec_weights_mat)
G = G_f(param_sig,wave,cosine_indices,dt=1,n_cos=n_cos)
Q = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,spec_weights_mat=spec_weights_mat)

for (i in 1:n_x){
  for (j in 1:n_y){
    
    current_iter = (i-1)*n_x + j
    remaining_iter = n_x*n_y - current_iter
    current_time1 = Sys.time()
    
    last_loc = c(x_locs[i],y_locs[j])
    grid = matrix(data = c(c((n-1)/5,(n-1)/5,(n-1)/5,4*(n-1)/5,4*(n-1)/5,(n-1)/5,4*(n-1)/5,4*(n-1)/5),last_loc),nrow=5,byrow=T)
    coords = lapply(seq_len(nrow(grid)), function(i) grid[i,])
    coords_indices = sapply(coords,function(x) x[1]+n*x[2]+1)
    n_obs = length(coords)
    obs_loc = lapply(1:n_obs,function(j) matrix(coords[[j]],nrow=t,ncol=2,byrow=T))
    
    ## Observation parameters ##
    m = 1; m_indices = list(coords_indices); tau2_vals = list(0.01)
    param_obs = param_obs_func(m,m_indices,tau2_vals)
    param_obs_t_dep =  param_t_dep_func(list(1),list(param_obs))
    
    ## Observation bias ##
    m1 = 0; m1_indices = list(coords_indices); bias_vals = list(0)
    param_bias = param_bias_func(m1,m1_indices,bias_vals)
    param_bias_t_dep = param_t_dep_func(list(1),list(param_bias))
    
    ## Simulate ##
    spde_sim = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                               param_bias_t_dep,seed = 1,n_obs = n_obs,
                               obs_loc = obs_loc, weights = weights)
    y_sim = spde_sim$y
    
    ## Kalman Filter ##
    Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs,coords)
    XBeta = XBeta_f(param_bias,n,t,n_obs)
    R = R_f(param_obs,n,n_obs)
    kf = KalmanFilter(y = y_sim, Phi = Phi, R = R, G = G, Q = Q, t = t, K = K, 
                      log_lik = TRUE, m0 = mu0, R0 = Sigma0)
    
    ## Optimal State Estimate ##
    alpha_hat = kf$m_i_i[t+1,]
    #x_hat = FFT_real_matrix(wave,cosine_indices,n_cos,n)%*%alpha_hat
    
    ## Objective Function ##
    sigma_hat = kf$R_i_i[t+1,,]
    #sigma_hat = objective_function(t,kf$R_i_i,diag(1,nrow=K))
    #sigma_hat_true = FFT_matrix%*%sigma_hat%*%t(FFT_matrix) [not needed]
    obj_func[i,j] = sum(diag(sigma_hat))
    
    current_time2 = Sys.time()
    
    loop_times[current_iter] = as.vector(current_time2-current_time1)
    average_loop_time = mean(loop_times[1:current_iter])
    
    if(current_iter%%50==0){
      cat(remaining_iter, "iterations remaining = ~", round(average_loop_time*remaining_iter/60,3), "mins\n")
    }
    
  }
}

## Save ##
setwd(fig_wd)
saveRDS(obj_func,"param_dep_simulation_3")

## Load ##
if(FALSE){
  obj_func = readRDS("param_dep_simulation_3")
}

## Match obj_func values with (x,y) coords ##
grid = expand.grid(x = 0:(n-1),y = 0:(n-1))/n
grid$obj_func = as.vector(obj_func)

## Plot ##
library("plot3D")
library("extrafont")

## Wireframe Plot ##
if(FALSE){
  setwd(fig_wd)
  #pdf("section4_sim2_plot3.pdf",height=7,width=7)
  par(mar=c(1,.5,.5,2),mfrow=c(1,1))
  #persp3D(z = obj_func,theta = 30,phi=20)
  wireframe(obj_func~x*y, data = grid, cuts = 100, main="", 
            col.regions = jet.col(1000)[1:1000],
            drape = TRUE,zlab="z",xlab="x",ylab="y", light.source = c(20,10,10), 
            aspect = c(1,1), screen = list(z = -45, x = -55), col = F, 
            zoom = 1, scales = list(arrows = F, tck = c(0.8, 0.6, 0.4), 
                                    distance =c(.7, .7, .5)))
  #dev.off()
}

## Image Plot ##
if(TRUE){
  pdf("section4_sim2_plot3.pdf",height=7,width=7)
  par(mar=c(0.5,0.5,0.5,0.5),mgp=c(1.5,1,0),mfrow=c(1,1),family="LM Roman 10")
  image2D(z = obj_func, x=seq(0,1-1/n,1/n), y = seq(0,1-1/n,1/n),main="",xlab="",ylab="",cex.axis=1.5,colkey=F,
          xaxs="i",yaxs="i",resfac=10,xaxt="n",yaxt="n")#colkey = list(cex.axis=2))
  points(x = grid$x[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         y = grid$y[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         col="white",pch=19,cex=4)
  points(x = grid$x[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         y = grid$y[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         col="black",pch=19,cex=2)
  dev.off()
}

##############################
}

if(FALSE){
##############################
#### Plot 4 (rho0 = 0.08) ####
##############################

## Parameters ##
param_sig = c(rho0 = 0.08, sigma2 = 0.2, zeta = 0.5,rho1 = 0.1, gamma = 5, 
              alpha = pi/4, mu_x = 0.3, mu_y = -0.3)
param_sig_i_t_dep = lapply(1:length(param_sig), 
                           function(i) param_t_dep_func(1,param_sig[i]))
param_sig_t_dep = 
  param_sig_t_dep_func(param_sig_i_t_dep[[1]],param_sig_i_t_dep[[2]],
                       param_sig_i_t_dep[[3]],param_sig_i_t_dep[[4]],
                       param_sig_i_t_dep[[5]],param_sig_i_t_dep[[6]],
                       param_sig_i_t_dep[[7]],param_sig_i_t_dep[[8]])


## Alternative: No Weighting Matrix ##
weights = NULL; spec_weights_mat = NULL

## FFT real matrix ##
FFT_matrix = FFT_real_matrix(wave,cosine_indices,n_cos,n)

## Loop Over Different Observation Locations ##
x_locs = seq(0,n-1,1)
y_locs = seq(0,n-1,1)
n_x = length(x_locs); n_y = length(y_locs)
obj_func = matrix(rep(0,n_x*n_y),nrow=n_x,ncol=n_y)
loop_times = rep(0,n_x*n_y)

## Kalman Filter (independent of obs. locations) ##
mu0 = rep(0,K); Sigma0 = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,
                             spec_weights_mat = spec_weights_mat)
G = G_f(param_sig,wave,cosine_indices,dt=1,n_cos=n_cos)
Q = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,spec_weights_mat=spec_weights_mat)

for (i in 1:n_x){
  for (j in 1:n_y){
    
    current_iter = (i-1)*n_x + j
    remaining_iter = n_x*n_y - current_iter
    current_time1 = Sys.time()
    
    last_loc = c(x_locs[i],y_locs[j])
    grid = matrix(data = c(c((n-1)/5,(n-1)/5,(n-1)/5,4*(n-1)/5,4*(n-1)/5,(n-1)/5,4*(n-1)/5,4*(n-1)/5),last_loc),nrow=5,byrow=T)
    coords = lapply(seq_len(nrow(grid)), function(i) grid[i,])
    coords_indices = sapply(coords,function(x) x[1]+n*x[2]+1)
    n_obs = length(coords)
    obs_loc = lapply(1:n_obs,function(j) matrix(coords[[j]],nrow=t,ncol=2,byrow=T))
    
    ## Observation parameters ##
    m = 1; m_indices = list(coords_indices); tau2_vals = list(0.01)
    param_obs = param_obs_func(m,m_indices,tau2_vals)
    param_obs_t_dep =  param_t_dep_func(list(1),list(param_obs))
    
    ## Observation bias ##
    m1 = 0; m1_indices = list(coords_indices); bias_vals = list(0)
    param_bias = param_bias_func(m1,m1_indices,bias_vals)
    param_bias_t_dep = param_t_dep_func(list(1),list(param_bias))
    
    ## Simulate ##
    spde_sim = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                               param_bias_t_dep,seed = 1,n_obs = n_obs,
                               obs_loc = obs_loc, weights = weights)
    y_sim = spde_sim$y
    
    ## Kalman Filter ##
    Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs,coords)
    XBeta = XBeta_f(param_bias,n,t,n_obs)
    R = R_f(param_obs,n,n_obs)
    kf = KalmanFilter(y = y_sim, Phi = Phi, R = R, G = G, Q = Q, t = t, K = K, 
                      log_lik = TRUE, m0 = mu0, R0 = Sigma0)
    
    ## Optimal State Estimate ##
    alpha_hat = kf$m_i_i[t+1,]
    #x_hat = FFT_real_matrix(wave,cosine_indices,n_cos,n)%*%alpha_hat
    
    ## Objective Function ##
    sigma_hat = kf$R_i_i[t+1,,]
    #sigma_hat = objective_function(t,kf$R_i_i,diag(1,nrow=K))
    #sigma_hat_true = FFT_matrix%*%sigma_hat%*%t(FFT_matrix) [not needed]
    obj_func[i,j] = sum(diag(sigma_hat))
    
    current_time2 = Sys.time()
    
    loop_times[current_iter] = as.vector(current_time2-current_time1)
    average_loop_time = mean(loop_times[1:current_iter])
    
    if(current_iter%%50==0){
      cat(remaining_iter, "iterations remaining = ~", round(average_loop_time*remaining_iter/60,3), "mins\n")
    }
    
  }
}

## Save ##
setwd(fig_wd)
saveRDS(obj_func,"param_dep_simulation_4")

## Load ##
if(FALSE){
  obj_func = readRDS("param_dep_simulation_4")
}

## Match obj_func values with (x,y) coords ##
grid = expand.grid(x = 0:(n-1),y = 0:(n-1))/n
grid$obj_func = as.vector(obj_func)

## Plot ##
library("plot3D")
library("extrafont")

## Wireframe Plot ##
if(FALSE){
  setwd(fig_wd)
  #pdf("section4_sim2_plot4.pdf",height=7,width=7)
  par(mar=c(1,.5,.5,2),mfrow=c(1,1))
  #persp3D(z = obj_func,theta = 30,phi=20)
  wireframe(obj_func~x*y, data = grid, cuts = 100, main="", 
            col.regions = jet.col(1000)[1:1000],
            drape = TRUE,zlab="z",xlab="x",ylab="y", light.source = c(20,10,10), 
            aspect = c(1,1), screen = list(z = -45, x = -55), col = F, 
            zoom = 1, scales = list(arrows = F, tck = c(0.8, 0.6, 0.4), 
                                    distance =c(.7, .7, .5)))
  #dev.off()
}

## Image Plot ##
if(TRUE){
  pdf("section4_sim2_plot4.pdf",height=7,width=7)
  par(mar=c(0.5,0.5,0.5,0.5),mgp=c(1.5,1,0),mfrow=c(1,1),family="LM Roman 10")
  image2D(z = obj_func, x=seq(0,1-1/n,1/n), y = seq(0,1-1/n,1/n),main="",xlab="",ylab="",cex.axis=1.5,colkey=F,
          xaxs="i",yaxs="i",resfac=10,xaxt="n",yaxt="n")#colkey = list(cex.axis=2))
  points(x = grid$x[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         y = grid$y[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         col="white",pch=19,cex=4)
  points(x = grid$x[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         y = grid$y[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         col="black",pch=19,cex=2)
  dev.off()
}

##############################
}

##############################
#### Plot 5 (rho0 = 0.10) ####
##############################

## Parameters ##
param_sig = c(rho0 = 0.10, sigma2 = 0.2, zeta = 0.5,rho1 = 0.1, gamma = 5, 
              alpha = pi/4, mu_x = 0.3, mu_y = -0.3)
param_sig_i_t_dep = lapply(1:length(param_sig), 
                           function(i) param_t_dep_func(1,param_sig[i]))
param_sig_t_dep = 
  param_sig_t_dep_func(param_sig_i_t_dep[[1]],param_sig_i_t_dep[[2]],
                       param_sig_i_t_dep[[3]],param_sig_i_t_dep[[4]],
                       param_sig_i_t_dep[[5]],param_sig_i_t_dep[[6]],
                       param_sig_i_t_dep[[7]],param_sig_i_t_dep[[8]])


## Alternative: No Weighting Matrix ##
weights = NULL; spec_weights_mat = NULL

## FFT real matrix ##
FFT_matrix = FFT_real_matrix(wave,cosine_indices,n_cos,n)

## Loop Over Different Observation Locations ##
x_locs = seq(0,n-1,1)
y_locs = seq(0,n-1,1)
n_x = length(x_locs); n_y = length(y_locs)
obj_func = matrix(rep(0,n_x*n_y),nrow=n_x,ncol=n_y)
loop_times = rep(0,n_x*n_y)

## Kalman Filter (independent of obs. locations) ##
mu0 = rep(0,K); Sigma0 = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,
                             spec_weights_mat = spec_weights_mat)
G = G_f(param_sig,wave,cosine_indices,dt=1,n_cos=n_cos)
Q = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,spec_weights_mat=spec_weights_mat)

for (i in 1:n_x){
  for (j in 1:n_y){
    
    current_iter = (i-1)*n_x + j
    remaining_iter = n_x*n_y - current_iter
    current_time1 = Sys.time()
    
    last_loc = c(x_locs[i],y_locs[j])
    grid = matrix(data = c(c((n-1)/5,(n-1)/5,(n-1)/5,4*(n-1)/5,4*(n-1)/5,(n-1)/5,4*(n-1)/5,4*(n-1)/5),last_loc),nrow=5,byrow=T)
    coords = lapply(seq_len(nrow(grid)), function(i) grid[i,])
    coords_indices = sapply(coords,function(x) x[1]+n*x[2]+1)
    n_obs = length(coords)
    obs_loc = lapply(1:n_obs,function(j) matrix(coords[[j]],nrow=t,ncol=2,byrow=T))
    
    ## Observation parameters ##
    m = 1; m_indices = list(coords_indices); tau2_vals = list(0.01)
    param_obs = param_obs_func(m,m_indices,tau2_vals)
    param_obs_t_dep =  param_t_dep_func(list(1),list(param_obs))
    
    ## Observation bias ##
    m1 = 0; m1_indices = list(coords_indices); bias_vals = list(0)
    param_bias = param_bias_func(m1,m1_indices,bias_vals)
    param_bias_t_dep = param_t_dep_func(list(1),list(param_bias))
    
    ## Simulate ##
    spde_sim = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                               param_bias_t_dep,seed = 1,n_obs = n_obs,
                               obs_loc = obs_loc, weights = weights)
    y_sim = spde_sim$y
    
    ## Kalman Filter ##
    Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs,coords)
    XBeta = XBeta_f(param_bias,n,t,n_obs)
    R = R_f(param_obs,n,n_obs)
    kf = KalmanFilter(y = y_sim, Phi = Phi, R = R, G = G, Q = Q, t = t, K = K, 
                      log_lik = TRUE, m0 = mu0, R0 = Sigma0)
    
    ## Optimal State Estimate ##
    alpha_hat = kf$m_i_i[t+1,]
    #x_hat = FFT_real_matrix(wave,cosine_indices,n_cos,n)%*%alpha_hat
    
    ## Objective Function ##
    sigma_hat = kf$R_i_i[t+1,,]
    #sigma_hat = objective_function(t,kf$R_i_i,diag(1,nrow=K))
    #sigma_hat_true = FFT_matrix%*%sigma_hat%*%t(FFT_matrix) [not needed]
    obj_func[i,j] = sum(diag(sigma_hat))
    
    current_time2 = Sys.time()
    
    loop_times[current_iter] = as.vector(current_time2-current_time1)
    average_loop_time = mean(loop_times[1:current_iter])
    
    if(current_iter%%50==0){
      cat(remaining_iter, "iterations remaining = ~", round(average_loop_time*remaining_iter/60,3), "mins\n")
    }
    
  }
}

## Save ##
setwd(fig_wd)
saveRDS(obj_func,"param_dep_simulation_5")

## Load ##
if(FALSE){
  obj_func = readRDS("param_dep_simulation_5")
}

## Match obj_func values with (x,y) coords ##
grid = expand.grid(x = 0:(n-1),y = 0:(n-1))/n
grid$obj_func = as.vector(obj_func)

## Plot ##
library("plot3D")
library("extrafont")

## Wireframe Plot ##
if(FALSE){
  setwd(fig_wd)
  #pdf("section4_sim2_plot5.pdf",height=7,width=7)
  par(mar=c(1,.5,.5,2),mfrow=c(1,1))
  #persp3D(z = obj_func,theta = 30,phi=20)
  wireframe(obj_func~x*y, data = grid, cuts = 100, main="", 
            col.regions = jet.col(1000)[1:1000],
            drape = TRUE,zlab="z",xlab="x",ylab="y", light.source = c(20,10,10), 
            aspect = c(1,1), screen = list(z = -45, x = -55), col = F, 
            zoom = 1, scales = list(arrows = F, tck = c(0.8, 0.6, 0.4), 
                                    distance =c(.7, .7, .5)))
  #dev.off()
}

## Image Plot ##
if(TRUE){
  pdf("section4_sim2_plot5.pdf",height=7,width=7)
  par(mar=c(0.5,0.5,0.5,0.5),mgp=c(1.5,1,0),mfrow=c(1,1),family="LM Roman 10")
  image2D(z = obj_func, x=seq(0,1-1/n,1/n), y = seq(0,1-1/n,1/n),main="",xlab="",ylab="",cex.axis=1.5,colkey=F,
          xaxs="i",yaxs="i",resfac=10,xaxt="n",yaxt="n")#colkey = list(cex.axis=2))
  points(x = grid$x[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         y = grid$y[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         col="white",pch=19,cex=4)
  points(x = grid$x[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         y = grid$y[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         col="black",pch=19,cex=2)
  dev.off()
}

##############################

if(FALSE){
##############################
#### Plot 6 (rho0 = 0.12) ####
##############################

## Parameters ##
param_sig = c(rho0 = 0.12, sigma2 = 0.2, zeta = 0.5,rho1 = 0.1, gamma = 5, 
              alpha = pi/4, mu_x = 0.3, mu_y = -0.3)
param_sig_i_t_dep = lapply(1:length(param_sig), 
                           function(i) param_t_dep_func(1,param_sig[i]))
param_sig_t_dep = 
  param_sig_t_dep_func(param_sig_i_t_dep[[1]],param_sig_i_t_dep[[2]],
                       param_sig_i_t_dep[[3]],param_sig_i_t_dep[[4]],
                       param_sig_i_t_dep[[5]],param_sig_i_t_dep[[6]],
                       param_sig_i_t_dep[[7]],param_sig_i_t_dep[[8]])


## Alternative: No Weighting Matrix ##
weights = NULL; spec_weights_mat = NULL

## FFT real matrix ##
FFT_matrix = FFT_real_matrix(wave,cosine_indices,n_cos,n)

## Loop Over Different Observation Locations ##
x_locs = seq(0,n-1,1)
y_locs = seq(0,n-1,1)
n_x = length(x_locs); n_y = length(y_locs)
obj_func = matrix(rep(0,n_x*n_y),nrow=n_x,ncol=n_y)
loop_times = rep(0,n_x*n_y)

## Kalman Filter (independent of obs. locations) ##
mu0 = rep(0,K); Sigma0 = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,
                             spec_weights_mat = spec_weights_mat)
G = G_f(param_sig,wave,cosine_indices,dt=1,n_cos=n_cos)
Q = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,spec_weights_mat=spec_weights_mat)

for (i in 1:n_x){
  for (j in 1:n_y){
    
    current_iter = (i-1)*n_x + j
    remaining_iter = n_x*n_y - current_iter
    current_time1 = Sys.time()
    
    last_loc = c(x_locs[i],y_locs[j])
    grid = matrix(data = c(c((n-1)/5,(n-1)/5,(n-1)/5,4*(n-1)/5,4*(n-1)/5,(n-1)/5,4*(n-1)/5,4*(n-1)/5),last_loc),nrow=5,byrow=T)
    coords = lapply(seq_len(nrow(grid)), function(i) grid[i,])
    coords_indices = sapply(coords,function(x) x[1]+n*x[2]+1)
    n_obs = length(coords)
    obs_loc = lapply(1:n_obs,function(j) matrix(coords[[j]],nrow=t,ncol=2,byrow=T))
    
    ## Observation parameters ##
    m = 1; m_indices = list(coords_indices); tau2_vals = list(0.01)
    param_obs = param_obs_func(m,m_indices,tau2_vals)
    param_obs_t_dep =  param_t_dep_func(list(1),list(param_obs))
    
    ## Observation bias ##
    m1 = 0; m1_indices = list(coords_indices); bias_vals = list(0)
    param_bias = param_bias_func(m1,m1_indices,bias_vals)
    param_bias_t_dep = param_t_dep_func(list(1),list(param_bias))
    
    ## Simulate ##
    spde_sim = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                               param_bias_t_dep,seed = 1,n_obs = n_obs,
                               obs_loc = obs_loc, weights = weights)
    y_sim = spde_sim$y
    
    ## Kalman Filter ##
    Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs,coords)
    XBeta = XBeta_f(param_bias,n,t,n_obs)
    R = R_f(param_obs,n,n_obs)
    kf = KalmanFilter(y = y_sim, Phi = Phi, R = R, G = G, Q = Q, t = t, K = K, 
                      log_lik = TRUE, m0 = mu0, R0 = Sigma0)
    
    ## Optimal State Estimate ##
    alpha_hat = kf$m_i_i[t+1,]
    #x_hat = FFT_real_matrix(wave,cosine_indices,n_cos,n)%*%alpha_hat
    
    ## Objective Function ##
    sigma_hat = kf$R_i_i[t+1,,]
    #sigma_hat = objective_function(t,kf$R_i_i,diag(1,nrow=K))
    #sigma_hat_true = FFT_matrix%*%sigma_hat%*%t(FFT_matrix) [not needed]
    obj_func[i,j] = sum(diag(sigma_hat))
    
    current_time2 = Sys.time()
    
    loop_times[current_iter] = as.vector(current_time2-current_time1)
    average_loop_time = mean(loop_times[1:current_iter])
    
    if(current_iter%%50==0){
      cat(remaining_iter, "iterations remaining = ~", round(average_loop_time*remaining_iter/60,3), "mins\n")
    }
    
  }
}

## Save ##
setwd(fig_wd)
saveRDS(obj_func,"param_dep_simulation_6")

## Load ##
if(FALSE){
  obj_func = readRDS("param_dep_simulation_6")
}

## Match obj_func values with (x,y) coords ##
grid = expand.grid(x = 0:(n-1),y = 0:(n-1))/n
grid$obj_func = as.vector(obj_func)

## Plot ##
library("plot3D")
library("extrafont")

## Wireframe Plot ##
if(FALSE){
  setwd(fig_wd)
  #pdf("section4_sim2_plot6.pdf",height=7,width=7)
  par(mar=c(1,.5,.5,2),mfrow=c(1,1))
  #persp3D(z = obj_func,theta = 30,phi=20)
  wireframe(obj_func~x*y, data = grid, cuts = 100, main="", 
            col.regions = jet.col(1000)[1:1000],
            drape = TRUE,zlab="z",xlab="x",ylab="y", light.source = c(20,10,10), 
            aspect = c(1,1), screen = list(z = -45, x = -55), col = F, 
            zoom = 1, scales = list(arrows = F, tck = c(0.8, 0.6, 0.4), 
                                    distance =c(.7, .7, .5)))
  #dev.off()
}

## Image Plot ##
if(TRUE){
  pdf("section4_sim2_plot6.pdf",height=7,width=7)
  par(mar=c(0.5,0.5,0.5,0.5),mgp=c(1.5,1,0),mfrow=c(1,1),family="LM Roman 10")
  image2D(z = obj_func, x=seq(0,1-1/n,1/n), y = seq(0,1-1/n,1/n),main="",xlab="",ylab="",cex.axis=1.5,colkey=F,
          xaxs="i",yaxs="i",resfac=10,xaxt="n",yaxt="n")#colkey = list(cex.axis=2))
  points(x = grid$x[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         y = grid$y[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         col="white",pch=19,cex=4)
  points(x = grid$x[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         y = grid$y[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         col="black",pch=19,cex=2)
  dev.off()
}

##############################
}

##############################
#### Plot 7 (rho0 = 0.15) ####
##############################

## Parameters ##
param_sig = c(rho0 = 0.15, sigma2 = 0.2, zeta = 0.5,rho1 = 0.1, gamma = 5, 
              alpha = pi/4, mu_x = 0.3, mu_y = -0.3)
param_sig_i_t_dep = lapply(1:length(param_sig), 
                           function(i) param_t_dep_func(1,param_sig[i]))
param_sig_t_dep = 
  param_sig_t_dep_func(param_sig_i_t_dep[[1]],param_sig_i_t_dep[[2]],
                       param_sig_i_t_dep[[3]],param_sig_i_t_dep[[4]],
                       param_sig_i_t_dep[[5]],param_sig_i_t_dep[[6]],
                       param_sig_i_t_dep[[7]],param_sig_i_t_dep[[8]])


## Alternative: No Weighting Matrix ##
weights = NULL; spec_weights_mat = NULL

## FFT real matrix ##
FFT_matrix = FFT_real_matrix(wave,cosine_indices,n_cos,n)

## Loop Over Different Observation Locations ##
x_locs = seq(0,n-1,1)
y_locs = seq(0,n-1,1)
n_x = length(x_locs); n_y = length(y_locs)
obj_func = matrix(rep(0,n_x*n_y),nrow=n_x,ncol=n_y)
loop_times = rep(0,n_x*n_y)

## Kalman Filter (independent of obs. locations) ##
mu0 = rep(0,K); Sigma0 = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,
                             spec_weights_mat = spec_weights_mat)
G = G_f(param_sig,wave,cosine_indices,dt=1,n_cos=n_cos)
Q = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,spec_weights_mat=spec_weights_mat)

for (i in 1:n_x){
  for (j in 1:n_y){
    
    current_iter = (i-1)*n_x + j
    remaining_iter = n_x*n_y - current_iter
    current_time1 = Sys.time()
    
    last_loc = c(x_locs[i],y_locs[j])
    grid = matrix(data = c(c((n-1)/5,(n-1)/5,(n-1)/5,4*(n-1)/5,4*(n-1)/5,(n-1)/5,4*(n-1)/5,4*(n-1)/5),last_loc),nrow=5,byrow=T)
    coords = lapply(seq_len(nrow(grid)), function(i) grid[i,])
    coords_indices = sapply(coords,function(x) x[1]+n*x[2]+1)
    n_obs = length(coords)
    obs_loc = lapply(1:n_obs,function(j) matrix(coords[[j]],nrow=t,ncol=2,byrow=T))
    
    ## Observation parameters ##
    m = 1; m_indices = list(coords_indices); tau2_vals = list(0.01)
    param_obs = param_obs_func(m,m_indices,tau2_vals)
    param_obs_t_dep =  param_t_dep_func(list(1),list(param_obs))
    
    ## Observation bias ##
    m1 = 0; m1_indices = list(coords_indices); bias_vals = list(0)
    param_bias = param_bias_func(m1,m1_indices,bias_vals)
    param_bias_t_dep = param_t_dep_func(list(1),list(param_bias))
    
    ## Simulate ##
    spde_sim = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                               param_bias_t_dep,seed = 1,n_obs = n_obs,
                               obs_loc = obs_loc, weights = weights)
    y_sim = spde_sim$y
    
    ## Kalman Filter ##
    Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs,coords)
    XBeta = XBeta_f(param_bias,n,t,n_obs)
    R = R_f(param_obs,n,n_obs)
    kf = KalmanFilter(y = y_sim, Phi = Phi, R = R, G = G, Q = Q, t = t, K = K, 
                      log_lik = TRUE, m0 = mu0, R0 = Sigma0)
    
    ## Optimal State Estimate ##
    alpha_hat = kf$m_i_i[t+1,]
    #x_hat = FFT_real_matrix(wave,cosine_indices,n_cos,n)%*%alpha_hat
    
    ## Objective Function ##
    sigma_hat = kf$R_i_i[t+1,,]
    #sigma_hat = objective_function(t,kf$R_i_i,diag(1,nrow=K))
    #sigma_hat_true = FFT_matrix%*%sigma_hat%*%t(FFT_matrix) [not needed]
    obj_func[i,j] = sum(diag(sigma_hat))
    
    current_time2 = Sys.time()
    
    loop_times[current_iter] = as.vector(current_time2-current_time1)
    average_loop_time = mean(loop_times[1:current_iter])
    
    if(current_iter%%50==0){
      cat(remaining_iter, "iterations remaining = ~", round(average_loop_time*remaining_iter/60,3), "mins\n")
    }
    
  }
}

## Save ##
setwd(fig_wd)
saveRDS(obj_func,"param_dep_simulation_7")

## Load ##
if(FALSE){
  obj_func = readRDS("param_dep_simulation_7")
}

## Match obj_func values with (x,y) coords ##
grid = expand.grid(x = 0:(n-1),y = 0:(n-1))/n
grid$obj_func = as.vector(obj_func)

## Plot ##
library("plot3D")
library("extrafont")

## Wireframe Plot ##
if(FALSE){
  setwd(fig_wd)
  #pdf("section4_sim2_plot7.pdf",height=7,width=7)
  par(mar=c(1,.5,.5,2),mfrow=c(1,1))
  #persp3D(z = obj_func,theta = 30,phi=20)
  wireframe(obj_func~x*y, data = grid, cuts = 100, main="", 
            col.regions = jet.col(1000)[1:1000],
            drape = TRUE,zlab="z",xlab="x",ylab="y", light.source = c(20,10,10), 
            aspect = c(1,1), screen = list(z = -45, x = -55), col = F, 
            zoom = 1, scales = list(arrows = F, tck = c(0.8, 0.6, 0.4), 
                                    distance =c(.7, .7, .5)))
  #dev.off()
}

## Image Plot ##
if(TRUE){
  pdf("section4_sim2_plot7.pdf",height=7,width=7)
  par(mar=c(0.5,0.5,0.5,0.5),mgp=c(1.5,1,0),mfrow=c(1,1),family="LM Roman 10")
  image2D(z = obj_func, x=seq(0,1-1/n,1/n), y = seq(0,1-1/n,1/n),main="",xlab="",ylab="",cex.axis=1.5,colkey=F,
          xaxs="i",yaxs="i",resfac=10,xaxt="n",yaxt="n")#colkey = list(cex.axis=2))
  points(x = grid$x[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         y = grid$y[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         col="white",pch=19,cex=4)
  points(x = grid$x[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         y = grid$y[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         col="black",pch=19,cex=2)
  dev.off()
}

##############################

if(FALSE){
##############################
#### Plot 8 (rho0 = 0.17) ####
##############################

## Parameters ##
param_sig = c(rho0 = 0.17, sigma2 = 0.2, zeta = 0.5,rho1 = 0.1, gamma = 5, 
              alpha = pi/4, mu_x = 0.3, mu_y = -0.3)
param_sig_i_t_dep = lapply(1:length(param_sig), 
                           function(i) param_t_dep_func(1,param_sig[i]))
param_sig_t_dep = 
  param_sig_t_dep_func(param_sig_i_t_dep[[1]],param_sig_i_t_dep[[2]],
                       param_sig_i_t_dep[[3]],param_sig_i_t_dep[[4]],
                       param_sig_i_t_dep[[5]],param_sig_i_t_dep[[6]],
                       param_sig_i_t_dep[[7]],param_sig_i_t_dep[[8]])


## Alternative: No Weighting Matrix ##
weights = NULL; spec_weights_mat = NULL

## FFT real matrix ##
FFT_matrix = FFT_real_matrix(wave,cosine_indices,n_cos,n)

## Loop Over Different Observation Locations ##
x_locs = seq(0,n-1,1)
y_locs = seq(0,n-1,1)
n_x = length(x_locs); n_y = length(y_locs)
obj_func = matrix(rep(0,n_x*n_y),nrow=n_x,ncol=n_y)
loop_times = rep(0,n_x*n_y)

## Kalman Filter (independent of obs. locations) ##
mu0 = rep(0,K); Sigma0 = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,
                             spec_weights_mat = spec_weights_mat)
G = G_f(param_sig,wave,cosine_indices,dt=1,n_cos=n_cos)
Q = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,spec_weights_mat=spec_weights_mat)

for (i in 1:n_x){
  for (j in 1:n_y){
    
    current_iter = (i-1)*n_x + j
    remaining_iter = n_x*n_y - current_iter
    current_time1 = Sys.time()
    
    last_loc = c(x_locs[i],y_locs[j])
    grid = matrix(data = c(c((n-1)/5,(n-1)/5,(n-1)/5,4*(n-1)/5,4*(n-1)/5,(n-1)/5,4*(n-1)/5,4*(n-1)/5),last_loc),nrow=5,byrow=T)
    coords = lapply(seq_len(nrow(grid)), function(i) grid[i,])
    coords_indices = sapply(coords,function(x) x[1]+n*x[2]+1)
    n_obs = length(coords)
    obs_loc = lapply(1:n_obs,function(j) matrix(coords[[j]],nrow=t,ncol=2,byrow=T))
    
    ## Observation parameters ##
    m = 1; m_indices = list(coords_indices); tau2_vals = list(0.01)
    param_obs = param_obs_func(m,m_indices,tau2_vals)
    param_obs_t_dep =  param_t_dep_func(list(1),list(param_obs))
    
    ## Observation bias ##
    m1 = 0; m1_indices = list(coords_indices); bias_vals = list(0)
    param_bias = param_bias_func(m1,m1_indices,bias_vals)
    param_bias_t_dep = param_t_dep_func(list(1),list(param_bias))
    
    ## Simulate ##
    spde_sim = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                               param_bias_t_dep,seed = 1,n_obs = n_obs,
                               obs_loc = obs_loc, weights = weights)
    y_sim = spde_sim$y
    
    ## Kalman Filter ##
    Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs,coords)
    XBeta = XBeta_f(param_bias,n,t,n_obs)
    R = R_f(param_obs,n,n_obs)
    kf = KalmanFilter(y = y_sim, Phi = Phi, R = R, G = G, Q = Q, t = t, K = K, 
                      log_lik = TRUE, m0 = mu0, R0 = Sigma0)
    
    ## Optimal State Estimate ##
    alpha_hat = kf$m_i_i[t+1,]
    #x_hat = FFT_real_matrix(wave,cosine_indices,n_cos,n)%*%alpha_hat
    
    ## Objective Function ##
    sigma_hat = kf$R_i_i[t+1,,]
    #sigma_hat = objective_function(t,kf$R_i_i,diag(1,nrow=K))
    #sigma_hat_true = FFT_matrix%*%sigma_hat%*%t(FFT_matrix) [not needed]
    obj_func[i,j] = sum(diag(sigma_hat))
    
    current_time2 = Sys.time()
    
    loop_times[current_iter] = as.vector(current_time2-current_time1)
    average_loop_time = mean(loop_times[1:current_iter])
    
    if(current_iter%%50==0){
      cat(remaining_iter, "iterations remaining = ~", round(average_loop_time*remaining_iter/60,3), "mins\n")
    }
    
  }
}

## Save ##
setwd(fig_wd)
saveRDS(obj_func,"param_dep_simulation_8")

## Load ##
if(FALSE){
  obj_func = readRDS("param_dep_simulation_8")
}

## Match obj_func values with (x,y) coords ##
grid = expand.grid(x = 0:(n-1),y = 0:(n-1))/n
grid$obj_func = as.vector(obj_func)

## Plot ##
library("plot3D")
library("extrafont")

## Wireframe Plot ##
if(FALSE){
  setwd(fig_wd)
  #pdf("section4_sim2_plot8.pdf",height=7,width=7)
  par(mar=c(1,.5,.5,2),mfrow=c(1,1))
  #persp3D(z = obj_func,theta = 30,phi=20)
  wireframe(obj_func~x*y, data = grid, cuts = 100, main="", 
            col.regions = jet.col(1000)[1:1000],
            drape = TRUE,zlab="z",xlab="x",ylab="y", light.source = c(20,10,10), 
            aspect = c(1,1), screen = list(z = -45, x = -55), col = F, 
            zoom = 1, scales = list(arrows = F, tck = c(0.8, 0.6, 0.4), 
                                    distance =c(.7, .7, .5)))
  #dev.off()
}

## Image Plot ##
if(TRUE){
  pdf("section4_sim2_plot8.pdf",height=7,width=7)
  par(mar=c(0.5,0.5,0.5,0.5),mgp=c(1.5,1,0),mfrow=c(1,1),family="LM Roman 10")
  image2D(z = obj_func, x=seq(0,1-1/n,1/n), y = seq(0,1-1/n,1/n),main="",xlab="",ylab="",cex.axis=1.5,colkey=F,
          xaxs="i",yaxs="i",resfac=10,xaxt="n",yaxt="n")#colkey = list(cex.axis=2))
  points(x = grid$x[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         y = grid$y[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         col="white",pch=19,cex=4)
  points(x = grid$x[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         y = grid$y[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         col="black",pch=19,cex=2)
  dev.off()
}

##############################
} 

##############################
#### Plot 9 (rho0 = 0.18) ####
##############################

## Parameters ##
param_sig = c(rho0 = 0.18, sigma2 = 0.2, zeta = 0.5,rho1 = 0.1, gamma = 5, 
              alpha = pi/4, mu_x = 0.3, mu_y = -0.3)
param_sig_i_t_dep = lapply(1:length(param_sig), 
                           function(i) param_t_dep_func(1,param_sig[i]))
param_sig_t_dep = 
  param_sig_t_dep_func(param_sig_i_t_dep[[1]],param_sig_i_t_dep[[2]],
                       param_sig_i_t_dep[[3]],param_sig_i_t_dep[[4]],
                       param_sig_i_t_dep[[5]],param_sig_i_t_dep[[6]],
                       param_sig_i_t_dep[[7]],param_sig_i_t_dep[[8]])


## Alternative: No Weighting Matrix ##
weights = NULL; spec_weights_mat = NULL

## FFT real matrix ##
FFT_matrix = FFT_real_matrix(wave,cosine_indices,n_cos,n)

## Loop Over Different Observation Locations ##
x_locs = seq(0,n-1,1)
y_locs = seq(0,n-1,1)
n_x = length(x_locs); n_y = length(y_locs)
obj_func = matrix(rep(0,n_x*n_y),nrow=n_x,ncol=n_y)
loop_times = rep(0,n_x*n_y)

## Kalman Filter (independent of obs. locations) ##
mu0 = rep(0,K); Sigma0 = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,
                             spec_weights_mat = spec_weights_mat)
G = G_f(param_sig,wave,cosine_indices,dt=1,n_cos=n_cos)
Q = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,spec_weights_mat=spec_weights_mat)

for (i in 1:n_x){
  for (j in 1:n_y){
    
    current_iter = (i-1)*n_x + j
    remaining_iter = n_x*n_y - current_iter
    current_time1 = Sys.time()
    
    last_loc = c(x_locs[i],y_locs[j])
    grid = matrix(data = c(c((n-1)/5,(n-1)/5,(n-1)/5,4*(n-1)/5,4*(n-1)/5,(n-1)/5,4*(n-1)/5,4*(n-1)/5),last_loc),nrow=5,byrow=T)
    coords = lapply(seq_len(nrow(grid)), function(i) grid[i,])
    coords_indices = sapply(coords,function(x) x[1]+n*x[2]+1)
    n_obs = length(coords)
    obs_loc = lapply(1:n_obs,function(j) matrix(coords[[j]],nrow=t,ncol=2,byrow=T))
    
    ## Observation parameters ##
    m = 1; m_indices = list(coords_indices); tau2_vals = list(0.01)
    param_obs = param_obs_func(m,m_indices,tau2_vals)
    param_obs_t_dep =  param_t_dep_func(list(1),list(param_obs))
    
    ## Observation bias ##
    m1 = 0; m1_indices = list(coords_indices); bias_vals = list(0)
    param_bias = param_bias_func(m1,m1_indices,bias_vals)
    param_bias_t_dep = param_t_dep_func(list(1),list(param_bias))
    
    ## Simulate ##
    spde_sim = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                               param_bias_t_dep,seed = 1,n_obs = n_obs,
                               obs_loc = obs_loc, weights = weights)
    y_sim = spde_sim$y
    
    ## Kalman Filter ##
    Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs,coords)
    XBeta = XBeta_f(param_bias,n,t,n_obs)
    R = R_f(param_obs,n,n_obs)
    kf = KalmanFilter(y = y_sim, Phi = Phi, R = R, G = G, Q = Q, t = t, K = K, 
                      log_lik = TRUE, m0 = mu0, R0 = Sigma0)
    
    ## Optimal State Estimate ##
    alpha_hat = kf$m_i_i[t+1,]
    #x_hat = FFT_real_matrix(wave,cosine_indices,n_cos,n)%*%alpha_hat
    
    ## Objective Function ##
    sigma_hat = kf$R_i_i[t+1,,]
    #sigma_hat = objective_function(t,kf$R_i_i,diag(1,nrow=K))
    #sigma_hat_true = FFT_matrix%*%sigma_hat%*%t(FFT_matrix) [not needed]
    obj_func[i,j] = sum(diag(sigma_hat))
    
    current_time2 = Sys.time()
    
    loop_times[current_iter] = as.vector(current_time2-current_time1)
    average_loop_time = mean(loop_times[1:current_iter])
    
    if(current_iter%%50==0){
      cat(remaining_iter, "iterations remaining = ~", round(average_loop_time*remaining_iter/60,3), "mins\n")
    }
    
  }
}

## Save ##
setwd(fig_wd)
saveRDS(obj_func,"param_dep_simulation_9")

## Load ##
if(FALSE){
  obj_func = readRDS("param_dep_simulation_9")
}

## Match obj_func values with (x,y) coords ##
grid = expand.grid(x = 0:(n-1),y = 0:(n-1))/n
grid$obj_func = as.vector(obj_func)

## Plot ##
library("plot3D")
library("extrafont")

## Wireframe Plot ##
if(FALSE){
  setwd(fig_wd)
  #pdf("section4_sim2_plot9.pdf",height=7,width=7)
  par(mar=c(1,.5,.5,2),mfrow=c(1,1))
  #persp3D(z = obj_func,theta = 30,phi=20)
  wireframe(obj_func~x*y, data = grid, cuts = 100, main="", 
            col.regions = jet.col(1000)[1:1000],
            drape = TRUE,zlab="z",xlab="x",ylab="y", light.source = c(20,10,10), 
            aspect = c(1,1), screen = list(z = -45, x = -55), col = F, 
            zoom = 1, scales = list(arrows = F, tck = c(0.8, 0.6, 0.4), 
                                    distance =c(.7, .7, .5)))
  #dev.off()
}

## Image Plot ##
if(TRUE){
  pdf("section4_sim2_plot9.pdf",height=7,width=7)
  par(mar=c(0.5,0.5,0.5,0.5),mgp=c(1.5,1,0),mfrow=c(1,1),family="LM Roman 10")
  image2D(z = obj_func, x=seq(0,1-1/n,1/n), y = seq(0,1-1/n,1/n),main="",xlab="",ylab="",cex.axis=1.5,colkey=F,
          xaxs="i",yaxs="i",resfac=10,xaxt="n",yaxt="n")#colkey = list(cex.axis=2))
  points(x = grid$x[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         y = grid$y[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         col="white",pch=19,cex=4)
  points(x = grid$x[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         y = grid$y[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         col="black",pch=19,cex=2)
  dev.off()
}

##############################

##############################
#### Plot 10 (rho0 = 0.20) ####
##############################

## Parameters ##
param_sig = c(rho0 = 0.20, sigma2 = 0.2, zeta = 0.5,rho1 = 0.1, gamma = 5, 
              alpha = pi/4, mu_x = 0.3, mu_y = -0.3)
param_sig_i_t_dep = lapply(1:length(param_sig), 
                           function(i) param_t_dep_func(1,param_sig[i]))
param_sig_t_dep = 
  param_sig_t_dep_func(param_sig_i_t_dep[[1]],param_sig_i_t_dep[[2]],
                       param_sig_i_t_dep[[3]],param_sig_i_t_dep[[4]],
                       param_sig_i_t_dep[[5]],param_sig_i_t_dep[[6]],
                       param_sig_i_t_dep[[7]],param_sig_i_t_dep[[8]])


## Alternative: No Weighting Matrix ##
weights = NULL; spec_weights_mat = NULL

## FFT real matrix ##
FFT_matrix = FFT_real_matrix(wave,cosine_indices,n_cos,n)

## Loop Over Different Observation Locations ##
x_locs = seq(0,n-1,1)
y_locs = seq(0,n-1,1)
n_x = length(x_locs); n_y = length(y_locs)
obj_func = matrix(rep(0,n_x*n_y),nrow=n_x,ncol=n_y)
loop_times = rep(0,n_x*n_y)

## Kalman Filter (independent of obs. locations) ##
mu0 = rep(0,K); Sigma0 = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,
                             spec_weights_mat = spec_weights_mat)
G = G_f(param_sig,wave,cosine_indices,dt=1,n_cos=n_cos)
Q = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,spec_weights_mat=spec_weights_mat)

for (i in 1:n_x){
  for (j in 1:n_y){
    
    current_iter = (i-1)*n_x + j
    remaining_iter = n_x*n_y - current_iter
    current_time1 = Sys.time()
    
    last_loc = c(x_locs[i],y_locs[j])
    grid = matrix(data = c(c((n-1)/5,(n-1)/5,(n-1)/5,4*(n-1)/5,4*(n-1)/5,(n-1)/5,4*(n-1)/5,4*(n-1)/5),last_loc),nrow=5,byrow=T)
    coords = lapply(seq_len(nrow(grid)), function(i) grid[i,])
    coords_indices = sapply(coords,function(x) x[1]+n*x[2]+1)
    n_obs = length(coords)
    obs_loc = lapply(1:n_obs,function(j) matrix(coords[[j]],nrow=t,ncol=2,byrow=T))
    
    ## Observation parameters ##
    m = 1; m_indices = list(coords_indices); tau2_vals = list(0.01)
    param_obs = param_obs_func(m,m_indices,tau2_vals)
    param_obs_t_dep =  param_t_dep_func(list(1),list(param_obs))
    
    ## Observation bias ##
    m1 = 0; m1_indices = list(coords_indices); bias_vals = list(0)
    param_bias = param_bias_func(m1,m1_indices,bias_vals)
    param_bias_t_dep = param_t_dep_func(list(1),list(param_bias))
    
    ## Simulate ##
    spde_sim = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                               param_bias_t_dep,seed = 1,n_obs = n_obs,
                               obs_loc = obs_loc, weights = weights)
    y_sim = spde_sim$y
    
    ## Kalman Filter ##
    Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs,coords)
    XBeta = XBeta_f(param_bias,n,t,n_obs)
    R = R_f(param_obs,n,n_obs)
    kf = KalmanFilter(y = y_sim, Phi = Phi, R = R, G = G, Q = Q, t = t, K = K, 
                      log_lik = TRUE, m0 = mu0, R0 = Sigma0)
    
    ## Optimal State Estimate ##
    alpha_hat = kf$m_i_i[t+1,]
    #x_hat = FFT_real_matrix(wave,cosine_indices,n_cos,n)%*%alpha_hat
    
    ## Objective Function ##
    sigma_hat = kf$R_i_i[t+1,,]
    #sigma_hat = objective_function(t,kf$R_i_i,diag(1,nrow=K))
    #sigma_hat_true = FFT_matrix%*%sigma_hat%*%t(FFT_matrix) [not needed]
    obj_func[i,j] = sum(diag(sigma_hat))
    
    current_time2 = Sys.time()
    
    loop_times[current_iter] = as.vector(current_time2-current_time1)
    average_loop_time = mean(loop_times[1:current_iter])
    
    if(current_iter%%50==0){
      cat(remaining_iter, "iterations remaining = ~", round(average_loop_time*remaining_iter/60,3), "mins\n")
    }
    
  }
}

## Save ##
setwd(fig_wd)
saveRDS(obj_func,"param_dep_simulation_10")

## Load ##
if(FALSE){
  obj_func = readRDS("param_dep_simulation_10")
}

## Match obj_func values with (x,y) coords ##
grid = expand.grid(x = 0:(n-1),y = 0:(n-1))/n
grid$obj_func = as.vector(obj_func)

## Plot ##
library("plot3D")
library("extrafont")

## Wireframe Plot ##
if(FALSE){
  setwd(fig_wd)
  #pdf("section4_sim2_plot10.pdf",height=7,width=7)
  par(mar=c(1,.5,.5,2),mfrow=c(1,1))
  #persp3D(z = obj_func,theta = 30,phi=20)
  wireframe(obj_func~x*y, data = grid, cuts = 100, main="", 
            col.regions = jet.col(1000)[1:1000],
            drape = TRUE,zlab="z",xlab="x",ylab="y", light.source = c(20,10,10), 
            aspect = c(1,1), screen = list(z = -45, x = -55), col = F, 
            zoom = 1, scales = list(arrows = F, tck = c(0.8, 0.6, 0.4), 
                                    distance =c(.7, .7, .5)))
  #dev.off()
}

## Image Plot ##
if(TRUE){
  pdf("section4_sim2_plot10.pdf",height=7,width=7)
  par(mar=c(0.5,0.5,0.5,0.5),mgp=c(1.5,1,0),mfrow=c(1,1),family="LM Roman 10")
  image2D(z = obj_func, x=seq(0,1-1/n,1/n), y = seq(0,1-1/n,1/n),main="",xlab="",ylab="",cex.axis=1.5,colkey=F,
          xaxs="i",yaxs="i",resfac=10,xaxt="n",yaxt="n")#colkey = list(cex.axis=2))
  points(x = grid$x[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         y = grid$y[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         col="white",pch=19,cex=4)
  points(x = grid$x[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         y = grid$y[which(grid$obj_func[1:(n^2/2)]==min(grid$obj_func[1:(n^2/2)]))],
         col="black",pch=19,cex=2)
  dev.off()
}

##############################

##################
#### RML_ROSP ####
##################

## Setup ##
if(TRUE){
  
  ## Prelims ##
  if(TRUE){
    n = 20; t = 20000; dt = 1; K = 21; nu = 1
  }
  
  ## RML & OSP Times ##
  if(TRUE){
    t_RML = 2; t_ROSP = 38
  }
  
  ## Initial Observation Locations ##
  if(TRUE){
    
    ## Coordinates ##
    obs0 = list(c(n/2+0.1,n/6+0.1),
                c((n-1)/5,(n-1)/5),
                c((n-1)/5,4*(n-1)/5),
                c(4*(n-1)/5,(n-1)/5),
                c(4*(n-1)/5,4*(n-1)/5))
    
    fixed_obs = list(c((n-1)/5,(n-1)/5),
                     c((n-1)/5,4*(n-1)/5),
                     c(4*(n-1)/5,(n-1)/5),
                     c(4*(n-1)/5,4*(n-1)/5))
    
    ## Indices ##
    obs0_indices = sapply(obs0,function(x) x[1]+n*x[2]+1)
    
    ## Number of observations ##
    n_obs = length(obs0)
    
    ## Observation location at all times ##
    obs_loc = lapply(1:n_obs,
                     function(j) matrix(obs0[[j]],nrow=t,ncol=2,byrow=T))
    
    ## Plot initial observation locations ##
    plot_obs = TRUE
    if(plot_obs){
      par(mfrow=c(1,1))
      plot(sapply(obs0,function(x) x)[1,],sapply(obs0,function(x) x)[2,],
           xlim = c(0,n-1),ylim=c(0,n-1),xlab=expression(x),ylab=expression(y),
           main="Sensor Locations")
      abline(v=0:(n-1),h=(0:n-1),col="gray92",lty=2)
    }
    
    ## 
  }
  
  ## True Parameters ##
  if(TRUE){
    
    ## Signal parameters ##
    param_sig = c(rho0 = 0.3, sigma2 = 0.2, zeta = 0.5,rho1 = 0.1, gamma = 5, 
                  alpha = pi/4, mu_x = 0.3, mu_y = -0.3)
    
    param_sig_i_t_dep = lapply(1:length(param_sig), 
                               function(i) param_t_dep_func(1,param_sig[i]))
    
    param_sig_t_dep = 
      param_sig_t_dep_func(param_sig_i_t_dep[[1]],param_sig_i_t_dep[[2]],
                           param_sig_i_t_dep[[3]],param_sig_i_t_dep[[4]],
                           param_sig_i_t_dep[[5]],param_sig_i_t_dep[[6]],
                           param_sig_i_t_dep[[7]],param_sig_i_t_dep[[8]])
    
    ## Observation parameters ##
    m = 1; m_indices = list(obs0_indices); tau2_vals = list(0.01)
    param_obs = param_obs_func(m,m_indices,tau2_vals)
    param_obs_t_dep =  param_t_dep_func(list(1),list(param_obs))
    
    ## Observation bias ##
    m1 = 0; m1_indices = list(obs0_indices); bias_vals = list(0)
    param_bias = param_bias_func(m1,m1_indices,bias_vals)
    param_bias_t_dep = param_t_dep_func(list(1),list(param_bias))
    
  }
  
  ## Number of parameters ##
  if(TRUE){
    n_param = length(param_sig)+param_obs$m+param_bias$m
  }
  
  ## Parameter names ##
  if(TRUE){
    param_sig_names = list(parse(text="rho[0]"),parse(text="sigma^2"),
                           parse(text="zeta"),parse(text="rho[1]"),
                           parse(text="gamma"),parse(text="alpha"),
                           parse(text="mu[x]"),parse(text="mu[y]"))
    param_obs_names = list(parse(text="tau^2"))
    param_bias_names = list(parse(text="beta"))
  }
  
  ## Initial parameters ##
  if(TRUE){
    param0 = c(rho0 = 0.01, sigma2 = 0.2, zeta = 0.5,rho1 = 0.1, gamma = 5, 
               alpha = pi/4, mu_x = 0.3, mu_y = -0.3, tau2 = 0.01)
    param_sig0 = param0[1:8]
    param_obs0 = param_obs_func(m,m_indices,list(param0[9]))
    param_bias0 = param_bias_func(m1,m1_indices,list(0))
  }
  
  ## Signal Noise Weights ##
  if(TRUE){
    spec_weights_mat = NULL
  }
  
  ## SPDE_FT object ##
  if(TRUE){
    SPDE_FT = spde_initialise(n,t,K)
    wave = SPDE_FT$wave
    cosine_indices = SPDE_FT$cosine_indices
    n_cos = SPDE_FT$n_cos
    n_cos_sin = length(cosine_indices)
    K = SPDE_FT$K
  }
  
  ## SPDE_SIM ##
  if(FALSE){
    spde_sim = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                               param_bias_t_dep,seed = 1,n_obs = n_obs,
                               obs_loc = obs_loc, alpha_sim = TRUE,
                               x_sim = TRUE, obs_noise_sim = TRUE,
                               obs_bias_sim = TRUE, y_sim = FALSE)
  }
  
  ## Log scale & log indices ##
  if(TRUE){
    log_scale = FALSE
    log_indices = c(1,2,3,4,5,8+(1:param_obs0$m))
  }
  
  ## Parameters to estimate ##
  if(TRUE){
    grad_params = 1
  }
  
  ## Observations to be optimised ##
  if(TRUE){
    grad_obs = 1
  }
  
  ## Hessian indices ##
  if(TRUE){
    hessian_indices = NULL
  }
  
  ## Number of iterations ##
  if(TRUE){
    n_iterations = 1
  }
  
  ## RML Step Sizes ##
  if(TRUE){
    
    ## Initial step-size vector ##
    step_size_vec = c(0.1,rep(1,8))
    
    ## Step size matrix ##
    step_size_RML = matrix(0,nrow=n_iterations*t,ncol=n_param)
    for (i in 1:dim(step_size_RML)[1]){
      step_size_RML[i,] = step_size_vec*i^(-0.55)
    }
    
  }
  
  ## ROSP Step sizes ##
  if(TRUE){
    
    ## Initial step size vector ##
    step_size_vec = c(0.1,0.1)
    
    ## Step size matrix ##
    step_size_mat = matrix(0,nrow=n_iterations*t,ncol=2)
    for (i in 1:(n_iterations*t)){
      step_size_mat[i,]=step_size_vec*i^(-0.51)
    }
    
    ## Step size matrices for all sensors ##
    step_size_ROSP = rep(list(step_size_mat),n_obs)
  }
  
  ## Weighting matrix ##
  if(TRUE){
    
    ## Weighting matrix (physical space) ##
    W = diag(1,n^2)
    W_coords = NULL
    
    ## Weighting matrix (spectral space) ##
    coords2 = lapply(1:n^2,function(i) as.numeric(expand.grid(0:(n-1),0:(n-1))[i,]))
    Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs=n^2,coords = coords2)
    W_fourier = t(Phi)%*%W%*%Phi
  }
  
  ## RML Plot? ##
  if(TRUE){
    RML_plot = TRUE
  }
  
  ## RML Plot Limits ##
  if(TRUE){
    plot_limits = cbind(0.5*sapply(1:8, function(i) min(abs(param_sig0[[i]]),abs(param_sig[[i]]))*abs(param_sig[[i]])/param_sig[[i]]),
                        1.5*sapply(1:8, function(i) max(abs(param_sig0[[i]]),abs(param_sig[[i]]))*abs(param_sig[[i]])/param_sig[[i]]))
    plot_limits = rbind(plot_limits,
                        cbind(0,1.5*max(unlist(param_obs$tau2_vals),unlist(param_obs0$tau2_vals))))    
    plot_limits = rbind(plot_limits,
                        cbind(0,1.5*max(unlist(param_bias$bias_vals),unlist(param_bias0$bias_vals)))) 
    plot_limits = t(apply(plot_limits,1,sort)) 
    plot_limits = lapply(1:dim(plot_limits)[1],function(i) plot_limits[i,])
  }
  
  ## Compute mean parameter estimate? ##
  if(TRUE){
    mean_param = FALSE
  }
  
  ## Parameter estimation period ##
  if(TRUE){
    param_est_starts = list(3000,2000,3000,3000,4000,2000,2000,2000,1000)
    param_est_ends = rep(list(n_iterations*t),n_param)
  }
  
  ## RML Plot Point Frequency ##
  if(TRUE){
    RML_plot_point_freq = 10
  }
  
  ## RML Plot Frequency ##
  if(TRUE){
    RML_plot_freq = 1000
  }
  
  ## Save RML plot? ##
  if(TRUE){
    save_RML_plot = TRUE
  }
  
  ## RML Plot Filename ##
  if(TRUE){
    RML_plot_filename = "RML_plot1.pdf"
  }
  
  ## ROSP Plots? ##
  if(TRUE){
    ROSP_plot = FALSE
    ROSP_plot2d = TRUE
  }
  
  ## ROSP Plot Point Frequency ##
  if(TRUE){
    ROSP_plot_point_freq = 10
  }
  
  ## ROSP Plot Frequency ##
  if(TRUE){
    ROSP_plot_freq = 2000
  }
  
  ## Save ROSP Plots? ##
  if(TRUE){
    save_ROSP_plot = FALSE
    save_ROSP_plot2d = TRUE
  }
  
  ## ROSP Plot Filenames ##
  if(TRUE){
    ROSP_plot_filename = "ROSP_plot1.pdf"
    ROSP_plot2d_filename = "ROSP_plot2d1.pdf"
  }
  
  ## Plot ROSP grid? ##
  if(TRUE){
    ROSP_grid = TRUE
  }
  
  ## Save RML? ##
  if(TRUE){
    save_RML = FALSE
  }
  
  ## RML Filename ##
  if(TRUE){
    RML_filename = "RML"
  }
  
  ## Save ROSP? ##
  if(TRUE){
    save_ROSP = FALSE
  }
  
  ## ROSP Filename ##
  if(TRUE){
    ROSP_filename = "ROSP"
  }
  
  ## RML Print Iteration Frequency ##
  if(TRUE){
    RML_print_iter_freq = 100
  }
  
  ## ROSP Print Iteration Frequency ##
  if(TRUE){
    ROSP_print_iter_freq = 100
  }
  
  ## Print Iteration Frequency ##
  if(TRUE){
    print_iter_freq = 1
  }
  
  ## RML Back Track Line Search? ##
  if(TRUE){
    RML_back_track = FALSE
  }
  
  ## RML Back Track Line Search Tolerance ##
  if(TRUE){
    RML_back_track_tol = 1e-3
  }
  
  ## RML Back Track Line Search for Each Param? ##
  if(TRUE){
    RML_back_track_independent = TRUE
    RML_back_track_param = 1:n_param
  }
  
  ## RML Back Track Line Search Time Independent? ##
  if(TRUE){
    RML_back_track_time_independent = FALSE
  }
  
  ## RML Back Track Scale ##
  if(TRUE){
    RML_back_track_scale = 0.5
  }
  
  ## ROSP Back Track Line Search? ##
  if(TRUE){
    ROSP_back_track = FALSE
  }
  
  ## ROSP Back Track Line Search Tolerance ##
  if(TRUE){
    ROSP_back_track_tol = 1e-7
  }
  
  ## ROSP Back Track Line Search for Each Param? ##
  if(TRUE){
    ROSP_back_track_independent = TRUE
    ROSP_back_track_obs = 1:n_obs
  }
  
  ## ROSP Back Track Line Search Time Independent? ##
  if(TRUE){
    ROSP_back_track_time_independent = TRUE
  }
  
  ## ROSP Back Track Scale ##
  if(TRUE){
    ROSP_back_track_scale = 0.5
  }
  
  ## Highlight ROSP initial points ##
  if(TRUE){
    ROSP_plot_initial = TRUE
  }
  
  ## Add ROSP legend? ##
  if(TRUE){
    ROSP_leg = TRUE
  }
  
}

## RML_ROSP ##
if(TRUE){
  RML_ROSP_sim2 = RML_ROSP(spde_sim = spde_sim, n = n, n_obs = n_obs, t = t, 
                            t_RML = t_RML, t_ROSP = t_ROSP, K = K, 
                            n_param = n_param, param_obs0 = param_obs0,
                            param_sig0 = param_sig0, param_bias0 = param_bias0,
                            log_scale = log_scale, log_indices = log_indices, 
                            grad_params = grad_params, 
                            hessian_indices = hessian_indices, obs0 = obs0, 
                            grad_obs = grad_obs, n_iterations = n_iterations, 
                            step_size_RML = step_size_RML, 
                            step_size_ROSP = step_size_ROSP, 
                            W_fourier = W_fourier, 
                            RML_plot = RML_plot, plot_limits = plot_limits, 
                            param_obs_true = param_obs_t_dep, 
                            param_sig_true = param_sig_t_dep, 
                            param_bias_true = param_bias_t_dep, 
                            param_obs_names = param_obs_names, 
                            param_sig_names = param_sig_names, 
                            param_bias_names = param_bias_names, 
                            mean_param = mean_param, 
                            param_est_starts = param_est_starts, 
                            param_est_ends = param_est_ends, 
                            RML_plot_point_freq = RML_plot_point_freq, 
                            RML_plot_freq = RML_plot_freq, 
                            save_RML_plot = save_RML_plot, 
                            RML_plot_filename = RML_plot_filename,
                            ROSP_plot = ROSP_plot, 
                            ROSP_plot2d = ROSP_plot2d, 
                            W_coords = W_coords, 
                            ROSP_plot_point_freq = ROSP_plot_point_freq, 
                            ROSP_plot_freq = ROSP_plot_freq,
                            save_ROSP_plot = save_ROSP_plot, 
                            save_ROSP_plot2d = save_ROSP_plot2d, 
                            ROSP_plot_filename = ROSP_plot_filename,
                            ROSP_plot2d_filename = ROSP_plot2d_filename,
                            ROSP_grid = ROSP_grid, save_RML = save_RML, 
                            RML_filename = RML_filename,
                            save_ROSP = save_ROSP, ROSP_filename = ROSP_filename,
                            RML_print_iter_freq = RML_print_iter_freq,
                            ROSP_print_iter_freq = ROSP_print_iter_freq,
                            print_iter_freq = print_iter_freq,
                            RML_back_track = RML_back_track,
                            RML_back_track_tol = RML_back_track_tol, 
                            RML_back_track_independent = RML_back_track_independent,
                            RML_back_track_time_independent = RML_back_track_time_independent, 
                            RML_back_track_scale = RML_back_track_scale,
                            RML_back_track_param = RML_back_track_param,
                            ROSP_back_track = ROSP_back_track,
                            ROSP_back_track_tol = ROSP_back_track_tol, 
                            ROSP_back_track_independent = ROSP_back_track_independent,
                            ROSP_back_track_time_independent = ROSP_back_track_time_independent, 
                            ROSP_back_track_scale = ROSP_back_track_scale,
                            ROSP_back_track_obs = ROSP_back_track_obs,
                            ROSP_plot_initial = ROSP_plot_initial,
                            ROSP_leg = ROSP_leg)
}

## Save ##
setwd(fig_wd)
saveRDS(RML_ROSP_sim2,"RML_ROSP_sim2")

## Load ##
if(FALSE){
  RML_ROSP_sim2 = readRDS("RML_ROSP_sim2")
}

##################

##############
#### PLOT ####
##############

## RML ##
if(TRUE){
  setwd(fig_wd)
  #pdf("RML_sim2_plot.pdf",width=7,height=7)
  tikz('RML_sim2_plot.tex', standAlone=TRUE, height = 7, width = 14)
  #par(mar=c(4,4.5,.6,.6),mgp=c(2.7,1,0),family="LM Roman 10",cex=2)
  par(mfrow=c(1,2))
  par(mar = c(4,4.5,6,3),mgp = c(2.7,1,0),xpd=F)
  
  param_est = RML_ROSP_sim2$param_est
  
  t_0 = 1
  t_1 = min(which(param_est>0.03))
  t_2 = min(which(param_est>0.10))
  t_3 = min(which(param_est>0.15))
  t_4 = min(which(param_est>0.20))
  t_5 = t
  
  plot_points = list(seq(t_0,t_1,5),seq(t_1+1,t_2,5),seq(t_2+1,t_3,5),
                     seq(t_3+1,t_4,5),seq(t_4+1,t_5,5))
  
  param_sig_indices = 1
  param_sig_true = param_sig_t_dep
  param_sig_names = list("${\\rho}_0(t)$")
  
  ## Signal Parameters ##
  for (k in param_sig_indices){
    
    plot(plot_points[[1]], param_est[plot_points[[1]],k], cex=0.8,
         xlim = c(1,n_iterations*t+1), xlab = "t",ylab = param_sig_names[[k]],
         ylim = c(0,0.32), main = "", cex.lab = 1.8, cex.axis = 1.8,
         pch = 19, col = jet.col(length(plot_points))[1])
    
    for (l in 2:length(plot_points)){
      points(plot_points[[l]], param_est[plot_points[[l]],k], cex = 0.8, 
             pch = 19,col=jet.col(length(plot_points))[l])
    }
    
    abline(h=param_sig_true$param_vals[[1]][k],col="red",lty=2)
    
    #legend(x=10000,y=0.2, legend = c("$0.00\\leq\\rho_0(t)<0.03$","$0.03\\leq\\rho_0(t)<0.10$","$0.10\\leq\\rho_0(t)<0.15$","$0.15\\leq\\rho_0(t)<0.20$","$0.20\\leq\\rho_0(t)<0.30$"),
    #       col = jet.col(length(plot_points))[1:length(plot_points)],pch = 19, cex = 1.8)
  }
  
  #dev.off()
  #tools::texi2dvi('RML_sim2_plot.tex',pdf=T)
#}

## ROSP ##
#if(TRUE){
  setwd(fig_wd)
  #pdf("ROSP_sim2_plot.pdf",width=7,height=7)
  #tikz('ROSP_sim2_plot.tex', standAlone=TRUE, height = 7, width = 7)
  #par(mar=c(4,4.5,.6,.6),mgp=c(2.7,1,0),family="LM Roman 10")
  #par(mfrow=c(1,1))
  par(mar = c(4,4.9,6,0.6),mgp = c(2.7,1,0),xpd=F)
  
  obs = RML_ROSP_sim2$obs
  param_est = RML_ROSP_sim2$param_est
  
  t_0 = 1
  t_1 = min(which(param_est>0.03))
  t_2 = min(which(param_est>0.10))
  t_3 = min(which(param_est>0.15))
  t_4 = min(which(param_est>0.20))
  t_5 = t
  
  t_freq_0 = 10
  t_freq_1 = 200
  t_freq_2 = 200
  t_freq_3 = 100
  t_freq_4 = 500
  
  plot_points = list(seq(t_0,t_1,t_freq_0),
                     seq(t_1+1,t_2,t_freq_1),
                     seq(t_2+1,t_3,t_freq_2),
                     seq(t_3+1,t_4,t_freq_3),
                     seq(t_4+1,t_5,t_freq_4))
  
  ## Initial Point ##
  plot(obs[[1]][1,1],obs[[1]][1,2],xlim = c(0,n),ylim = c(0,n),
       cex = 4, pch = 21, col=brewer.pal(8,"Greens")[7],
       bg = brewer.pal(8,"Greens")[2], xlab="$x$",
       ylab="$y$",cex.lab = 1.8,cex.axis = 1.8,
       xaxt ="n", yaxt = "n", xaxs = "i", yaxs = "i")
  
  ## Fixed Points ##
  for(k in 1:length(fixed_obs)){
    points(fixed_obs[[k]][1],fixed_obs[[k]][2],cex=4,
           pch=21,col=brewer.pal(8,"Reds")[7],
           bg = brewer.pal(8,"Reds")[2])
  }
  
  ## Final Point ##
  for (k in grad_obs){
    points(obs[[k]][n_iterations*t,1],obs[[k]][n_iterations*t,2],
           cex = 4, pch = 21, col = brewer.pal(8,"Blues")[7],
           bg = brewer.pal(8,"Blues")[2])
  }
  
  ## Plot remaining points ##
  for (l in 1:length(plot_points)){
    points(obs[[1]][plot_points[[l]],1],obs[[1]][plot_points[[l]],2],
           col=jet.col(length(plot_points))[l], cex=1.2, pch = 19)
  }
  
  ## Axes ##
  axis(side=1,at=seq(0,n,length.out=6),labels = c(0,0.2,0.4,0.6,0.8,1.0),
       cex.lab = 1.8,cex.axis = 1.8)
  axis(side=2,at=seq(0,n,length.out=6),labels = c(0,0.2,0.4,0.6,0.8,1.0),
       cex.lab = 1.8, cex.axis = 1.8)
  
  ## Grid Lines ##
  abline(h=seq(1,(n-1),2),v=seq(1,(n-1),2),col="gray90",lty=2)
  
  ## Legend ##
  legend(x = 5.0, y = 19.5, pt.cex = c(4,4,4,4),
         pch = c(21,21,21),
         legend=c("Initial Sensor Locations",
                  "Fixed Sensor Locations",
                  "Final Sensor Locations"),
         pt.bg = c(brewer.pal(8,"Greens")[2],
                   brewer.pal(8,"Reds")[2],
                   brewer.pal(8,"Blues")[2]),
         col = c(brewer.pal(8,"Greens")[7],
                 brewer.pal(8,"Reds")[7],
                 brewer.pal(8,"Blues")[7]),
         cex = 1.6, y.intersp = 1.2)
  
  par(xpd=NA)
  legend(x = -23.1, y = 23, legend = c("$0.00\\leq\\rho_0<0.03$","$0.03\\leq\\rho_0<0.10$","$0.10\\leq\\rho_0<0.15$","$0.15\\leq\\rho_0<0.20$","$0.20\\leq\\rho_0<0.30$"),
         col = jet.col(length(plot_points))[1:length(plot_points)],pch = 19, cex = 1.8, ncol=5,
         x.intersp = 0.6, text.width = c(7.35,7.35,7.35,7.35,7.35))
  
  dev.off()
  tools::texi2dvi('RML_sim2_plot.tex',pdf=T)
}

##############


