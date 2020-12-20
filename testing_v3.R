## spde_v3.R Testing ##

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

#############
#### BML ####
#############

## Setup ##
if(TRUE){
  
  ## Prelims ##
  if(TRUE){
    n = 6; t = 100; dt = 0.1; K = 7; nu = 1
  }
  
  ## Observation Locations ##
  if(TRUE){
    
    ## Observation locations 
    coords = list(c(1.4,1.6),c(4.3,4.7),c(1.3,4.2),c(4.7,1.9),c(3.3,3.2),
                  c(0.1,2.0),c(5.0,4.8),c(0.1,0.2),c(2.8,5.4),c(3.6,0.2))
    
    coords = lapply(1:n^2,
                    function(i) as.numeric(expand.grid(0:(n-1),0:(n-1))[i,]))
    
    ## Observation locations 'indices'
    coords_indices = sapply(coords,function(x) x[1]+n*x[2]+1)
    
    ## Number of observations 
    n_obs = length(coords)
    
    ## Observation location at all times
    obs_loc = lapply(1:n_obs,function(j) matrix(coords[[j]],nrow=t,ncol=2,byrow=T))
    
    ## Plot observation locations
    plot_obs = TRUE
    if(plot_obs){
      par(mfrow=c(1,1))
      plot(sapply(coords,function(x) x)[1,],sapply(coords,function(x) x)[2,],
           xlim = c(0,n-1),ylim=c(0,n-1),xlab=expression(x),ylab=expression(y),
           main="Sensor Locations")
      abline(v=0:(n-1),h=(0:n-1),col="gray92",lty=2)
    }
  }
  
  ## Parameters ##
  if(TRUE){
    
    ## Signal parameters ##
    param_sig = c(rho0 = 0.5, sigma2 = 0.2, zeta = 0.5,rho1 = 0.1, gamma = 2, 
                  alpha = pi/4, mu_x = 0.3, mu_y = -0.3)
    param_sig_i_t_dep = lapply(1:length(param_sig), 
                               function(i) param_t_dep_func(1,param_sig[i]))
    
    param_sig_t_dep = 
      param_sig_t_dep_func(param_sig_i_t_dep[[1]],param_sig_i_t_dep[[2]],
                           param_sig_i_t_dep[[3]],param_sig_i_t_dep[[4]],
                           param_sig_i_t_dep[[5]],param_sig_i_t_dep[[6]],
                           param_sig_i_t_dep[[7]],param_sig_i_t_dep[[8]])
    
    ## Observation parameters ##
    m = 1; m_indices = list(coords_indices); tau2_vals = list(0.01)
    param_obs = param_obs_func(m,m_indices,tau2_vals)
    param_obs_t_dep =  param_t_dep_func(list(1),list(param_obs))
    
    ## Observation bias ##
    m1 = 0; m1_indices = list(coords_indices); bias_vals = list(0)
    param_bias = param_bias_func(m1,m1_indices,bias_vals)
    param_bias_t_dep = param_t_dep_func(list(1),list(param_bias))
    
    ## Parameter names ##
    param_sig_names = list(parse(text="rho[0]"),parse(text="sigma^2"),
                           parse(text="zeta"),parse(text="rho[1]"),
                           parse(text="gamma"),parse(text="alpha"),
                           parse(text="mu[x]"),parse(text="mu[y]"))
    param_obs_names = list(parse(text="tau^2"))
    param_bias_names = list(parse(text="beta"))
    
    ## Number of parameters ##
    n_param = length(param_sig)+param_obs$m+param_bias$m
  }
  
  ## Signal Noise Weights ##
  if(TRUE){
    
    ## Weighting Matrix (Physical Space) ##
    grid = expand.grid(x = 0:(n-1),y = 0:(n-1))/n
    centre_x = 6/12; centre_y = 7/12; sigma_x = 0.3; sigma_y = 0.3; rho = 0
    weights = with(grid,sech2d(x = x,y = y, centre_x = centre_x,
                               centre_y = centre_y, sigma_x = sigma_x,
                               sigma_y=sigma_y, rho = rho))
    weights_mat = weights_mat_func(weights)
    
    ## Weighting Matrix (Spectral Space) ##
    spec_weights_mat = spec_weights_mat_func(n = n, K = K, 
                                             weights_mat = weights_mat)
    
    ## Alternative: No Weighting Matrix ##
    weights = NULL; spec_weights_mat = NULL
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
  
  ## Simulate Y_t ##
  if(TRUE){
    spde_sim = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                               param_bias_t_dep,seed = 1,n_obs = n_obs,
                               obs_loc = obs_loc, weights = weights, dt = dt)
    y_sim = spde_sim$y
  }
  
  ## Weighting matrix ##
  if(TRUE){
    W_coords = NULL; W = diag(1,n^2)
    coords2 = lapply(1:n^2,function(i) as.numeric(expand.grid(0:(n-1),0:(n-1))[i,]))
    Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs=n^2,coords = coords2)
    W_fourier = t(Phi)%*%W%*%Phi
  }
  
  ## Initial parameters ##
  if(TRUE){
    param0 = c(rho0 = 0.25, sigma2 = 0.5, zeta = 0.2,rho1 = 0.2, gamma = 1.5, 
               alpha = pi/3, mu_x = 0.1, mu_y = -0.2, tau2 = 0.05)
    param_sig0 = param0[1:8]
    param_obs0 = param_obs_func(m,m_indices,list(param0[9]))
    param_bias0 = param_bias_func(m1,m1_indices,list(0))
  }
  
  ## Log indices ##
  if(TRUE){
    log_indices = c(1,2,3,4,5,8+(1:param_obs0$m))
  }
  
  ## Parameter plot limits ##
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
  
  ## Parameters to estimate ##
  if(TRUE){
    grad_params = 1:n_param
  }
  
  ## Hessian indices ##
  if(TRUE){
    hessian_indices = NULL
  }
  
  ## Number of iterations ##
  if(TRUE){
    n_iterations = 1000
  }
  
  ## Step sizes ##
  if(TRUE){
    
    ## Initial step-size vector ##
    step_size_vec = 0.005*c(0.01,0.1,0.1,0.005,0.05,0.01,0.01,0.01,0.001)
    
    ## Step size matrix ##
    step_size = matrix(0,nrow=n_iterations,ncol=n_param)
    for (i in 1:dim(step_size)[1]){
      step_size[i,] = step_size_vec 
    }
  }
  
  ## Parameter estimation period ##
  if(TRUE){
    param_est_starts = list(3000,2000,3000,3000,4000,2000,2000,2000,1000)
    param_est_ends = rep(list(n_iterations),n_param)
  }
  
  ## Plot vertical (burn in) lines? ##
  if(TRUE){
    plot_vertical = FALSE
  }
  
  ## Plot Point Frequency ##
  if(TRUE){
    plot_point_freq = 5
  }
  
  ## Plot Frequency ##
  if(TRUE){
    plot_freq = 50
  }
  
  ## Filename ##
  if(TRUE){
    filename = paste("n",n,"_t",t,"_nobs",n_obs,sep="")
    plot_filename = paste(filename,".pdf",sep="")
  }
  
  ## Save ##
  if(TRUE){
    save = FALSE 
    save_plot = FALSE
  }
  
  ## Print Iteration Frequency ##
  if(TRUE){
    print_iter_freq = 10
  }
  
  ## Back Track Line Search? ##
  if(TRUE){
    back_track = FALSE
  }
  
  ## Back Track Line Search Tolerance ##
  if(TRUE){
    back_track_tol = 1e-10
  }
  
  ## Back Track Line Search for Each Param? ##
  if(TRUE){
    back_track_independent = TRUE
    back_track_param = 1:n_param
  }
  
  ## Back Track Line Search Time Independent? ##
  if(TRUE){
    back_track_time_independent = FALSE
  }
  
  ## Back Track Scale ##
  if(TRUE){
    back_track_scale = 0.2
  }
}

## BML ##
if(TRUE){
  if(n_obs == n^2){
    BML_test2 = BML_spectral(w = y_sim, n = n, t = t, K = K, 
                            param_obs0 = param_obs0, param_sig0 = param_sig0, 
                            n_param = n_param,grad_params = grad_params, 
                            step_size = step_size,log_scale = FALSE, 
                            log_indices = log_indices,
                            hessian_indices = hessian_indices, 
                            n_iterations = n_iterations,
                            plot = TRUE, plot_limits = plot_limits,
                            param_obs_true = param_obs_t_dep,
                            param_sig_true = param_sig_t_dep, 
                            param_obs_names = param_obs_names,
                            param_sig_names = param_sig_names, 
                            param_est_starts = param_est_starts,
                            param_est_ends = param_est_ends,
                            plot_point_freq = plot_point_freq, 
                            plot_freq = plot_freq,save_plot = save_plot, 
                            plot_filename = plot_filename, save = save, 
                            filename = filename, mean_param = TRUE,
                            print_iter_freq = print_iter_freq,
                            back_track = back_track, 
                            back_track_tol = back_track_tol,
                            back_track_independent = back_track_independent,
                            back_track_time_independent = back_track_time_independent,
                            back_track_scale = back_track_scale,
                            back_track_param = back_track_param,
                            plot_vertical = plot_vertical,
                            W_fourier = W_fourier, dt = dt)
  }
  if(n_obs < n^2){
    BML_test = BML(y = y_sim, n = n, t = t, K = K, param_obs0 = param_obs0, 
                   param_sig0 = param_sig0, param_bias0 = param_bias0, 
                   n_param = n_param, grad_params = grad_params, 
                   step_size = step_size, log_scale = FALSE, 
                   log_indices = log_indices, hessian_indices = hessian_indices, 
                   n_obs = n_obs, coords = coords, n_iterations = n_iterations,
                   spec_weights_mat = spec_weights_mat, plot = TRUE, 
                   plot_limits = plot_limits,
                   param_obs_true = param_obs_t_dep,
                   param_sig_true = param_sig_t_dep, 
                   param_bias_true = param_bias_t_dep,
                   param_obs_names = param_obs_names,
                   param_sig_names = param_sig_names, 
                   param_bias_names = param_bias_names,
                   param_est_starts = param_est_starts,
                   param_est_ends = param_est_ends,
                   plot_point_freq = plot_point_freq, plot_freq = plot_freq,
                   save_plot = save_plot, plot_filename = plot_filename, 
                   save = save, filename = filename, mean_param = TRUE,
                   print_iter_freq = print_iter_freq,
                   back_track = back_track, back_track_tol = back_track_tol,
                   back_track_independent = back_track_independent,
                   back_track_time_independent = back_track_time_independent,
                   back_track_scale = back_track_scale,
                   back_track_param = back_track_param,
                   plot_vertical = plot_vertical,
                   W_fourier = W_fourier, dt = dt)
  }
}

#############

#############
#### RML ####
#############

## Setup ##
if(TRUE){
  
  ## Prelims ##
  if(TRUE){
    n = 6; t = 10000; dt = 1; K = 10; nu = 1
  }
  
  ## Observation Locations ##
  if(TRUE){
    
    all_obs = FALSE
    
    ## Observation locations 
    if(!all_obs){
      coords = list(c(1.4,1.6),c(4.3,4.7),c(1.3,4.2),c(4.7,1.9),c(3.3,3.2),
                    c(0.1,2.0),c(5.0,4.8),c(0.1,0.2),c(2.8,5.4),c(3.6,0.2),
                    c(2.5,1.6),c(2.5,3.0),c(2.8,4.2),c(3.02,2.0))
    }
    
    if(all_obs){
      coords = lapply(1:n^2,
                      function(i) as.numeric(expand.grid(0:(n-1),0:(n-1))[i,]))
    }
    
    ## Observation locations 'indices'
    coords_indices = sapply(coords,function(x) x[1]+n*x[2]+1)
    
    ## Number of observations 
    n_obs = length(coords)
    
    ## Observation location at all times
    obs_loc = lapply(1:n_obs,function(j) matrix(coords[[j]],nrow=t,ncol=2,byrow=T))
    
    ## Plot observation locations
    plot_obs = FALSE
    if(plot_obs){
      par(mfrow=c(1,1))
      plot(sapply(coords,function(x) x)[1,],sapply(coords,function(x) x)[2,],
           xlim = c(0,n-1),ylim=c(0,n-1),xlab=expression(x),ylab=expression(y),
           main="Sensor Locations")
      abline(v=0:(n-1),h=(0:n-1),col="gray92",lty=2)
    }
  }
  
  ## Parameters ##
  if(TRUE){
      
    ## Signal parameters ##
    param_sig = c(rho0 = 0.5, sigma2 = 0.2, zeta = 0.5,rho1 = 0.1, gamma = 2, 
                  alpha = pi/4, mu_x = 0.3, mu_y = -0.3)
    
    time_dep_sig_param = FALSE
    
    if(!time_dep_sig_param){
      param_sig_i_t_dep = lapply(1:length(param_sig), 
                                 function(i) param_t_dep_func(1,param_sig[i]))
    
    param_sig_t_dep = 
      param_sig_t_dep_func(param_sig_i_t_dep[[1]],param_sig_i_t_dep[[2]],
                           param_sig_i_t_dep[[3]],param_sig_i_t_dep[[4]],
                           param_sig_i_t_dep[[5]],param_sig_i_t_dep[[6]],
                           param_sig_i_t_dep[[7]],param_sig_i_t_dep[[8]])
    }
    
    if(time_dep_sig_param){
      rho0_t_dep = param_t_dep_func(list(1,16000),list(param_sig[1],0.5*param_sig[1]))
      sigma2_t_dep = param_t_dep_func(list(1,10000),list(param_sig[2],2*param_sig[2]))
      zeta_t_dep = param_t_dep_func(list(1,27000),list(param_sig[3],0.2))
      rho1_t_dep = param_t_dep_func(list(1,22000),list(param_sig[4],2*param_sig[4]))
      gamma_t_dep = param_t_dep_func(list(1,31000),list(param_sig[5],1.2))
      alpha_t_dep = param_t_dep_func(list(1,25000),list(param_sig[6],1.1))
      mu_x_t_dep = param_t_dep_func(list(1,13000),list(param_sig[7],.3*param_sig[7]))
      mu_y_t_dep = param_t_dep_func(list(1,19000),list(param_sig[8],param_sig[8]+0.2))
      
      param_sig_t_dep = param_sig_t_dep_func(rho0_t_dep,sigma2_t_dep,zeta_t_dep,
                                             rho1_t_dep,gamma_t_dep,alpha_t_dep,
                                             mu_x_t_dep,mu_y_t_dep)
    }
    
    ## Observation parameters ##
    mult_noise = FALSE
    time_dep_obs_noise = FALSE
    
    if(mult_noise){
      m = 3; 
      m_indices = list(coords_indices[1:3],coords_indices[4:6],coords_indices[7:10])
      tau2_vals = list(0.01,0.1,0.3)
    }
    
    if(!mult_noise){
      m = 1; m_indices = list(coords_indices[1:n_obs]); tau2_vals = list(0.01)
    }
    
    param_obs = param_obs_func(m,m_indices,tau2_vals)
    param_obs_t_dep =  param_t_dep_func(list(1),list(param_obs))
    
    if(time_dep_obs_noise){
      param_obs_t_dep =  param_t_dep_func(list(1,4000,7000,15000,19000),
                                          list(param_obs,
                                               param_obs_func(1,m_indices,list(0.05)),
                                               param_obs_func(1,m_indices,list(0.03)),
                                               param_obs_func(1,m_indices,list(0.1)),
                                               param_obs_func(1,m_indices,list(0.02))))
    }
    
    
    ## Observation bias ##
    mult_bias = FALSE
    if(mult_bias){
      m1 = 2
      m1_indices = list(coords_indices[1:5],coords_indices[6:10])
      bias_vals = list(0,2)
    }
    if(!mult_bias){
      m1 = 0; m1_indices = list(coords_indices[1:n_obs]); bias_vals = list(0)
    }
    param_bias = param_bias_func(m1,m1_indices,bias_vals)
    param_bias_t_dep = param_t_dep_func(list(1),list(param_bias))
    
    ## Parameter names ##
    param_sig_names = list(parse(text="rho[0]"),parse(text="sigma^2"),
                           parse(text="zeta"),parse(text="rho[1]"),
                           parse(text="gamma"),parse(text="alpha"),
                           parse(text="mu[x]"),parse(text="mu[y]"))
    param_obs_names = list(parse(text="tau^2"))
    param_bias_names = list(parse(text="beta"))
    
    ## Number of parameters ##
    n_param = length(param_sig)+param_obs$m+param_bias$m
  }
  
  ## Signal Noise Weights ##
  if(TRUE){
    
    ## Weighting Matrix (Physical Space) ##
    grid = expand.grid(x = 0:(n-1),y = 0:(n-1))/n
    centre_x = 6/12; centre_y = 7/12; sigma_x = 0.3; sigma_y = 0.3; rho=0
    weights = with(grid,sech2d(x = x,y = y, centre_x = centre_x,
                               centre_y = centre_y, sigma_x = sigma_x,
                               sigma_y=sigma_y, rho=rho))
    weights_mat = weights_mat_func(weights)
    
    ## Weighting Matrix (Spectral Space) ##
    spec_weights_mat = spec_weights_mat_func(n = n, K = K, 
                                             weights_mat = weights_mat)
    
    ## Alternative: No Weighting Matrix ##
    weights = NULL; spec_weights_mat = NULL
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
  
  ## Simulate Y_t ##
  if(TRUE){
    spde_sim = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                               param_bias_t_dep,seed = 2,n_obs = n_obs,
                               obs_loc = obs_loc, weights = weights,dt=dt)
    y_sim = spde_sim$y
  }
  
  ## Weighting matrix ##
  if(TRUE){
    W_coords = NULL; W = diag(1,n^2)
    coords2 = lapply(1:n^2,function(i) as.numeric(expand.grid(0:(n-1),0:(n-1))[i,]))
    Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs=n^2,coords = coords2)
    W_fourier = t(Phi)%*%W%*%Phi
  }
  
  ## Initial parameters ##
  if(TRUE){
    param0 = c(rho0 = 0.25, sigma2 = 0.5, zeta = 0.3,rho1 = 0.2, gamma = 1.5, 
               alpha = pi/3, mu_x = 0.1, mu_y = -0.15, tau2 = 0.1)
    param_sig0 = param0[1:8]
    #param_sig0[6:8] = param_sig[6:8]
    
    if(mult_noise){
      param_obs0 = param_obs_func(m,m_indices,list(param0[9],0.3,0.1))
    }
    if(!mult_noise){
      param_obs0 = param_obs_func(m,m_indices,list(param0[9]))
      #param_obs0 = param_obs
    }
    if(mult_bias){
      param_bias0 = param_bias_func(m1,m1_indices,list(0.9,1.1))
    }
    if(!mult_bias){
      param_bias0 = param_bias_func(m1,m1_indices,list(0))
    }
  }
  
  ## Log indices ##
  if(TRUE){
    log_indices = c(1,2,3,4,5,8+(1:param_obs0$m))
  }
  
  ## Parameter plot limits ##
  if(TRUE){
    
    plot_limits_auto = FALSE
    if(plot_limits_auto){
      plot_limits = cbind(0.5*sapply(1:8, function(i) min(abs(param_sig0[[i]]),abs(param_sig[[i]]))*abs(param_sig[[i]])/param_sig[[i]]),
                          1.5*sapply(1:8, function(i) max(abs(param_sig0[[i]]),abs(param_sig[[i]]))*abs(param_sig[[i]])/param_sig[[i]]))
      plot_limits = rbind(plot_limits,
                          cbind(0,1.5*max(unlist(param_obs$tau2_vals),unlist(param_obs0$tau2_vals))))    
      plot_limits = rbind(plot_limits,
                          cbind(0,1.5*max(unlist(param_bias$bias_vals),unlist(param_bias0$bias_vals)))) 
      plot_limits = t(apply(plot_limits,1,sort)) 
      plot_limits = lapply(1:dim(plot_limits)[1],function(i) plot_limits[i,])
    }
    
    if(!plot_limits_auto){
      plot_limits = rep(list(rep(0,2)),9)
      plot_limits[[1]] = c(0.2,0.6)
      plot_limits[[2]] = c(0.1,0.7)
      plot_limits[[3]] = c(0.1,0.6)
      plot_limits[[4]] = c(0.05,0.4)
      plot_limits[[5]] = c(0.9,2.3)
      plot_limits[[6]] = c(0.6,1.2)
      plot_limits[[7]] = c(0.05,0.35)
      plot_limits[[8]] = c(-0.35,0)
      plot_limits[[9]] = c(0,0.5)
      plot_limits[[10]] = c(0,2.5)
    }
    
    
  }
  
  ## Parameters to estimate ##
  if(TRUE){
    grad_params = 1:9
    
    #grad_params = 9:13
  }
  
  ## Hessian indices ##
  if(TRUE){
    hessian_indices = NULL
  }
  
  ## Number of iterations ##
  if(TRUE){
    n_iterations = 1
  }
  
  ## Step sizes ##
  if(TRUE){
    
    ## Initial step-size vector ##
    if(mult_noise && mult_bias){
      step_size_vec = c(0.002,0.005,0.005,0.001,0.015,0.005,0.001,0.001,0.0004,0.0006,0.0003,0.0004,0.0004)
    }
    if(!mult_noise || !mult_bias){
      step_size_vec = c(0.001,0.005,0.003,0.002,0.02,0.005,0.0003,0.0003,0.0005) # use for time_indep sim
      #step_size_vec = c(0.001,0.002,0.0015,0.0015,0.006,0.0015,0.0002,0.0002,0.0001) # use for time_dep sim
    
      step_size_vec[1] = 1/(dt^0.2)*step_size_vec[1]
      step_size_vec[2] = 1/(dt^0.1)*step_size_vec[2]
      step_size_vec[3] = 1/(dt^0.4)*step_size_vec[3]
      step_size_vec[4] = 1/(dt^0.1)*step_size_vec[4]
      step_size_vec[5] = 1/(dt^0.2)*step_size_vec[5]
      step_size_vec[6] = 1/(dt^0.2)*step_size_vec[6]
      step_size_vec[7:8] = 1/(dt^0.1)*step_size_vec[7:8]
    }
    
    ## Contant step size? ##
    constant_step_size = TRUE
    
    ## Step size matrix ##
    step_size = matrix(0,nrow=n_iterations*t,ncol=n_param)
    for (i in 1:dim(step_size)[1]){
      
      if(constant_step_size){
        step_size[i,] = step_size_vec
      }
      if(!constant_step_size){
        step_size[i,] = 10*step_size_vec*rep(i^(-0.5),n_param) 
      }
    }
    
    step_size[500:10000,9] = 0.1*step_size[500:10000,9]
  }
  
  ## Parameter estimation period ##
  if(TRUE){
    if(constant_step_size){
      if(mult_noise && mult_bias){
        param_est_starts = list(3000,2000,3000,3000,4000,2000,2000,2000,3000,3000,3000,3000,3000)
      }
      if(!mult_noise || !mult_bias){
        param_est_starts = list(8000,8000,8000,8000,8000,8000,8000,8000,8000,1000)
      }
    }
    if(!constant_step_size){
      if(mult_noise && mult_bias){
        param_est_starts = list(5000,5000,10000,5000,10000,4000,4000,4000,3000,3000,3000,3000,3000)
      }
      if(!mult_noise || !mult_bias){
        param_est_starts = list(10000,10000,20000,5000,15000,15000,10000,10000,1000)
      }
    }
    param_est_ends = rep(list(n_iterations*t),n_param)
    
    if(time_dep_sig_param && time_dep_obs_noise){
      param_est_starts = list(c(4000,19000),c(4000,12000),c(5000,31000),c(4000,24000),
                              c(5000,33000),c(5000,27000),c(4000,15000),c(4000,21000),
                              c(500,4500,7500,15500,19500))
      param_est_ends = list(c(16000,t*n_iterations),c(10000,t*n_iterations),
                            c(27000,t*n_iterations),c(22000,t*n_iterations),
                            c(31000,t*n_iterations),c(25000,t*n_iterations),
                            c(13000,t*n_iterations),c(19000,t*n_iterations),
                            c(4000,7000,15000,19000,t*n_iterations))
    }
  }
  
  ## Plot vertical (burn in) lines? ##
  if(TRUE){
    plot_vertical = FALSE
  }
  
  ## Plot Point Frequency ##
  if(TRUE){
    plot_point_freq = 20
  }
  
  ## Plot Frequency ##
  if(TRUE){
    plot_freq = 500
  }
  
  ## Filename ##
  if(TRUE){
    filename = paste("n",n,"_t",t,"_nobs",n_obs,sep="")
    plot_filename = paste(filename,".pdf",sep="")
  }
  
  ## Save ##
  if(TRUE){
    save = TRUE 
    save_plot = TRUE
  }
  
  ## Print Iteration Frequency ##
  if(TRUE){
    print_iter_freq = 100
  }
  
  ## Back Track Line Search? ##
  if(TRUE){
    back_track = FALSE
  }
  
  ## Back Track Line Search Tolerance ##
  if(TRUE){
    back_track_tol = 1e-3
  }
  
  ## Back Track Line Search for Each Param? ##
  if(TRUE){
    back_track_independent = FALSE
    back_track_param = 1:n_param
  }
  
  ## Back Track Line Search Time Independent? ##
  if(TRUE){
    back_track_time_independent = FALSE
  }
  
  ## Back Track Scale ##
  if(TRUE){
    back_track_scale = 0.5
  }
}

## RML ##
if(TRUE){
  if(K==n^2 && n_obs == n^2){
    RML_test = RML_spectral(w = y_sim, n = n, t = t, K = K, 
                            param_obs0 = param_obs0, param_sig0 = param_sig0, n_param = n_param,
                            grad_params = grad_params, step_size = step_size,
                            log_scale = FALSE, log_indices = log_indices,
                            hessian_indices = hessian_indices, 
                            n_iterations = n_iterations,
                            plot = TRUE, plot_limits = plot_limits, 
                            param_obs_true = param_obs_t_dep,
                            param_sig_true = param_sig_t_dep, 
                            param_obs_names = param_obs_names,
                            param_sig_names = param_sig_names,
                            param_est_starts = param_est_starts,
                            param_est_ends = param_est_ends,
                            plot_point_freq = plot_point_freq, plot_freq = plot_freq,
                            save_plot = save_plot, plot_filename = plot_filename, 
                            save = save, filename = filename, mean_param = TRUE,
                            print_iter_freq = print_iter_freq,
                            back_track = back_track, 
                            back_track_tol = back_track_tol,
                            back_track_independent = back_track_independent,
                            back_track_time_independent = back_track_time_independent,
                            back_track_scale = back_track_scale,
                            back_track_param = back_track_param,
                            plot_vertical = plot_vertical, 
                            W_fourier = W_fourier, dt = dt)
  }
  if(K<n^2 || n_obs < n^2){
    RML_test = RML(y = y_sim, n = n, t = t, K = K, 
                   param_obs0 = param_obs0, param_sig0 = param_sig0,
                   param_bias0 = param_bias0, n_param = n_param,
                   grad_params = grad_params, step_size = step_size,
                   log_scale = FALSE, log_indices = log_indices,
                   hessian_indices = hessian_indices, 
                   n_iterations = n_iterations,n_obs = n_obs, coords = coords,
                   spec_weights_mat = spec_weights_mat, plot = TRUE, 
                   plot_limits = plot_limits, 
                   param_obs_true = param_obs_t_dep,
                   param_sig_true = param_sig_t_dep, 
                   param_bias_true = param_bias_t_dep,
                   param_obs_names = param_obs_names,
                   param_sig_names = param_sig_names,
                   param_bias_names = param_bias_names, 
                   param_est_starts = param_est_starts,
                   param_est_ends = param_est_ends,
                   plot_point_freq = plot_point_freq, plot_freq = plot_freq,
                   save_plot = save_plot, plot_filename = plot_filename, 
                   save = save, filename = filename, mean_param = TRUE,
                   print_iter_freq = print_iter_freq,
                   back_track = back_track, back_track_tol = back_track_tol,
                   back_track_independent = back_track_independent,
                   back_track_time_independent = back_track_time_independent,
                   back_track_scale = back_track_scale,
                   back_track_param = back_track_param,
                   plot_vertical = plot_vertical,
                   W_fourier = W_fourier, dt = dt)
    
  }
}

#############

######################
#### BatchOptSens ####
######################

## Setup ##
if(TRUE){
  
  ## Prelims ##
  if(TRUE){
    n = 12; t = 15; dt = 1; K = 10; nu = 1
  }
  
  ## Initial Observation Locations ##
  if(TRUE){
    
    ## Coordinates ##
    obs0 = list(c(8,5),c(4,6),c(5,4),c(7,4))
    obs0 = list(c(12,12),c(1,1),c(1,11),c(11,1),c(11,11))
    
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
           xlim = c(0,n),ylim=c(0,n),xlab=expression(x),ylab=expression(y),
           main="Sensor Locations")
      abline(v=0:(n),h=(0:n),col="gray92",lty=2)
    }
  }
  
  ## Parameters ##
  if(TRUE){
    
    ## Signal parameters ##
    param_sig = c(rho0 = 0.5, sigma2 = 0.2, zeta = 0.5,rho1 = 0.1, gamma = 2, 
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
    
    ## Parameter names ##
    param_sig_names = list(parse(text="rho[0]"),parse(text="sigma^2"),
                           parse(text="zeta"),parse(text="rho[1]"),
                           parse(text="gamma"),parse(text="alpha"),
                           parse(text="mu[x]"),parse(text="mu[y]"))
    param_obs_names = list(parse(text="tau^2"))
    param_bias_names = list(parse(text="beta"))
    
    ## Number of parameters ##
    n_param = length(param_sig)+param_obs$m+param_bias$m
  }
  
  ## Signal Noise Weights ##
  if(TRUE){
    grid = expand.grid(x = 0:(n-1),y = 0:(n-1))/n
    centre_x = 4/12; centre_y = 5/12; sigma_x = 0.2; sigma_y = 0.2; rho = 0 
    weights = with(grid,sech2d(x = x,y = y, centre_x = centre_x,
                                    centre_y = centre_y, sigma_x = sigma_x,
                                    sigma_y=sigma_y, rho = rho))
    weights_mat = weights_mat_func(weights)
    spec_weights_mat = spec_weights_mat_func(n = n, K = K, 
                                             weights_mat = weights_mat)
    weights = NULL; spec_weights_mat = NULL
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
  
  ## Simulate Y_t ##
  if(TRUE){
    spde_sim = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                               param_bias_t_dep,seed = 1,n_obs = n_obs,
                               obs_loc = obs_loc, weights = weights, dt = dt)
    y_sim = spde_sim$y
  }
  
  ## Weighting matrix ##
  if(TRUE){
    
    ## Weighting coordinates ##
    W_coords = list(c(10,11),c(7,10),c(6,8),c(9,6))
    
    ## Weighting indices ##
    W_indices = lapply(W_coords, function(A) A[2]*n+A[1]+1)
    
    ## Weighting matrix (physical space) ##
    W = matrix(0,nrow=n^2,ncol=n^2)
    diag(W)[unlist(W_indices)] = rep(10,length(W_coords))
    
    ## Alternative: no weighting ##
    W_coords = NULL; diag(W) = rep(1,n^2); 
    
    ## Weighting matrix (spectral space) ##
    coords = lapply(1:n^2,function(i) as.numeric(expand.grid(0:(n-1),0:(n-1))[i,]))
    Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs=n^2,coords = coords)
    W_fourier = t(Phi)%*%W%*%Phi
  }
  
  ## Number of iterations ##
  if(TRUE){
    n_iterations = 100
  }
  
  ## Step sizes ##
  if(TRUE){
    
    ## Initial step size vector ##
    step_size_vec = c(0.05,0.05)
    
    ## Step size matrix ##
    step_size_mat = matrix(0,nrow=n_iterations,ncol=2)
    for (i in 1:n_iterations){
      step_size_mat[i,]=step_size_vec
    }
    
    ## Step size matrices for all sensors ##
    step_size = rep(list(step_size_mat),n_obs)
  }
  
  ## Observations to be optimised ##
  if(TRUE){
    grad_obs = 1
  }
  
  ## Plot Point Frequency ##
  if(TRUE){
    plot_point_freq = 2
  }
  
  ## Plot Frequency ##
  if(TRUE){
    plot_freq = 50
  }
  
  ## Filenames ##
  if(TRUE){
    filename_plot = "offline_opt_sens.pdf"
    filename_plot2d = "offline_opt_sens2d.pdf"
  }
  
  ## Save plots? ##
  if(TRUE){
    save_plot = FALSE
    save_plot2d = FALSE
  }
  
  ## Plot grid? ##
  if(TRUE){
    grid = TRUE
  }
  
  ## Print Iteration Frequency ##
  if(TRUE){
    print_iter_freq = 10
  }
  
  ## Back Track Line Search? ##
  if(TRUE){
    back_track = FALSE
  }
  
  ## Back Track Line Search Tolerance ##
  if(TRUE){
    back_track_tol = 1e-4
  }
  
  ## Back Track Line Search for Each Param? ##
  if(TRUE){
    back_track_independent = TRUE
    back_track_obs = 1:n_obs
  }
  
  ## Back Track Line Search Time Independent? ##
  if(TRUE){
    back_track_time_independent = TRUE
  }
  
  ## Back Track Scale ##
  if(TRUE){
    back_track_scale = 0.8
  }
  
  ## Highlight initial point ##
  if(TRUE){
    plot_initial = TRUE
  }
  
  ## Add legend? ##
  if(TRUE){
    leg = TRUE
  }
  
  
}

## BOSP ##
if(TRUE){
  BatchOptSens_test = BatchOptSens(y = y_sim,n_obs = n_obs,t = t,K = K,
                                   param_obs = param_obs, 
                                   param_sig = param_sig,
                                   param_bias = param_bias,
                                   obs0 = obs0,
                                   n_iterations = n_iterations,
                                   step_size = step_size,
                                   W_fourier = W_fourier,
                                   spec_weights_mat=spec_weights_mat,
                                   plot=TRUE,W_coords=W_coords,
                                   grad_obs=grad_obs, plot2d=TRUE,
                                   plot_point_freq = plot_point_freq,
                                   plot_freq = plot_freq,save_plot=save_plot,
                                   filename_plot=filename_plot,
                                   save_plot2d=save_plot2d,
                                   filename_plot2d=filename_plot2d,
                                   grid = grid,
                                   print_iter_freq = print_iter_freq,
                                   back_track = back_track, 
                                   back_track_tol = back_track_tol,
                                   back_track_independent = back_track_independent,
                                   back_track_time_independent = back_track_time_independent,
                                   back_track_scale = back_track_scale,
                                   back_track_obs = back_track_obs,
                                   plot_initial = plot_initial,
                                   leg = leg, dt = dt)
}

######################

#######################
#### RecursOptSens ####
#######################

## Setup ##
if(TRUE){
  
  ## Prelims ##
  if(TRUE){
    n = 12; t = 2000; dt = 1; K = 20; nu = 1
  }
  
  ## Initial Observation Locations ##
  if(TRUE){
    
    ## Coordinates (weighted coordinates simulation) ##
    obs0 = list(c(10.1,7.8),c(4.1,6.01),c(5.2,3.75),c(7.2,4.02),c(3.2,3.1),
                c(6.1,2.1),c(1.01,2.8),c(3,1))
    
    ## Coordinates (weighted signal noise simulation) ##
    #obs0 = list(c(1.7,2.2),c(1.5,9.6),c(11.1,3.7),c(10.1,9.8))
    #obs0 = list(c(3.595071,1.032352),
    #            c(18.60910,11.09808),
    #            c(13.557721,1.618061),
    #            c(8.638112,11.551511))
    
    #obs0 = list(c(5,5),c(15,5),c(5,15),c(17,15),c(10,10))
    #obs0 = list(c(5,5),c(15,5),c(5,15),c(15,15),c(10,20))
    
    #obs0 = list(c(5,5),c(5.5,5.5))#,c(2,5),c(6.5,6.5),c(7,7))
    
    ## Indices ##
    obs0_indices = sapply(obs0,function(x) x[1]+n*x[2]+1)
    
    ## Number of observations ##
    n_obs = length(obs0)
    
    ## Observation location at all times ##
    obs_loc = lapply(1:n_obs,
                     function(j) matrix(obs0[[j]],nrow=t,ncol=2,byrow=T))
    
    ## Plot initial observation locations ##
    plot_obs = FALSE
    if(plot_obs){
      par(mfrow=c(1,1))
      plot(sapply(obs0,function(x) x)[1,],sapply(obs0,function(x) x)[2,],
           xlim = c(0,n),ylim=c(0,n),xlab=expression(x),ylab=expression(y),
           main="Sensor Locations")
      abline(v=0:(n),h=(0:n),col="gray92",lty=2)
    }
    
    ## 
  }
  
  ## Parameters ##
  if(TRUE){
    
    ## Signal parameters ##
    param_sig = c(rho0 = 1.00, sigma2 = 0.2, zeta = 0.5, rho1 = 0.1, gamma = 2, 
                  alpha = pi/2, mu_x = 0.3, mu_y = -0.3)
    
    param_sig_i_t_dep = lapply(1:length(param_sig), 
                               function(i) param_t_dep_func(1,param_sig[i]))
    
    param_sig_t_dep = 
      param_sig_t_dep_func(param_sig_i_t_dep[[1]],param_sig_i_t_dep[[2]],
                           param_sig_i_t_dep[[3]],param_sig_i_t_dep[[4]],
                           param_sig_i_t_dep[[5]],param_sig_i_t_dep[[6]],
                           param_sig_i_t_dep[[7]],param_sig_i_t_dep[[8]])
    
    ## Observation parameters ##
    m = 1; m_indices = list(obs0_indices); tau2_vals = list(0.1)
    param_obs = param_obs_func(m,m_indices,tau2_vals)
    param_obs_t_dep =  param_t_dep_func(list(1),list(param_obs))
    
    ## Observation bias ##
    m1 = 0; m1_indices = list(obs0_indices); bias_vals = list(0)
    param_bias = param_bias_func(m1,m1_indices,bias_vals)
    param_bias_t_dep = param_t_dep_func(list(1),list(param_bias))
    
    ## Parameter names ##
    param_sig_names = list(parse(text="rho[0]"),parse(text="sigma^2"),
                           parse(text="zeta"),parse(text="rho[1]"),
                           parse(text="gamma"),parse(text="alpha"),
                           parse(text="mu[x]"),parse(text="mu[y]"))
    param_obs_names = list(parse(text="tau^2"))
    param_bias_names = list(parse(text="beta"))
    
    ## Number of parameters ##
    n_param = length(param_sig)+param_obs$m+param_bias$m
  }
  
  ## Signal Noise Weights ##
  if(TRUE){
    
    ## Weighting Matrix (Physical Space) ##
    grid = expand.grid(x = 0:(n-1),y = 0:(n-1))/n
    centre_x = 6/12; centre_y = 6/12; sigma_x = 0.2; sigma_y = 0.2; rho = 0
    weights = with(grid,sech2d(x = x,y = y, centre_x = centre_x,
                               centre_y = centre_y, sigma_x = sigma_x,
                               sigma_y=sigma_y, rho=rho))
    grid$weights = weights
    weights_mat = weights_mat_func(weights)
    
    ## Weighting Matrix (Spectral Space) ##
    spec_weights_mat = spec_weights_mat_func(n = n, K = K, 
                                             weights_mat = weights_mat)
    
    
    ## Plot ##
    if(FALSE){
      pdf("signal_weights.pdf",width=6,height=6)
      wireframe(weights~x*y, data = grid, cuts = 100, main="", 
                col.regions = heat.colors(1000)[1:length(heat.colors(1000))],
                drape = TRUE,zlab="z",xlab="x",ylab="y", light.source = c(10,10,10), 
                aspect = c(1,1), screen = list(z = 20, x = -60), col = F, 
                zoom = 1, scales = list(arrows = F, tck = c(0.8, 0.6, 0.4), 
                                        distance =c(.7, .7, .5)))
      dev.off()
    }
    
    ## Alternative: No Weighting Matrix ##
    weights = NULL; weights_mat = NULL; spec_weights_mat = NULL
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
  
  ## Simulate Y_t ##
  if(TRUE){
    spde_sim = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                               param_bias_t_dep,seed = 1,n_obs = n_obs,
                               obs_loc = obs_loc, weights = weights,dt=dt)
    y_sim = spde_sim$y
  }
  
  ## Weighting matrix ##
  if(TRUE){
    
    ## Weighting coordinates ##
    W_coords = list(c(9,6),c(7,10),c(6,8),c(0,7),c(4,4),c(10,11),c(1,1),c(3,10))
    
    ## Weighting indices ##
    W_indices = lapply(W_coords, function(A) A[2]*n+A[1]+1)
    
    ## Weighting matrix (physical space) ##
    W = matrix(0,nrow=n^2,ncol=n^2)
    diag(W)[unlist(W_indices)] = rep(100,length(W_coords))
    
    ## Alternative: no weighting ##
    if(FALSE){
      W_coords = NULL; diag(W) = rep(10,n^2)
    }
    
    ## Weighting matrix (spectral space) ##
    coords = lapply(1:n^2,function(i) as.numeric(expand.grid(0:(n-1),0:(n-1))[i,]))
    Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs=n^2,coords = coords)
    W_fourier = t(Phi)%*%W%*%Phi
  }
  
  ## Number of iterations ##
  if(TRUE){
    n_iterations = 1
  }
  
  ## Step sizes ##
  if(TRUE){
    
    ## Initial step size vector ##
    #step_size_vec = 0.1*c(1,1)
    step_size_vec = c(10,10)
    
    ## Contant step size? ##
    constant_step_size = TRUE
    
    ## Step size matrix ##
    step_size_mat = matrix(0,nrow=n_iterations*t,ncol=2)
    for (i in 1:(dim(step_size_mat)[1])){
      
      if(constant_step_size){
        step_size_mat[i,] = step_size_vec
      }
      if(!constant_step_size){
        step_size_mat[i,] = step_size_vec*rep(i^(-0.5),2) 
      }
      
    }
    
    ## Step size matrices for all sensors ##
    step_size = rep(list(step_size_mat),n_obs)
  }
  
  ## Plots? ##
  if(TRUE){
    plot = TRUE
    plot2d = FALSE
  }
  ## Observations to be optimised ##
  if(TRUE){
    grad_obs = 1:n_obs
  }
  
  ## Plot Point Frequency ##
  if(TRUE){
    plot_point_freq = 5
  }
  
  ## Plot Frequency ##
  if(TRUE){
    plot_freq = 500
  }
  
  ## Filenames ##
  if(TRUE){
    filename_plot = "ROSP_plot.pdf"
    filename_plot2d = "ROSP_plot2d.pdf"
  }
  
  
  ## Save plots? ##
  if(TRUE){
    save_plot = TRUE
    save_plot2d = TRUE
  }
  
  ## Plot grid? ##
  if(TRUE){
    grid = TRUE
  }
  
  ## Print Iteration Frequency ##
  if(TRUE){
    print_iter_freq = 100
  }
  
  
  ## Back Track Line Search? ##
  if(TRUE){
    back_track = FALSE
  }
  
  ## Back Track Line Search Tolerance ##
  if(TRUE){
    back_track_tol = 1e-7
  }
  
  ## Back Track Line Search for Each Param? ##
  if(TRUE){
    back_track_independent = TRUE
    back_track_obs = 1:n_obs
  }
  
  ## Back Track Line Search Time Independent? ##
  if(TRUE){
    back_track_time_independent = TRUE
  }
  
  ## Back Track Scale ##
  if(TRUE){
    back_track_scale = 0.5
  }
  
  ## Highlight initial point ##
  if(TRUE){
    plot_initial = TRUE
  }
  
  ## Add legend? ##
  if(TRUE){
    leg = TRUE
  }
  
}

## ROSP ##
if(TRUE){
  RecursOptSens_test = RecursOptSens(y_sim,n_obs,t,K,param_obs,param_sig,
                                    param_bias,obs0,n_iterations,step_size,
                                    W_fourier = W_fourier,
                                    spec_weights_mat=spec_weights_mat,
                                    plot=plot,W_coords=W_coords,
                                    grad_obs=grad_obs, plot2d=plot2d,
                                    plot_point_freq = plot_point_freq,
                                    plot_freq = plot_freq,save_plot=save_plot,
                                    filename_plot=filename_plot,
                                    save_plot2d=save_plot2d,
                                    filename_plot2d=filename_plot2d,
                                    grid = grid,
                                    print_iter_freq = print_iter_freq,
                                    back_track = back_track, 
                                    back_track_tol = back_track_tol,
                                    back_track_independent = back_track_independent,
                                    back_track_time_independent = back_track_time_independent,
                                    back_track_scale = back_track_scale,
                                    back_track_obs = back_track_obs,
                                    plot_initial = plot_initial,
                                    leg = leg,
                                    spde_sim=spde_sim, dt = dt)
}

#######################

##################
#### RML_ROSP ####
##################

## Setup ##
if(TRUE){
  
  ## Prelims ##
  if(TRUE){
    n = 12; t = 5000; dt = 1; K = 10; nu = 1
  }
  
  ## RML & OSP Times ##
  if(TRUE){
    t_RML = 8; t_ROSP = 2
  }
  
  ## Initial Observation Locations ##
  if(TRUE){
    
    ## Coordinates ##
    obs0 = list(c(10.1,7.8),c(4.1,6.01),c(5.2,3.75),c(7.2,4.02),c(3.2,3.1),
                c(6.1,2.1),c(1.01,2.8),c(3,1))
    
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
    param_sig = c(rho0 = 0.5, sigma2 = 0.2, zeta = 0.5,rho1 = 0.1, gamma = 2, 
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
    param0 = c(rho0 = 0.25, sigma2 = 0.8, zeta = 0.1,rho1 = 0.2, gamma = 1.2, 
               alpha = pi/2, mu_x = 0.1, mu_y = -0.15, tau2 = 0.1)
    param_sig0 = param0[1:8]
    param_obs0 = param_obs_func(m,m_indices,list(param0[9]))
    param_bias0 = param_bias_func(m1,m1_indices,list(0))
  }
  
  ## Signal Noise Weights ##
  if(TRUE){
    
    ## Weighting Matrix (Physical Space) ##
    grid = expand.grid(x = 0:(n-1),y = 0:(n-1))/n
    centre_x = 6/12; centre_y = 7/12; sigma_x = 0.3; sigma_y = 0.3; rho = 0
    weights = with(grid,sech2d(x = x,y = y, centre_x = centre_x,
                               centre_y = centre_y, sigma_x = sigma_x,
                               sigma_y=sigma_y, rho = rho))
    weights_mat = weights_mat_func(weights)
    
    ## Weighting Matrix (Spectral Space) ##
    spec_weights_mat = spec_weights_mat_func(n = n, K = K, 
                                             weights_mat = weights_mat)
    
    ## Alternative: No Weighting Matrix ##
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
  if(TRUE){
    spde_sim = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                               param_bias_t_dep,seed = 1,n_obs = n_obs,
                               obs_loc = obs_loc, alpha_sim = TRUE,
                               x_sim = TRUE, obs_noise_sim = TRUE,
                               obs_bias_sim = TRUE, y_sim = FALSE, dt = dt)
  }
  
  ## Log scale & log indices ##
  if(TRUE){
    log_scale = FALSE
    log_indices = c(1,2,3,4,5,8+(1:param_obs0$m))
  }
  
  ## Parameters to estimate ##
  if(TRUE){
    grad_params = 1:n_param
  }
  
  ## Observations to be optimised ##
  if(TRUE){
    grad_obs = 1:n_obs
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
    
    if(dt==1){
      ## Initial step-size vector ##
      step_size_vec = c(0.005,0.001,0.01,0.005,0.015,0.02,0.1,0.1,0.00005)#*20
      step_size_vec = c(0.002,0.002,0.01,0.001,0.02,0.01,0.001,0.001,0.0001)
      
      ## Step size matrix ##
      step_size_RML = matrix(0,nrow=n_iterations*t,ncol=n_param)
      for (i in 1:dim(step_size_RML)[1]){
        step_size_RML[i,] = step_size_vec
        step_size_RML[i,] = 20*step_size_vec*i^(-0.51)
        step_size_RML[i,4] = 10*step_size_vec[4]*i^(-0.51)
        step_size_RML[i,5] = 1*step_size_vec[5]*i^(-0.51)
        step_size_RML[i,6] = 100*step_size_vec[6]*i^(-0.51)
        step_size_RML[i,7:8] = 100*step_size_vec[7:8]*i^(-0.51)
      }
    
      for (i in 1000:dim(step_size_RML)[1]){
        step_size_RML[i,4] = 40*step_size_vec[4]*(i)^(-0.51)
        step_size_RML[i,5] = 40*step_size_vec[5]*(i)^(-0.51)
        step_size_RML[i,7:8] = 50*step_size_vec[7:8]*(i)^(-0.51)
      }
    }
    
    if(dt==0.1){
      ## Initial step-size vector ##
      step_size_vec = c(0.005,0.001,0.01,0.005,0.015,0.02,0.1,0.1,0.00005)#*20
      step_size_vec = c(0.002,0.002,0.01,0.001,0.02,0.01,0.001,0.001,0.0005)
      
      ## Step size matrix ##
      step_size_RML = matrix(0,nrow=n_iterations*t,ncol=n_param)
      for (i in 1:dim(step_size_RML)[1]){
        step_size_RML[i,] = step_size_vec
        step_size_RML[i,] = 20*step_size_vec*i^(-0.51)
        step_size_RML[i,3] = 80*step_size_vec[3]*i^(-0.51)
        step_size_RML[i,4] = 30*step_size_vec[4]*i^(-0.51)
        step_size_RML[i,5] = 40*step_size_vec[5]*i^(-0.51)
        step_size_RML[i,6] = 100*step_size_vec[6]*i^(-0.51)
        step_size_RML[i,7:8] = 100*step_size_vec[7:8]*i^(-0.51)
      }
      
      for (i in 2000:dim(step_size_RML)[1]){
        step_size_RML[i,4] = 40*step_size_vec[4]*(i)^(-0.51)
        step_size_RML[i,5] = 40*step_size_vec[5]*(i)^(-0.51)
        step_size_RML[i,7:8] = 50*step_size_vec[7:8]*(i)^(-0.51)
      }
    }
    
  }
  
  ## ROSP Step sizes ##
  if(TRUE){
    
    ## Initial step size vector ##
    step_size_vec = c(10,10)
    step_size_vec = c(1,1)
    
    ## Step size matrix ##
    step_size_mat = matrix(0,nrow=n_iterations*t,ncol=2)
    for (i in 1:(n_iterations*t)){
      step_size_mat[i,]=step_size_vec*i^(-0.505)
    }
    
    ## Step size matrices for all sensors ##
    step_size_ROSP = rep(list(step_size_mat),n_obs)
  }
  
  ## Weighting matrix ##
  if(TRUE){
    
    ## Weighting coordinates ##
    if(dt==1){
      W_coords = list(c(10,11),c(3,10),c(6,8),c(9,6),c(4,4),c(7,10),c(1,1),c(1,7))
    }
    
    if(dt==0.1){
      W_coords = list(c(10,11), c(6,8), c(4,4), c(1,1), c(3,10), c(9,6),
                      c(1,7),c(7,10))
    }
    
    ## Weighting coordinates at all times ##
    target_loc = lapply(1:n_obs,
                        function(j) matrix(W_coords[[j]],nrow=t,ncol=2,byrow=T))
    
    ## Weighting indices ##
    W_indices = lapply(W_coords, function(A) A[2]*n+A[1]+1)
    
    ## Weighting matrix (physical space) ##
    W = matrix(0,nrow=n^2,ncol=n^2)
    diag(W)[unlist(W_indices)] = rep(10,length(W_coords))
    
    ## Alternative: no weighting matrix ##
    if(FALSE){
      W = diag(1,n^2)
      W_coords = NULL
    }
    
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
    RML_plot_freq = 500
  }
  
  ## Save RML plot? ##
  if(TRUE){
    save_RML_plot = TRUE
  }
  
  ## RML Plot Filename ##
  if(TRUE){
    RML_plot_filename = "RML_plot.pdf"
  }
  
  ## ROSP Plots? ##
  if(TRUE){
    ROSP_plot = TRUE
    ROSP_plot2d = TRUE
  }
  
  ## ROSP Plot Point Frequency ##
  if(TRUE){
    ROSP_plot_point_freq = 20
  }
  
  ## ROSP Plot Frequency ##
  if(TRUE){
    ROSP_plot_freq = 2500
  }
  
  ## Save ROSP Plots? ##
  if(TRUE){
    save_ROSP_plot = TRUE
    save_ROSP_plot2d = TRUE
  }
  
  ## ROSP Plot Filenames ##
  if(TRUE){
    ROSP_plot_filename = "ROSP_plot.pdf"
    ROSP_plot2d_filename = "ROSP_plot2d.pdf"
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
  RML_ROSP_test = RML_ROSP(spde_sim = spde_sim, n = n, n_obs = n_obs, t = t, 
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
                               ROSP_leg = ROSP_leg,
                               y_sim_t_dep = TRUE, dt = dt)
}

##################

#######################
#### RML_ROSP (v2) ####
#######################

## Setup ##
if(TRUE){
  
  ## Prelims ##
  if(TRUE){
    n = 12; t = 5000; dt = 1; K = 10; nu = 1
  }
  
  ## RML & OSP Times ##
  if(TRUE){
    t_RML = 4; t_ROSP = 1
  }
  
  ## Initial Observation Locations ##
  if(TRUE){
    
    ## Coordinates ##
    obs0 = list(c(10.1,7.8),c(4.1,6.01),c(5.2,3.75),c(7.2,4.02),c(3.2,3.1),
                c(6.1,2.1),c(1.01,2.8),c(3,1))
    
    ## Indices ##
    obs0_indices = sapply(obs0,function(x) x[1]+n*x[2]+1)
    
    ## Number of observations ##
    n_obs = length(obs0)
    
    ## Observation location at all times ##
    obs_loc = lapply(1:n_obs,
                     function(j) matrix(obs0[[j]],nrow=t,ncol=2,byrow=T))
    
    ## Plot initial observation locations ##
    plot_obs = FALSE
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
    param_sig = c(rho0 = 0.5, sigma2 = 0.2, zeta = 0.5,rho1 = 0.1, gamma = 2, 
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
    param0 = c(rho0 = 0.25, sigma2 = 0.8, zeta = 0.1,rho1 = 0.2, gamma = 1.2, 
               alpha = pi/2, mu_x = 0.1, mu_y = -0.15, tau2 = 0.1)
    param_sig0 = param0[1:8]
    param_obs0 = param_obs_func(m,m_indices,list(param0[9]))
    param_bias0 = param_bias_func(m1,m1_indices,list(0))
  }
  
  ## Signal Noise Weights ##
  if(TRUE){
    
    ## Weighting Matrix (Physical Space) ##
    grid = expand.grid(x = 0:(n-1),y = 0:(n-1))/n
    centre_x = 6/12; centre_y = 7/12; sigma_x = 0.3; sigma_y = 0.3; rho = 0
    weights = with(grid,sech2d(x = x,y = y, centre_x = centre_x,
                               centre_y = centre_y, sigma_x = sigma_x,
                               sigma_y=sigma_y, rho = rho))
    weights_mat = weights_mat_func(weights)
    
    ## Weighting Matrix (Spectral Space) ##
    spec_weights_mat = spec_weights_mat_func(n = n, K = K, 
                                             weights_mat = weights_mat)
    
    ## Alternative: No Weighting Matrix ##
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
  if(TRUE){
    spde_sim = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                               param_bias_t_dep,seed = 1,n_obs = n_obs,
                               obs_loc = obs_loc, alpha_sim = TRUE,
                               x_sim = TRUE, obs_noise_sim = TRUE,
                               obs_bias_sim = TRUE, y_sim = FALSE, dt = dt)
  }
  
  ## Log scale & log indices ##
  if(TRUE){
    log_scale = FALSE
    log_indices = c(1,2,3,4,5,8+(1:param_obs0$m))
  }
  
  ## Parameters to estimate ##
  if(TRUE){
    grad_params = 1:n_param
  }
  
  ## Observations to be optimised ##
  if(TRUE){
    grad_obs = 1:n_obs
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
    step_size_vec = c(0.005,0.001,0.01,0.005,0.015,0.02,0.1,0.1,0.00005)#*20
    step_size_vec = c(0.002,0.002,0.01,0.001,0.02,0.01,0.001,0.001,0.0001)
    
    ## Step size matrix ##
    step_size_RML = matrix(0,nrow=n_iterations*t,ncol=n_param)
    for (i in 1:dim(step_size_RML)[1]){
      step_size_RML[i,] = step_size_vec
      step_size_RML[i,] = 20*step_size_vec*i^(-0.51)
      step_size_RML[i,4] = 10*step_size_vec[4]*i^(-0.51)
      step_size_RML[i,5] = 1*step_size_vec[5]*i^(-0.51)
      step_size_RML[i,6] = 100*step_size_vec[6]*i^(-0.51)
      step_size_RML[i,7:8] = 100*step_size_vec[7:8]*i^(-0.51)
    }
    
    for (i in 2000:dim(step_size_RML)[1]){
      step_size_RML[i,4] = 40*step_size_vec[4]*(i)^(-0.51)
      step_size_RML[i,5] = 40*step_size_vec[5]*(i)^(-0.51)
      step_size_RML[i,7:8] = 50*step_size_vec[7:8]*(i)^(-0.51)
    }
    
  }
  
  ## ROSP Step sizes ##
  if(TRUE){
    
    ## Initial step size vector ##
    step_size_vec = c(10,10)
    step_size_vec = c(1,1)
    
    ## Step size matrix ##
    step_size_mat = matrix(0,nrow=n_iterations*t,ncol=2)
    for (i in 1:(n_iterations*t)){
      step_size_mat[i,]=step_size_vec*i^(-0.505)
    }
    
    ## Step size matrices for all sensors ##
    step_size_ROSP = rep(list(step_size_mat),n_obs)
  }
  
  ## Weighting matrix ##
  if(TRUE){
    
    ## Weighting coordinates ##
    W_coords = list(c(10,11),c(4,4),c(9,6),c(6,8),c(1,7),c(7,10),
                    c(1,1),c(3,10))
    
    ## Weighting coordinates at all times ##
    target_loc = lapply(1:n_obs,
                     function(j) matrix(W_coords[[j]],nrow=t,ncol=2,byrow=T))
    
    ## Weighting indices ##
    W_indices = lapply(W_coords, function(A) A[2]*n+A[1]+1)
    
    ## Weighting matrix (physical space) ##
    W = matrix(0,nrow=n^2,ncol=n^2)
    diag(W)[unlist(W_indices)] = rep(10,length(W_coords))
    
    ## Alternative: no weighting matrix ##
    if(FALSE){
      W = diag(1,n^2)
      W_coords = NULL
    }
    
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
    RML_plot_filename = "RML_plot.pdf"
  }
  
  ## ROSP Plots? ##
  if(TRUE){
    ROSP_plot = TRUE
    ROSP_plot2d = FALSE
  }

  ## ROSP Plot Point Frequency ##
  if(TRUE){
    ROSP_plot_point_freq = 20
  }
  
  ## ROSP Plot Frequency ##
  if(TRUE){
    ROSP_plot_freq = 1000
  }
  
  ## Save ROSP Plots? ##
  if(TRUE){
    save_ROSP_plot = TRUE
    save_ROSP_plot2d = TRUE
  }
  
  ## ROSP Plot Filenames ##
  if(TRUE){
    ROSP_plot_filename = "ROSP_plot.pdf"
    ROSP_plot2d_filename = "ROSP_plot2d.pdf"
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

## RML_ROSP (v2) ##
if(TRUE){
  
  ## Additional setup for v2 ##
  if(TRUE){
  RML_print_iter_freq = NULL; ROSP_print_iter_freq = NULL
  RML_ROSP_print_iter_freq = 100
  y_sim_t_dep = TRUE
  }
  
  RML_ROSP_test2 = RML_ROSP_v2(spde_sim = spde_sim, n = n, n_obs = n_obs, t = t, 
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
                           ROSP_leg = ROSP_leg,
                           RML_ROSP_print_iter_freq = RML_ROSP_print_iter_freq,
                           y_sim_t_dep = y_sim_t_dep, dt = dt)
}

#######################
