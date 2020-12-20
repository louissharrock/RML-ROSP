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

##################
#### RML_ROSP ####
##################

## Initial Locations ##
n = 12
min_x = 0; space_x = 0.05; max_x = 1 - space_x
grid_x = seq(min_x,max_x,space_x); grid_x = n*grid_x
grid_y = grid_x
len_x = length(grid_x)
len_y = length(grid_y)
init_loc = lapply(1:(len_x*len_y),function(i) as.numeric(expand.grid(grid_x,grid_y)[i,]))
fig_wd = "/Users/ls616/Desktop/Plot5"
t = 2500
all_param_est_mux = matrix(0,ncol=len_x*len_y,nrow=t+1)
all_param_est_muy = matrix(0,ncol=len_x*len_y,nrow=t+1)
all_param_est_tau2 = matrix(0,ncol=len_x*len_y,nrow=t+1)
all_obs_x = matrix(0,ncol=len_x*len_y,nrow=t+1)
all_obs_y = matrix(0,ncol=len_x*len_y,nrow=t+1)

## Run RML_ROSP ##
for (j in 1:length(init_loc)){
  ## Setup ##
  if(TRUE){
    
    ## Prelims ##
    if(TRUE){
      n = n; t = t; dt = 1; K = 13; nu = 1
    }
    
    ## RML & OSP Times ##
    if(TRUE){
      t_RML = 2; t_ROSP = 8
    }
    
    ## Initial Observation Locations ##
    if(TRUE){
      
      ## Coordinates ##
      #obs0 = list(c(10.1,7.8),c(4.1,6.01),c(5.2,3.75),c(7.2,4.02),c(3.2,3.1),
      #            c(6.1,2.1),c(1.01,2.8),c(3,1))
      
      ## Coordinates (weighted signal noise simulation) ##
      obs0 = list(c(1,1),c(1,11),c(11,1),c(11,11),#c(2,3),c(2.5,11.5),c(11,2.5),c(9,9),
                  #c(1,1),
                  c(1,3),c(1,5),c(1,7),c(1,9),#c(1,11),#c(1,2),c(1,4),c(1,6),c(1,8),c(1,10),
                  c(3,1),c(5,1),c(7,1),c(9,1),#c(11,1),#,c(2,1),#c(4,1),c(6,1),c(8,1),c(10,1),
                  c(11,1),c(11,3),c(11,5),c(11,7),c(11,9),c(11,11),#c(11,2),c(11,4),c(11,6),c(11,8),c(11,10)
                  c(3,11),c(5,11),c(7,11),c(9,11))#,c(11,11))#,c(2,11),#c(4,11),c(6,11),c(8,11),c(10,11),
      
      obs0 = list(init_loc[[j]],#,c(6,8),
                  c(2,2),c(2,6),c(2,10),
                  c(6,2),c(6,6),c(6,10),
                  c(10,2),c(10,6),c(10,10))
                  #c(1.5,1.5),c(4.5,1.5),c(7.5,1.5),c(10.5,1.5),
                  #c(1.5,4.5),c(4.5,4.5),c(7.5,4.5),c(10.5,4.5),
                  #c(1.5,7.5),c(4.5,7.5),c(7.5,7.5),c(10.5,7.5),
                  #c(1.5,10.5),c(4.5,10.5),c(7.5,10.5),c(10.5,10.5))
      
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
        abline(v=0:(n-1),h=(0:n-1),col="gray92",lty=2)
      }
      
      ## 
    }
    
    ## True Parameters ##
    if(TRUE){
      
      ## Signal parameters ##
      param_sig = c(rho0 = 0.5, sigma2 = 0.2, zeta = 0.5,rho1 = 0.1, gamma = 2, 
                    alpha = pi/4, mu_x = 0.1, mu_y = -0.1)
      
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
      param0 = c(rho0 = 0.5, sigma2 = 0.2, zeta = 0.5,rho1 = 0.1, gamma = 2, 
                 alpha = pi/4, mu_x = 0.39, mu_y = -0.41, tau2 = 0.5)
      param_sig0 = param0[1:8]
      param_obs0 = param_obs_func(m,m_indices,list(param0[9]))
      param_bias0 = param_bias_func(m1,m1_indices,list(0))
    }
    
    ## Signal Noise Weights ##
    if(TRUE){
      
      ## Weighting Matrix (Physical Space) ##
      grid = expand.grid(x = 0:(n-1),y = 0:(n-1))/n
      centre_x = (5)/12; centre_y = (5)/12; sigma_x = 0.4; sigma_y = 0.4; rho = 0
      weights = with(grid,sech2d(x = x,y = y, centre_x = centre_x,
                                 centre_y = centre_y, sigma_x = sigma_x,
                                 sigma_y=sigma_y, rho = rho))
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
      #spec_weights_mat = NULL
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
                                 param_bias_t_dep,seed = j,n_obs = n_obs,
                                 obs_loc = obs_loc, alpha_sim = TRUE,
                                 x_sim = TRUE, obs_noise_sim = TRUE,
                                 obs_bias_sim = TRUE, y_sim = FALSE, 
                                 weights = weights)
    }
    
    ## Log scale & log indices ##
    if(TRUE){
      log_scale = FALSE
      log_indices = c(1,2,3,4,5,8+(1:param_obs0$m))
    }
    
    ## Parameters to estimate ##
    if(TRUE){
      grad_params = 7:9
    }
    
    ## Observations to be optimised ##
    if(TRUE){
      grad_obs = 1:1
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
        step_size_RML[i,7:8] = 100*step_size_vec[7:8]*i^(-0.6) #unique to this sim
        step_size_RML[i,7:8] = 500*step_size_vec[7:8]*i^(-1.00) #unique to this sim
        step_size_RML[i,9] = 300*step_size_vec[9]*i^(-0.55) #unique to this sim
      }
      
      if(FALSE){
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
      step_size_vec = c(.5,.5)
      step_size_vec = c(0.2,0.2)
      #step_size_vec = c(0.01,0.01)
      
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
      W_coords = list(c(10,11),c(6,8),c(4,4),c(9,6),c(1,7),c(7,10),
                      c(1,1),c(3,10))
      
      ## Weighting indices ##
      W_indices = lapply(W_coords, function(A) A[2]*n+A[1]+1)
      
      ## Weighting matrix (physical space) ##
      W = matrix(0,nrow=n^2,ncol=n^2)
      diag(W)[unlist(W_indices)] = rep(10,length(W_coords))
      
      ## Alternative: no weighting matrix ##
      if(TRUE){
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
      RML_plot = FALSE
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
      #plot_limits[7,] = c(0.2,1.0)
      #plot_limits[8,] = c(-1.0,-0.2)
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
      RML_plot_freq = t
    }
    
    ## Save RML plot? ##
    if(TRUE){
      save_RML_plot = FALSE
    }
    
    ## RML Plot Filename ##
    if(TRUE){
      RML_plot_filename = paste("RML_plot5_",j,".pdf",sep="")
    }
    
    ## ROSP Plots? ##
    if(TRUE){
      ROSP_plot = FALSE
      ROSP_plot2d = FALSE
    }
    
    ## ROSP Plot Point Frequency ##
    if(TRUE){
      ROSP_plot_point_freq = 10
    }
    
    ## ROSP Plot Frequency ##
    if(TRUE){
      ROSP_plot_freq = t
    }
    
    ## Save ROSP Plots? ##
    if(TRUE){
      save_ROSP_plot = FALSE
      save_ROSP_plot2d = FALSE
    }
    
    ## ROSP Plot Filenames ##
    if(TRUE){
      ROSP_plot_filename = paste("ROSP_plot5_",j,".pdf",sep="")
      ROSP_plot2d_filename = paste("ROSP_plot2d5_",j,".pdf",sep="")
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
    RML_ROSP_plot5 = RML_ROSP(spde_sim = spde_sim, n = n, n_obs = n_obs, t = t, 
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
                              spec_weights_mat = spec_weights_mat,
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
  saveRDS(RML_ROSP_plot5,paste("RML_ROSP_plot5_",j,sep=""))
  
  ## Add to matrices ##
  all_param_est_mux[,j] = RML_ROSP_plot5$param_est[,7]
  all_param_est_muy[,j] = RML_ROSP_plot5$param_est[,8]
  all_param_est_tau2[,j] = RML_ROSP_plot5$param_est[,9]
  all_obs_x[,j] = RML_ROSP_plot5$obs[[1]][,1]
  all_obs_y[,j] = RML_ROSP_plot5$obs[[1]][,2]
  
  ## Print Iterations ##
  print(paste(length(init_loc)-j," iterations remaining"))
}

saveRDS(all_param_est_mux,"all_sims_mux")
saveRDS(all_param_est_muy,"all_sims_muy")
saveRDS(all_param_est_tau2,"all_sims_tau2")
saveRDS(all_obs_x,"all_sims_obs_x")
saveRDS(all_obs_y,"all_sims_obs_y")

if(FALSE){
all_param_est_mux = readRDS("all_sims_mux")
all_param_est_muy = readRDS("all_sims_muy")
all_param_est_tau2 = readRDS("all_sims_tau2")
all_obs_x = readRDS("all_sims_obs_x")
all_obs_y = readRDS("all_sims_obs_y")
}

##################

#########################
#### Parameter Plots ####
#########################

## mu_x ##
#pdf("mu_x.pdf",height = 7, width = 9)
require(tikzDevice)
tikz('mu_x.tex', standAlone=TRUE, height = 7, width = 10)
par(mfrow=c(1,1),mar=c(2.7,2.7,0.3,0.3),mgp=c(1.4,0.6,0),family="LM Roman 10",cex=1.5)
ymid=apply(all_param_est_mux,1,mean)
ylow=apply(all_param_est_mux,1,mean)-apply(all_param_est_mux,1,sd)
yhigh=apply(all_param_est_mux,1,mean)+apply(all_param_est_mux,1,sd)
plot(ymid,type = 'n', ylim = c(0.05, 0.43),xlab="$t$",ylab="$\\theta(t)$")
lines(ylow, col = 'red',lwd=3);lines(yhigh, col = 'red',lwd=3)
polygon(c(1:(t+1),rev(1:(t+1))), c(yhigh, rev(ylow)), 
        col = rgb(1, 0, 0,0.3), border = NA)
lines(ymid,col="red",lwd=3)
abline(h=param_sig[7],lty=2,lwd=3)
dev.off()
tools::texi2dvi('mu_x.tex',pdf=T)

## mu_y ##
#cairo_pdf("mu_y.pdf",height = 7, width = 9)
tikz('mu_y.tex', standAlone=TRUE, height = 7, width = 10)
par(mfrow=c(1,1),mar=c(2.7,2.7,0.3,0.3),mgp=c(1.4,0.6,0),family="LM Roman 10",cex=1.5)
ymid=apply(all_param_est_muy,1,mean)
ylow=apply(all_param_est_muy,1,mean)-apply(all_param_est_muy,1,sd)
yhigh=apply(all_param_est_muy,1,mean)+apply(all_param_est_muy,1,sd)
plot(ymid,type = 'n', ylim = c(-0.43, -0.06),xlab=expression(t),ylab='$\\theta(t)$')
lines(ylow, col = 'blue',lwd=3);lines(yhigh, col = 'blue',lwd=3)
polygon(c(1:(t+1),rev(1:(t+1))), c(yhigh, rev(ylow)), 
        col = rgb(0, 0, 1,0.3), border = NA)
lines(ymid,col="blue",lwd=3)
abline(h=param_sig[8],lty=2,lwd=3)
dev.off()
tools::texi2dvi('mu_y.tex',pdf=T)

## tau2 ##
#pdf("tau2.pdf",height = 7, width = 9)
tikz('tau2_1.tex', standAlone=TRUE, height = 7, width = 10)
par(mfrow=c(1,1),mar=c(2.7,2.7,0.3,0.3),mgp=c(1.4,0.6,0),family="LM Roman 10",cex=1.5)
ymid=apply(all_param_est_tau2,1,mean)
ylow=apply(all_param_est_tau2,1,mean)-apply(all_param_est_tau2,1,sd)
yhigh=apply(all_param_est_tau2,1,mean)+apply(all_param_est_tau2,1,sd)
plot(ymid,type = 'n', ylim = c(0.0, 0.50),xlab="$t$",ylab="$\\theta$",lwd=3)
lines(ylow, col = 'green',lwd=3);lines(yhigh, col = 'green',lwd=3)
polygon(c(1:(t+1),rev(1:(t+1))), c(yhigh, rev(ylow)), 
        col = rgb(0, 1, 0,0.3), border = NA)
lines(ymid,col="green",lwd=3)
abline(h=param_obs$tau2_vals[[1]],lty=2)
dev.off()
tools::texi2dvi('tau2_1.tex',pdf=T)


## all parameters ##
require(tikzDevice)
tikz('all_params.tex', standAlone=TRUE, height = 7, width = 14)
par(mfrow=c(1,1),mar=c(2.7,2.7,0.3,0.3),mgp=c(1.6,0.6,0),family="LM Roman 10",cex=3.5)
ymid=apply(all_param_est_mux,1,mean)
ylow=apply(all_param_est_mux,1,mean)-apply(all_param_est_mux,1,sd)
yhigh=apply(all_param_est_mux,1,mean)+apply(all_param_est_mux,1,sd)
plot(ymid,type = 'n', ylim = c(-0.43, 0.60),xlab="$t$",ylab="$\\theta(t)$")
lines(ylow, col = 'red',lwd=3);lines(yhigh, col = 'red',lwd=3)
polygon(c(1:(t+1),rev(1:(t+1))), c(yhigh, rev(ylow)), 
        col = rgb(1, 0, 0,0.3), border = NA)
lines(ymid,col="red",lwd=3)
abline(h=param_sig[7],lty=2,lwd=3)

ymid=apply(all_param_est_muy,1,mean)
ylow=apply(all_param_est_muy,1,mean)-apply(all_param_est_muy,1,sd)
yhigh=apply(all_param_est_muy,1,mean)+apply(all_param_est_muy,1,sd)
lines(ylow, col = 'blue',lwd=3);lines(yhigh, col = 'blue',lwd=3)
polygon(c(1:(t+1),rev(1:(t+1))), c(yhigh, rev(ylow)), 
        col = rgb(0, 0, 1,0.3), border = NA)
lines(ymid,col="blue",lwd=3)
abline(h = param_sig[8],lty=2,lwd=3)

ymid=apply(all_param_est_tau2,1,mean)
ylow=apply(all_param_est_tau2,1,mean)-apply(all_param_est_tau2,1,sd)
yhigh=apply(all_param_est_tau2,1,mean)+apply(all_param_est_tau2,1,sd)
lines(ylow, col = 'green',lwd=3);lines(yhigh, col = 'green',lwd=3)
polygon(c(1:(t+1),rev(1:(t+1))), c(yhigh, rev(ylow)), 
        col = rgb(0, 1, 0,0.3), border = NA)
lines(ymid,col="green",lwd=3)
abline(h=param_obs$tau2_vals[[1]],lty=2,lwd=3)

legend(x=600,y=0.6,ncol=3,legend=c("$\\mu_{x}(t)$","$\\mu_{y}(t)$","$\\tau^2(t)$"),lwd=3,col=c("red","blue","green"))
dev.off()
tools::texi2dvi('all_params.tex',pdf=T)

#########################

###########################
#### Observation Plots ####
###########################

## Spatial Weighting Function ##
require(tikzDevice)
tikz('signal_weights.tex', standAlone=TRUE, height = 7, width = 7)
par(mar=c(2.8,2.6,0,0),mfrow=c(1,1),mgp=c(1.2,1.0,0),cex=1.5)
library(lattice)
theme.novpadding <- list(
  layout.heights = list(
    top.padding = -4.0,
    main.key.padding = 0,
    key.axis.padding = 0,
    axis.xlab.padding = 0,
    xlab.key.padding = 0,
    key.sub.padding = 0,
    bottom.padding = 2.0
  ),
  layout.widths = list(
    left.padding = 0.5,
    key.ylab.padding = 0,
    ylab.axis.padding = 0,
    axis.key.padding = -0.5,
    right.padding = 0.5
  )
)
wireframe(weights~x*y, data = grid, cuts = 100, main="", 
          col.regions = heat.colors(1000)[1:length(heat.colors(1000))],
          drape = TRUE, zlab=list("",cex=1.5),
          xlab=list("$x$",cex=1.5,hjust="1",vjust="-5"),
          ylab=list("$y$",cex=1.5,hjust="1",vjust="-5"), 
          light.source = c(10,10,10),xlim=c(0,1),ylim=c(0,1),
          aspect = c(1,1), screen = list(z = 45, x = -60), col = F, 
          zoom = 0.95, scales = list(arrows = F, tck = 1.0, 
                                  distance =c(2.0, 2.0, 1.5),
                                  cex=1.5, at = c(seq(0,1,0.2))),
          par.settings = theme.novpadding,
          colorkey=list(labels = list(cex=1.5)))
dev.off()
tools::texi2dvi('signal_weights.tex',pdf=T)




## 2D Plot ##
require(tikzDevice)
tikz('2d_plot.tex', standAlone=TRUE, height = 7, width = 7)

## Plot points ##
plot_points = c(seq(1,200,5),seq(201,300,10),seq(301,400,20),seq(401,500,50),
                seq(501,600,50),seq(601,1000,100),seq(1001,t,200))

## Plot parameters ##
par(mar=c(2.8,2.6,0.7,0.7),mfrow=c(1,1),mgp=c(1.5,0.5,0),cex=1.5)

## Initialise Plot ##
plot(-1,-1,xlim=c(0,n),ylim=c(0,n),xlab="$x$",xaxt="n",yaxt="n",
     ylab="$y$",xaxs="i",yaxs="i")

axis(side=1,at=seq(0,n,length.out=6),labels = c(0,0.2,0.4,0.6,0.8,1.0))
axis(side=2,at=seq(0,n,length.out=6),labels = c(0,0.2,0.4,0.6,0.8,1.0))

## Legend ##
legend(x="bottomright", inset = c(0.035,0.6), pt.cex = c(2,2,2),
       pch = c(21,21,21),
       legend=c("Initial Sensor Locations",
                "Fixed Sensor Locations", 
                "Final Sensor Locations"),
       pt.bg = c(brewer.pal(8,"Greens")[2],
                 brewer.pal(8,"Reds")[2],
                 brewer.pal(8,"Blues")[2]),
       col = c(brewer.pal(8,"Greens")[7],
               brewer.pal(8,"Reds")[7],
               brewer.pal(8,"Blues")[7]))

## Add grid ##
abline(h=12*seq(0,1,length.out=11),v=12*seq(0,1,length.out=11),col="gray90",lty=2)


counter = 1

## Loop Over All Initial Conditions ##
for (i in c(10,102,177,370)){
  
  ## Plot Initial Location ##
  points(all_obs_x[1,i],all_obs_y[1,i][[1]],xlim = c(0,n),ylim = c(0,n),
       cex = 2, pch = 21, col=brewer.pal(8,"Greens")[9-1],
       bg = brewer.pal(8,"Greens")[2])
  
  ## Plot Fixed Sensor Locations ##
  fixed_obs = list(c(2,2),c(2,6),c(2,10),
                   c(6,2),c(6,6),c(6,10),
                   c(10,2),c(10,6),c(10,10))
  
  for (k in 1:length(fixed_obs)){
    points(fixed_obs[[k]][1],fixed_obs[[k]][2],
           cex = 2, pch = 21, col = brewer.pal(8,"Reds")[7],
           bg = brewer.pal(8,"Reds")[2])
  }
      
  for (l in plot_points){
    points(all_obs_x[l,i],all_obs_y[l,i],
           col=brewer.pal(9,"Greens")[9-1],cex=.6,type="b",pch=1,lwd=1,lty=1)
  }
  
  ## Plot Final Location ##
  points(all_obs_x[t,i],all_obs_y[t,i],
         cex = 2, pch = 21, col = brewer.pal(8,"Blues")[7],
         bg = brewer.pal(8,"Blues")[2])
  
  counter = counter + 1
}

dev.off()
tools::texi2dvi('2d_plot.tex',pdf=T)




## Trial-Averaged Plots ##

## x coord ##
tikz('obs_x.tex', standAlone=TRUE, height = 3.5, width = 10)
par(mfrow=c(1,1),mar=c(2.7,2.7,0.3,0.3),mgp=c(1.4,0.6,0),family="LM Roman 10",cex=1.5)
ymid=1/12*apply(all_obs_x,1,mean)
ylow=1/12*(apply(all_obs_x,1,mean)-apply(all_obs_x,1,sd))
yhigh=1/12*(apply(all_obs_x,1,mean)+apply(all_obs_x,1,sd))
plot(ymid,type = 'n', ylim = c(0,1),xlab=expression(t),ylab='$o(t)$')
lines(ylow, col = 'red',lwd=3);lines(yhigh, col = 'red',lwd=3)
polygon(c(1:(t+1),rev(1:(t+1))), c(yhigh, rev(ylow)), 
        col = rgb(1, 0, 0, 0.3), border = NA)
lines(ymid,col="red",lwd=3)
abline(h=5,lty=2,lwd=3)
dev.off()
tools::texi2dvi('obs_x.tex',pdf=T)

## y coord ##
tikz('obs_y.tex', standAlone=TRUE, height = 3.5, width = 10)
par(mfrow=c(1,1),mar=c(2.7,2.7,0.3,0.3),mgp=c(1.4,0.6,0),family="LM Roman 10",cex=1.5)
ymid=1/12*apply(all_obs_y,1,mean)
ylow=1/12*(apply(all_obs_y,1,mean)-apply(all_obs_y,1,sd))
yhigh=1/12*(apply(all_obs_y,1,mean)+apply(all_obs_y,1,sd))
plot(ymid,type = 'n', ylim = c(0,1),xlab=expression(t),ylab='$o(t)$')
lines(ylow, col = 'blue',lwd=3);lines(yhigh, col = 'blue',lwd=3)
polygon(c(1:(t+1),rev(1:(t+1))), c(yhigh, rev(ylow)), 
        col = rgb(0, 0, 1,0.3), border = NA)
lines(ymid,col="blue",lwd=3)
abline(h=5,lty=2,lwd=3)
dev.off()
tools::texi2dvi('obs_y.tex',pdf=T)

## all coords ##
tikz('obs.tex', standAlone=TRUE, height = 7, width = 14)
par(mfrow=c(1,1),mar=c(2.7,2.7,0.3,0.3),mgp=c(1.4,0.6,0),family="LM Roman 10",cex=3.5)
ymid=1/12*apply(all_obs_x,1,mean)
ylow=1/12*(apply(all_obs_x,1,mean)-apply(all_obs_x,1,sd))
yhigh=1/12*(apply(all_obs_x,1,mean)+apply(all_obs_x,1,sd))
plot(ymid,type = 'n', ylim = c(0.18,0.8),xlab=expression(t),ylab='$o(t)$')
lines(ylow, col = 'red',lwd=3);lines(yhigh, col = 'red',lwd=3)
polygon(c(1:(t+1),rev(1:(t+1))), c(yhigh, rev(ylow)), 
        col = rgb(1, 0, 0, 0.3), border = NA)
lines(ymid,col="red",lwd=3)

ymid=1/12*apply(all_obs_y,1,mean)
ylow=1/12*(apply(all_obs_y,1,mean)-apply(all_obs_y,1,sd))
yhigh=1/12*(apply(all_obs_y,1,mean)+apply(all_obs_y,1,sd))
lines(ylow, col = 'blue',lwd=3);lines(yhigh, col = 'blue',lwd=3)
polygon(c(1:(t+1),rev(1:(t+1))), c(yhigh, rev(ylow)), 
        col = rgb(0, 0, 1,0.3), border = NA)
lines(ymid,col="blue",lwd=3)
abline(h=5/12,lty=2,lwd=3)

legend(x=1800,y=0.78,legend=c("$o_x(t)$    ","$o_y(t)$    "),lwd=3,col=c("red","blue"))

dev.off()
tools::texi2dvi('obs.tex',pdf=T)

###########################





