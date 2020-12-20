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

## Setup ##
if(TRUE){
  
  ## Prelims ##
  if(TRUE){
    n = 12; t = 100000; dt = 1; K = 21; nu = 1
  }
  
  ## RML & OSP Times ##
  if(TRUE){
    t_RML = 49; t_ROSP = 1
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
                               obs_bias_sim = TRUE, y_sim = FALSE)
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
    step_size_vec = c(.1,.1)
    
    ## Step size matrix ##
    step_size_mat = matrix(0,nrow=n_iterations*t,ncol=2)
    for (i in 1:(n_iterations*t)){
      step_size_mat[i,]=step_size_vec*i^(-0.9)
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
    RML_plot_freq = 50000
  }
  
  ## Save RML plot? ##
  if(TRUE){
    save_RML_plot = TRUE
  }
  
  ## RML Plot Filename ##
  if(TRUE){
    RML_plot_filename = "RML_plot2.pdf"
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
    ROSP_plot_freq = 500
  }
  
  ## Save ROSP Plots? ##
  if(TRUE){
    save_ROSP_plot = TRUE
    save_ROSP_plot2d = TRUE
  }
  
  ## ROSP Plot Filenames ##
  if(TRUE){
    ROSP_plot_filename = "ROSP_plot2.pdf"
    ROSP_plot2d_filename = "ROSP_plot2d2.pdf"
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
  RML_ROSP_plot2 = RML_ROSP(spde_sim = spde_sim, n = n, n_obs = n_obs, t = t, 
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
saveRDS(RML_ROSP_plot2,"RML_ROSP_plot2")

##################

##############
#### PLOT ####
##############

## RML ##
if(TRUE){
  setwd(fig_wd)
  #pdf("RML_plot2.pdf",width=7,height=7)
  tikz('RML_plot2.tex', standAlone=TRUE, height = 7, width = 7)
  par(mar=c(3,4.3,.6,2),mgp=c(2.7,1,0),family="LM Roman 10",cex=2)
  par(mfrow=c(3,3))
  plot_points = seq(1,t,100)
  if(FALSE){
    setwd("/Users/ls616/Google Drive/MPE CDT/PhD/Year 2/Application_Paper/Current/Draft5/Figures/Section4/RML_ROSP_plot2")
    RML_ROSP_plot2 = readRDS("RML_ROSP_plot2")
  }
  param_est = RML_ROSP_plot2$param_est
  param_sig_indices = 1:8
  param_sig_true = param_sig_t_dep
  param_obs_indices = 9
  param_obs_true = param_obs_t_dep
  param_sig_names = list("$\\hat{\\rho}_0$","$\\hat{\\sigma}^2$",
                         "$\\hat{\\zeta}$","$\\hat{\\rho}_1$",
                         "$\\hat{\\gamma}$","$\\hat{\\alpha}$",
                         "$\\hat{\\mu}_x$","$\\hat{\\mu}_y$")
  param_obs_names = list("$\\hat{\\tau}^2$")
  
  ## Signal Parameters ##
  for (k in param_sig_indices){
    
    y_lim = NULL
    y_lab = expression(theta)
    main = NULL
    
    if(!is.null(plot_limits)){
      y_lim = plot_limits[[k]]
    }
    
    if(!is.null(param_sig_names)){
      y_lab = param_sig_names[[k]]
      main = param_sig_names[[k]]
    }
    
    plot(plot_points, param_est[plot_points,k], cex=0.4,
         xlim = c(1,n_iterations*t+1), xlab = "",ylab = y_lab,
         ylim = y_lim, main = "", cex.lab = 1.8, cex.axis = 1.8)
    
    if (!is.null(param_sig_true)){
      if(length(param_sig_true$t_vals)==1){
        abline(h=param_sig_true$param_vals[[1]][k],col="red",lty=2)
      }
      else{
        for(l in 1:(length(param_sig_true$t_vals)-1)){
          segments(param_sig_true$t_vals[[l]],
                   param_sig_true$param_vals[[l]][k],
                   param_sig_true$t_vals[[l+1]],
                   param_sig_true$param_vals[[l]][k],
                   col="red",lty=2)
          segments(param_sig_true$t_vals[[l+1]],
                   param_sig_true$param_vals[[l]][k],
                   param_sig_true$t_vals[[l+1]],
                   param_sig_true$param_vals[[l+1]][k],
                   col="red",lty=2)
        }
        segments(param_sig_true$t_vals[[length(param_sig_true$t_vals)]],
                 param_sig_true$param_vals[[length(param_sig_true$t_vals)]][k],
                 n_iterations*t+1,
                 param_sig_true$param_vals[[length(param_sig_true$t_vals)]][k],
                 col="red",lty=2)
      }
    }
    
    if(mean_param){
      if(t_end_RML>(param_est_starts[[k]][1]+1)){
        for (m in 1:(which_start[[k]])){
          abline(v=param_est_starts[[k]][m],col="black",lty=2)
          if(length(param_sig_true$t_vals)==1){
            abline(h=mean_param_est[[k]][1],col="blue",lty=2)
          }
          else{
            segments(mean_param_est_plot_start[[k]][m],mean_param_est[[k]][m],
                     param_est_ends[[k]][m],mean_param_est[[k]][m],
                     col="blue",lty=2)
          }
        }
      }
    }
  }
  
  ## Observation Parameters ##
  for (k in param_obs_indices){
    if(k==9){
      
      y_lim = NULL
      y_lab = expression(theta)
      main = NULL
      
      if(!is.null(plot_limits)){
        y_lim = plot_limits[[9]]
      }
      
      if(!is.null(param_obs_names)){
        y_lab = param_obs_names[[k-8]]
        main = param_obs_names[[k-8]]
      }
      
      plot(plot_points, param_est[plot_points,k], cex=0.4,
           xlim = c(1,n_iterations*t+1) ,xlab="", ylab=y_lab,
           ylim = y_lim, main = "",cex.lab = 1.8, cex.axis = 1.8)
    }
    
    if(k>9){
      points(plot_points,param_est[plot_points,k],pch=19,cex=0.4)
    }
    
    if (!is.null(param_obs_true)){
      if(length(param_obs_true$t_vals)==1){
        abline(h=param_obs_true$param_vals[[1]]$tau2_vals[[k-8]],
               col="red",lty=2)
      }
      else{
        for(l in 1:(length(param_obs_true$t_vals)-1)){
          segments(param_obs_true$t_vals[[l]],
                   param_obs_true$param_vals[[l]]$tau2_vals[[k-8]],
                   param_obs_true$t_vals[[l+1]],
                   param_obs_true$param_vals[[l]]$tau2_vals[[k-8]],
                   col="red",lty=2)
          segments(param_obs_true$t_vals[[l+1]],
                   param_obs_true$param_vals[[l]]$tau2_vals[[k-8]],
                   param_obs_true$t_vals[[l+1]],
                   param_obs_true$param_vals[[l+1]]$tau2_vals[[k-8]],
                   col="red",lty=2)
        }
        segments(param_obs_true$t_vals[[length(param_obs_true$t_vals)]],
                 param_obs_true$param_vals[[length(param_obs_true$t_vals)]]$tau2_vals[[k-8]],
                 n_iterations*t+1,
                 param_obs_true$param_vals[[length(param_obs_true$t_vals)]]$tau2_vals[[k-8]],
                 col="red",lty=2)
      }
    }
  }
  
  dev.off()
  tools::texi2dvi('RML_plot2.tex',pdf=T)
}

## ROSP ##
if(TRUE){
  
  setwd(fig_wd)
  #pdf("ROSP_plot2.pdf",width=7,height=7)
  tikz('ROSP_plot2.tex', standAlone=TRUE, height = 7, width = 7)
  par(mar=c(3,3,2,2),mgp=c(2.7,1,0),family="LM Roman 10",cex=2)
  par(mfrow=c(2,4))
  plot_points = c(seq(1,10000,20),seq(10001,t,500))
  if(FALSE){
    setwd("/Users/ls616/Google Drive/MPE CDT/PhD/Year 2/Application_Paper/Current/Draft5/Figures/Section4/RML_ROSP_plot2")
    RML_ROSP_plot2 = readRDS("RML_ROSP_plot2")
  }
  obs = RML_ROSP_plot2$obs
  
  
  ## Plot ##
  for(k in 1:n_obs){
    plot(plot_points,obs[[k]][plot_points,1],cex=.4,ylim=c(0,n),
         xlim=c(0,n_iterations*t),xlab="",
         ylab="$\\hat{\\mathbf{o}}$",col=brewer.pal(8,"Greys")[6],
         main=paste("Sensor",k),cex.axis=2.0,cex.main=2.0,cex.lab=2.0)
    points(plot_points,obs[[k]][plot_points,2],cex=.4,
           col=brewer.pal(8,"Greys")[8])
    legend(x=45000,y=3.2,
           legend=c("$\\hat{\\mathbf{o}}_x$","$\\hat{\\mathbf{o}}_x$"),
           col=brewer.pal(8,"Greys")[c(6,8)],pch=19,ncol=1,cex=1.8)
    
    ## Add dashed lines for weighted coordinates ##
    if(!is.null(W_coords)){
      abline(h=W_coords[[k]][1],lty=2,col=brewer.pal(8,"Reds")[3])
      abline(h=W_coords[[k]][2],lty=2,col=brewer.pal(8,"Reds")[7])
    }
  }
  
  dev.off()
  tools::texi2dvi('ROSP_plot2.tex',pdf=T)
}

##############

#####################
#### RML_no_ROSP ####
#####################

## Run 'RML' by setting 'ROSP' step-sizes to 0 in RML_ROSP ##
## Speed up by setting 'grad_obs = 1' ##

## Setup ##
if(TRUE){
  
  ## Parameters to estimate ##
  if(TRUE){
    grad_params = 1:n_param
  }
  
  ## Observations to be optimised ##
  if(TRUE){
    grad_obs = 1
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
    step_size_vec = c(0,0)
    
    ## Step size matrix ##
    step_size_mat = matrix(0,nrow=n_iterations*t,ncol=2)
    for (i in 1:(n_iterations*t)){
      step_size_mat[i,]=step_size_vec*i^(-0.505)
    }
    
    ## Step size matrices for all sensors ##
    step_size_ROSP = rep(list(step_size_mat),n_obs)
  }
  
  ## RML Plot? ##
  if(TRUE){
    RML_plot = TRUE
  }
  
  ## Save RML plot? ##
  if(TRUE){
    save_RML_plot = TRUE
  }
  
  ## RML Plot Filename ##
  if(TRUE){
    RML_plot_filename = "RML_no_ROSP_plot2.pdf"
  }
  
  ## ROSP Plots? ##
  if(TRUE){
    ROSP_plot = FALSE
    ROSP_plot2d = FALSE
  }
  
  ## Save ROSP Plots? ##
  if(TRUE){
    save_ROSP_plot = FALSE
    save_ROSP_plot2d = FALSE
  }
  
}

## RML_ROSP ##
if(TRUE){
  RML_no_ROSP_plot2 = RML_ROSP(spde_sim = spde_sim, n = n, n_obs = n_obs, t = t, 
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
saveRDS(RML_no_ROSP_plot2,"RML_no_ROSP_plot2")

#####################

#####################
#### ROSP_no_RML ####
#####################

## Run 'ROSP' by setting 'RML' step-sizes to 0 in RML_ROSP ##
## Speed up by setting 'grad_param = 1' ##

## Setup ##
if(TRUE){
  
  ## Parameters to estimate ##
  if(TRUE){
    grad_params = 1
  }
  
  ## Observations to be optimised ##
  if(TRUE){
    grad_obs = 1:n_obs
  }
  
  ## RML Step Sizes ##
  if(TRUE){
    
    ## Initial step-size vector ##
    step_size_vec = rep(0,9)
    
    ## Step size matrix ##
    step_size_RML = matrix(0,nrow=n_iterations*t,ncol=n_param)
    for (i in 1:dim(step_size_RML)[1]){
      step_size_RML[i,] = step_size_vec
    }
    
  }
  
  ## ROSP Step sizes ##
  if(TRUE){
    
    ## Initial step size vector ##
    step_size_vec = c(10,10)
    step_size_vec = c(1,1)
    step_size_vec = c(.1,.1)
    
    ## Step size matrix ##
    step_size_mat = matrix(0,nrow=n_iterations*t,ncol=2)
    for (i in 1:(n_iterations*t)){
      step_size_mat[i,]=step_size_vec*i^(-0.9)
    }
    
    ## Step size matrices for all sensors ##
    step_size_ROSP = rep(list(step_size_mat),n_obs)
  }
  
  ## RML Plot? ##
  if(TRUE){
    RML_plot = FALSE
  }
  
  ## Save RML plot? ##
  if(TRUE){
    save_RML_plot = FALSE
  }
  
  ## ROSP Plots? ##
  if(TRUE){
    ROSP_plot = TRUE
    ROSP_plot2d = FALSE
  }
  
  ## Save ROSP Plots? ##
  if(TRUE){
    save_ROSP_plot = TRUE
    save_ROSP_plot2d = TRUE
  }
  
  ## ROSP Plot Filenames ##
  if(TRUE){
    ROSP_plot_filename = "ROSP_no_RML_plot2.pdf"
    ROSP_plot2d_filename = "ROSP_no_RML_plot2d2.pdf"
  }
  
}

## RML_ROSP ##
if(TRUE){
  ROSP_no_RML_plot2 = RML_ROSP(spde_sim = spde_sim, n = n, n_obs = n_obs, t = t, 
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
saveRDS(ROSP_no_RML_plot2,"ROSP_no_RML_plot2")

#####################

##################
#### PLOT MSE ####
##################

## Simulate spde on n x n grid ##
spde_sim_all = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                               param_bias_t_dep,seed = 1,alpha_sim = TRUE,
                               x_sim = TRUE, obs_noise_sim = TRUE,
                               obs_bias_sim = TRUE, y_sim = FALSE)

## Simulate spde at initial locations ##
spde_sim_init_loc = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                                    param_bias_t_dep,seed = 1,n_obs = n_obs,
                                    obs_loc = obs_loc,alpha_sim = TRUE,
                                    x_sim = TRUE, obs_noise_sim = TRUE,
                                    obs_bias_sim = TRUE, y_sim = TRUE)

## Simulate spde at target locations ##
spde_sim_target_loc = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                                      param_bias_t_dep,seed = 1,n_obs = n_obs,
                                      obs_loc = target_loc,alpha_sim = TRUE,
                                      x_sim = TRUE, obs_noise_sim = TRUE,
                                      obs_bias_sim = TRUE, y_sim = TRUE)

## Kalman Filter (true parameters, target sensor locations) ##
mu0 = rep(0,K)
Sigma0 = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,
             spec_weights_mat = spec_weights_mat)
Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs,W_coords)
XBeta = XBeta_f(param_bias,n,t,n_obs)
R = R_f(param_obs,n,n_obs)
G = G_f(param_sig,wave,cosine_indices,dt=1,n_cos=n_cos)
Q = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,
        spec_weights_mat = spec_weights_mat)
kf_true_param_target_loc = KalmanFilter(y = spde_sim_target_loc$y,Phi = Phi, R = R, G = G, Q = Q)

## Kalman Filter (true parameters, initial sensor locations) ##
mu0 = rep(0,K)
Sigma0 = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,
             spec_weights_mat = spec_weights_mat)
Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs,obs0)
XBeta = XBeta_f(param_bias,n,t,n_obs)
R = R_f(param_obs,n,n_obs)
G = G_f(param_sig,wave,cosine_indices,dt=1,n_cos=n_cos)
Q = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,
        spec_weights_mat = spec_weights_mat)
kf_true_param_init_loc = KalmanFilter(y = spde_sim_init_loc$y,Phi = Phi, R = R, G = G, Q = Q)

## Kalman Filter (initial parameters, target sensor locations) ##
mu0 = rep(0,K)
Sigma0 = Q_f(param_sig0,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,
             spec_weights_mat = spec_weights_mat)
Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs,W_coords)
XBeta = XBeta_f(param_bias0,n,t,n_obs)
R = R_f(param_obs0,n,n_obs)
G = G_f(param_sig0,wave,cosine_indices,dt=1,n_cos=n_cos)
Q = Q_f(param_sig0,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,
        spec_weights_mat = spec_weights_mat)
kf_init_param_target_loc = KalmanFilter(y = spde_sim_target_loc$y,Phi = Phi, R = R, G = G, Q = Q)

## Kalman Filter (initial parameters, initial sensor locations) ##
mu0 = rep(0,K)
Sigma0 = Q_f(param_sig0,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,
             spec_weights_mat = spec_weights_mat)
Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs,obs0)
XBeta = XBeta_f(param_bias0,n,t,n_obs)
R = R_f(param_obs0,n,n_obs)
G = G_f(param_sig0,wave,cosine_indices,dt=1,n_cos=n_cos)
Q = Q_f(param_sig0,wave,n,n_cos=n_cos,nu=1,dt=1,norm=TRUE,
        spec_weights_mat = spec_weights_mat)
kf_init_param_init_loc = KalmanFilter(y = spde_sim_init_loc$y,Phi = Phi, R = R, G = G, Q = Q)

## FFT matrix ##
fft_matrix = FFT_real_matrix(wave,cosine_indices,n_cos,n)

## Compute MSE ##
mean_square_errors_RML_ROSP = rep(0,t)
mean_square_errors_RML_no_ROSP = rep(0,t)
mean_square_errors_ROSP_no_RML = rep(0,t)
mean_square_errors_true_param_init_loc = rep(0,t)
mean_square_errors_true_param_target_loc = rep(0,t)
mean_square_errors_init_param_init_loc = rep(0,t)
mean_square_errors_init_param_target_loc = rep(0,t)

for (i in 1:t){
  mean_square_errors_RML_ROSP[i] = sum((as.vector(fft_matrix%*%RML_ROSP_plot2$kf_means[i,])
                                        -spde_sim_all$x[i,])^2)
  mean_square_errors_RML_no_ROSP[i] = sum((as.vector(fft_matrix%*%RML_no_ROSP_plot2$kf_means[i,])
                                           -spde_sim_all$x[i,])^2)
  mean_square_errors_ROSP_no_RML[i] = sum((as.vector(fft_matrix%*%ROSP_no_RML_plot2$kf_means[i,])
                                           -spde_sim_all$x[i,])^2)
  mean_square_errors_true_param_init_loc[i] = sum((as.vector(fft_matrix%*%kf_true_param_init_loc$m_i_i[i+1,])
                                                   -spde_sim_all$x[i,])^2)
  mean_square_errors_true_param_target_loc[i] = sum((as.vector(fft_matrix%*%kf_true_param_target_loc$m_i_i[i+1,])
                                                     -spde_sim_all$x[i,])^2)
  mean_square_errors_init_param_init_loc[i] = sum((as.vector(fft_matrix%*%kf_init_param_init_loc$m_i_i[i+1,])
                                                   -spde_sim_all$x[i,])^2)
  mean_square_errors_init_param_target_loc[i] = sum((as.vector(fft_matrix%*%kf_init_param_target_loc$m_i_i[i+1,])
                                                     -spde_sim_all$x[i,])^2)
  
}

## Plot MSE ##
library("zoo")
library("plotrix")
setwd(fig_wd)
pdf("MSE_plot2.pdf",width=13,height=4)
par(mfrow=c(1,1),mar=c(2.1,2.3,.8,.3),mgp = c(.9,.5,0),xpd=FALSE,family="LM Roman 10",cex=1.5)
plot_t_start = 150
plot_t_end = 100000
plot_t_freq = 100
plot(x = seq(plot_t_start,plot_t_end,plot_t_freq),y=rollmean(mean_square_errors_true_param_target_loc,1000,align="left")[seq(plot_t_start,plot_t_end,plot_t_freq)],lwd=2,col="blue",ylim = c(0.98,1.48),xlim = c(plot_t_start,plot_t_end+41000),xlab="t",type="l",xaxt="n",ylab="")
xtick = seq(0,100000,by = 20000)
axis(side = 1,at=xtick)
title(mgp=c(1.3, .5, 0),ylab="MSE")
lines(x = seq(plot_t_start,plot_t_end,plot_t_freq),y=rollmean(mean_square_errors_true_param_init_loc,1000,align="left")[seq(plot_t_start,plot_t_end,plot_t_freq)],lwd=2,col="red",type="l")
lines(x = seq(plot_t_start,plot_t_end,plot_t_freq),y=rollmean(mean_square_errors_init_param_init_loc,1000,align="left")[seq(plot_t_start,plot_t_end,plot_t_freq)],lwd=2,col="green",type="l")
lines(x = seq(plot_t_start,plot_t_end,plot_t_freq),y=rollmean(mean_square_errors_init_param_target_loc,1000,align="left")[seq(plot_t_start,plot_t_end,plot_t_freq)],lwd=2,col="orange",type="l")
lines(x = seq(plot_t_start,plot_t_end,plot_t_freq),y=rollmean(mean_square_errors_RML_ROSP,1000,align="left")[seq(plot_t_start,plot_t_end,plot_t_freq)],lwd=2,col="black",type="l")
lines(x = seq(plot_t_start,plot_t_end,plot_t_freq),y=rollmean(mean_square_errors_RML_no_ROSP,1000,align="left")[seq(plot_t_start,plot_t_end,plot_t_freq)],lwd=2,col="grey",type="l")
#lines(x = seq(plot_t_start,plot_t_end,plot_t_freq),y=rollmean(mean_square_errors_ROSP_no_RML,1000,align="left")[seq(plot_t_start,plot_t_end,plot_t_freq)],lwd=2,col="grey",type="l")
ablineclip(h = mean(mean_square_errors_true_param_target_loc), x2 = 100000,lty=2, lwd = 2)
ablineclip(h = mean(mean_square_errors_true_param_init_loc), x2 = 100000, lty=2,col="grey", lwd = 2)
#abline(v = c(20000,40000,60000,80000,100000),lty=2,col="lightgrey",lwd=0.5)
#abline(h = c(1.0,1.1,1.2,1.3,1.4,1.5),lty=2,col="lightgrey",lwd=0.5)

par(xpd=TRUE)
legend(x = 104000,y=1.05+0.03,legend = expression(paste(theta==theta(t),", ",bold(o)==bold(o)(t))) ,col="black",lwd=3,box.lty=0,bg="transparent") #legend = expression(paste(theta==theta[t]," , ",bold(o)==bold(o)[t], " (Learned Parameters, Learned Locations)"))
legend(x = 104000,y=1.11+0.03,legend = expression(paste(theta==theta[textstyle("*")],",   ",bold(o)==bold(o)[textstyle("*")])),col="blue",lwd=3,box.lty=0,bg="transparent") #legend = expression(paste(theta==theta["*"]," , ",bold(o)==bold(o)["*"], " (True Parameters, Target Locations)"))
legend(x = 104000,y=1.17+0.03,legend = expression(paste(theta==theta[0 ],",   ",bold(o)==bold(o)[textstyle("*")])), col="orange",lwd=3,box.lty=0,bg="transparent") #legend = expression(paste(theta==theta[0]," , ",bold(o)==bold(o)["*"], " (Initial Parameters, Target Locations)"))
legend(x = 104000,y=1.293+0.05,legend = expression(paste(theta==theta(t),", ",bold(o)==bold(o)[0])), col="grey",lwd=3,box.lty=0,bg="transparent") #legend = expression(paste(theta==theta[t]," , ",bold(o)==bold(o)[0], " (Learned Parameters, Initial Locations)"))
legend(x = 104000,y=1.353+0.05,legend = expression(paste(theta==theta[textstyle("*")],",   ",bold(o)==bold(o)[0])),col="red",lwd=3,box.lty=0,bg="transparent") #legend = expression(paste(theta==theta["*"]," , ",bold(o)==bold(o)[0], " (True Parameters, Initial Locations)"))
legend(x = 104000,y=1.413+0.05,legend = expression(paste(theta==theta[0 ],",   ",bold(o)==bold(o)[0])),col="green",lwd=3,box.lty=0,bg="transparent") #legend = expression(paste(theta==theta[0]," , ",bold(o)==bold(o)[0], " (Initial Parameters, Initial Locations)"))
dev.off()

## Additional Plots (start and end) ##
if(FALSE){
  par(mfrow=c(1,2))
  plot(rollmean(mean_square_errors_true_param_target_loc,100,align="right")[seq(1,1000,1)],col="blue",ylim = c(0.9,1.9),xlab="t",ylab="MSE",type="l",lwd=2)
  lines(rollmean(mean_square_errors_true_param_init_loc,100,align="right")[seq(1,1000,1)],col="red",lwd=2)
  lines(rollmean(mean_square_errors_init_param_init_loc,100,align="right")[seq(1,1000,1)],col="green",lwd=2)
  lines(rollmean(mean_square_errors_init_param_target_loc,100,align="right")[seq(1,1000,1)],col="orange",lwd=2)
  lines(rollmean(mean_square_errors_RML_ROSP,100,align="right")[seq(1,1000,1)],col="black",lwd=2)
  lines(rollmean(mean_square_errors_RML_no_ROSP,100,align="right")[seq(1,1000,1)],col="grey",lwd=2)
  lines(rollmean(mean_square_errors_ROSP_no_RML,100,align="right")[seq(1,1000,1)],col="grey",lwd=2)
  abline(h=mean(mean_square_errors_true_param_target_loc,align="right"),lty=2)
  
  plot(rollmean(mean_square_errors_true_param_target_loc,100,align="right")[seq(98901,99901,1)],col="blue",ylim = c(0.9,1.9),xlab="t",ylab="MSE",type="l",lwd=2)
  lines(rollmean(mean_square_errors_true_param_init_loc,100)[seq(98901,99901,1)],col="red",lwd=2)
  lines(rollmean(mean_square_errors_init_param_init_loc,100)[seq(98901,99901,1)],col="green",lwd=2)
  lines(rollmean(mean_square_errors_init_param_target_loc,100)[seq(98901,99901,1)],col="orange",lwd=2)
  lines(rollmean(mean_square_errors_RML_ROSP,100)[seq(98901,99901,1)],col="black",lwd=2)
  lines(rollmean(mean_square_errors_RML_no_ROSP,100,align="right")[seq(1,1000,1)],col="grey",lwd=2)
  lines(rollmean(mean_square_errors_ROSP_no_RML,100,align="right")[seq(1,1000,1)],col="grey",lwd=2)
  abline(h=mean(mean_square_errors_true_param_target_loc),lty=2)
}

##################

#################
#### PLOT KF ####
#################

## All Time (multiple fourier coefficients) ##
setwd(fig_wd)
cairo_pdf("KF_plot2_all_time_mult_coeff.pdf",width=13,height=9)
par(mfrow=c(2,2),mar=c(2.8,1.8,.5,.8),mgp = c(1.5,.5,0),xpd=FALSE,family="LM Roman 10",cex=1.5)
kf_plot_start = 1; kf_plot_end = 100000; kf_plot_freq = 50
alpha_to_plot = c(1,5,12,17)
for(i in alpha_to_plot){
  plot(seq(kf_plot_start,kf_plot_end,kf_plot_freq),kf_true_param_target_loc$m_i_i[seq(kf_plot_start+1,kf_plot_end+1,kf_plot_freq),i],type="l",xlab="t",ylab="",xaxt="n")#,yaxt="n")#ylab=bquote(alpha[.(i)](t)))
  xtick = seq(0,kf_plot_end,by = kf_plot_end/4)
  axis(side = 1,at=xtick)
  lines(seq(kf_plot_start,kf_plot_end,kf_plot_freq),RML_ROSP_plot2$kf_means[seq(kf_plot_start,kf_plot_end,kf_plot_freq),i],col="red")
}
dev.off()

## All Time (one fourier coefficient) ##
setwd(fig_wd)
cairo_pdf("KF_plot2_all_time_one_coeff.pdf",width=13,height=4)
par(mfrow=c(1,1),mar=c(2.1,2,.2,.4),mgp =c(1,.5,0),xpd=FALSE,family="LM Roman 10",cex=1.5)
kf_plot_start = 1; kf_plot_end = 100000; kf_plot_freq = 50
alpha_to_plot_new = alpha_to_plot[4]
plot(seq(kf_plot_start,kf_plot_end,kf_plot_freq),kf_true_param_target_loc$m_i_i[seq(kf_plot_start+1,kf_plot_end+1,kf_plot_freq),alpha_to_plot_new],type="l",xlab="t",ylab="",xaxt="n")#,yaxt="n")#ylab=bquote(alpha[.(i)](t)))
xtick = seq(0,kf_plot_end,by = kf_plot_end/5)
axis(side = 1,at=xtick)
lines(seq(kf_plot_start,kf_plot_end,kf_plot_freq),RML_ROSP_plot2$kf_means[seq(kf_plot_start,kf_plot_end,kf_plot_freq),alpha_to_plot_new],col="red")
dev.off()

## Start (one fourier coefficient) ##
setwd(fig_wd)
cairo_pdf("KF_plot2_start_one_coeff.pdf",width=6.5,height=4)
par(mfrow=c(1,1),mar=c(2.2,1.5,.5,.6),mgp =c(1.2,.5,0),xpd=FALSE,family="LM Roman 10",cex=1.5)
kf_plot2_start = 1; kf_plot2_end = 200; kf_plot2_freq = 1
alpha_to_plot2 = alpha_to_plot[4]
plot(seq(kf_plot2_start,kf_plot2_end,kf_plot2_freq),kf_true_param_target_loc$m_i_i[seq(kf_plot2_start+1,kf_plot2_end+1,kf_plot2_freq),alpha_to_plot2],type="l",xlab="t",ylab="")#ylab=bquote(alpha[.(i)](t)))
lines(seq(kf_plot2_start,kf_plot2_end,kf_plot2_freq),RML_ROSP_plot2$kf_means[seq(kf_plot2_start,kf_plot2_end,kf_plot2_freq),alpha_to_plot2],col="red")
dev.off()

## End (one fourier coefficient) ##
setwd(fig_wd)
cairo_pdf("KF_plot2_end_one_coeff.pdf",width=6.5,height=4)
par(mfrow=c(1,1),mar=c(2.2,1.5,.5,.6),mgp =c(1.2,.5,0),xpd=FALSE,family="LM Roman 10",cex=1.5)
kf_plot3_start = 99800
kf_plot3_end = 100000
kf_plot3_freq = kf_plot2_freq
alpha_to_plot3 = alpha_to_plot2
plot(seq(kf_plot3_start,kf_plot3_end,kf_plot3_freq),kf_true_param_target_loc$m_i_i[seq(kf_plot3_start+1,kf_plot3_end+1,kf_plot3_freq),alpha_to_plot3],type="l",xlab="t",ylab="")#ylab=bquote(alpha[.(i)](t)))
lines(seq(kf_plot3_start,kf_plot3_end,kf_plot3_freq),RML_ROSP_plot2$kf_means[seq(kf_plot3_start,kf_plot3_end,kf_plot3_freq),alpha_to_plot3],col="red")
dev.off()


## Start & End (one fourier coefficient) ##
setwd(fig_wd)
cairo_pdf("KF_plot2_start_end_one_coeff.pdf",width=13,height=4)
par(mfrow=c(1,2),mar=c(2.2,1.5,.5,.6),mgp =c(1.2,.5,0),xpd=FALSE,family="LM Roman 10",cex=1.5)
kf_plot2_start = 1; kf_plot2_end = 200; kf_plot2_freq = 1
alpha_to_plot2 = alpha_to_plot[4]
plot(seq(kf_plot2_start,kf_plot2_end,kf_plot2_freq),kf_true_param_target_loc$m_i_i[seq(kf_plot2_start+1,kf_plot2_end+1,kf_plot2_freq),alpha_to_plot2],type="l",xlab="t",ylab="")#ylab=bquote(alpha[.(i)](t)))
lines(seq(kf_plot2_start,kf_plot2_end,kf_plot2_freq),RML_ROSP_plot2$kf_means[seq(kf_plot2_start,kf_plot2_end,kf_plot2_freq),alpha_to_plot2],col="red")

kf_plot3_start = 99800
kf_plot3_end = 100000
kf_plot3_freq = kf_plot2_freq
alpha_to_plot3 = alpha_to_plot2
plot(seq(kf_plot3_start,kf_plot3_end,kf_plot3_freq),kf_true_param_target_loc$m_i_i[seq(kf_plot3_start+1,kf_plot3_end+1,kf_plot3_freq),alpha_to_plot3],type="l",xlab="t",ylab="")#ylab=bquote(alpha[.(i)](t)))
lines(seq(kf_plot3_start,kf_plot3_end,kf_plot3_freq),RML_ROSP_plot2$kf_means[seq(kf_plot3_start,kf_plot3_end,kf_plot3_freq),alpha_to_plot3],col="red")
dev.off()

#################
