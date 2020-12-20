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
    n = 12; t = 5000; dt = 1; K = 21; nu = 1
  }
  
  ## RML & OSP Times ##
  if(TRUE){
    t_RML = 48; t_ROSP = 2
  }
  
  ## Initial Observation Locations ##
  if(TRUE){
    
    ## Coordinates ##
    obs0 = list(c(10.1,7.8),
      c(3.1,5.01),c(8.2,9.25),
      c(10.2,4.02),c(3.2,9.1),
      c(6.1,2.1))#,c(1.01,2.8))
    #,c(3,1))
    
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
           main="Sensor Locations",xaxs="i",yaxs="i")
      abline(v=0:(n),h=(0:n),col="gray92",lty=2)
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
    m = 3; 
    m_indices = list(obs0_indices[1:2],obs0_indices[3:4],obs0_indices[5:6])
    tau2_vals = list(0.01,0.03,0.10)
    param_obs = param_obs_func(m,m_indices,tau2_vals)
    param_obs_t_dep =  param_t_dep_func(list(1),list(param_obs))
    
    ## Observation bias ##
    m1 = 2
    m1_indices = list(obs0_indices[1],obs0_indices[2:3])
    bias_vals = list(0,2)
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
    param0 = c(rho0 = 0.25, sigma2 = 0.5, zeta = 0.3,rho1 = 0.2, gamma = 1.5, 
               alpha = pi/3, mu_x = 0.1, mu_y = -0.15, tau2 = 0.1)
    param_sig0 = param_sig[1:8]
    param_obs0 = param_obs_func(m,m_indices,list(0.10,0.30,0.20))
    param_bias0 = param_bias_func(m1,m1_indices,list(0.9,1.1))
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
    grad_params = 9:n_param
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
    
    # Step Sizes #
    step_size_vec = 20*c(0.002,0.005,0.01,0.001,0.015,0.005,0.001,0.001,0.0004,0.0006,0.0004,0.0004,0.0004)
    
    step_size_RML = matrix(0,nrow=n_iterations*t,ncol=n_param)
    for (i in 1:dim(step_size_RML)[1]){
      step_size_RML[i,] = step_size_vec*i^(-0.505)
    }
    
  }
  
  ## ROSP Step sizes ##
  if(TRUE){
    
    ## Initial step size vector ##
    step_size_vec = c(.5,.5)
    
    ## Step size matrix ##
    step_size_mat = matrix(0,nrow=n_iterations*t,ncol=2)
    for (i in 1:(n_iterations*t)){
      step_size_mat[i,]=step_size_vec*i^(-0.7)
    }
    
    ## Step size matrices for all sensors ##
    step_size_ROSP = rep(list(step_size_mat),n_obs)
  }
  
  ## Weighting matrix ##
  if(TRUE){
    
    ## Weighting coordinates ##
    W_coords = list(c(3,2),
                    c(5,10),
                    c(11,6))
                    #c(10,11),
                    #c(7,5),
                    #c(1,7))
                    #c(6,8),
                    #c(4,4),
                    #c(9,6),#,c(1,7),c(7,10),
                    #c(1,1))
                    #c(3,10))
    
    plot_obs = TRUE
    if(plot_obs){
      par(mfrow=c(1,1))
      plot(sapply(obs0,function(x) x)[1,],sapply(obs0,function(x) x)[2,],
           xlim = c(0,n),ylim=c(0,n),xlab=expression(x),ylab=expression(y),
           main="Sensor Locations",xaxs="i",yaxs="i")
      points(sapply(W_coords,function(x) x)[1,], sapply(W_coords,function(x) x)[2,],pch=19)
      abline(v=0:(n),h=(0:n),col="gray92",lty=2)
    }
    
    ## Weighting coordinates at all times ##
    target_loc = lapply(1:length(W_coords),
                        function(j) matrix(W_coords[[j]],nrow=t,ncol=2,byrow=T))
    
    ## Weighting indices ##
    W_indices = lapply(W_coords, function(A) A[2]*n+A[1]+1)
    
    ## Weighting matrix (physical space) ##
    W = matrix(0,nrow=n^2,ncol=n^2)
    diag(W)[unlist(W_indices)] = rep(10,length(W_coords))
    
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
    plot_limits = rep(list(rep(0,2)),9)
    plot_limits[[1]] = c(0.2,0.6)
    plot_limits[[2]] = c(0.1,0.7)
    plot_limits[[3]] = c(0.1,0.6)
    plot_limits[[4]] = c(0.05,0.4)
    plot_limits[[5]] = c(0.9,2.3)
    plot_limits[[6]] = c(0.6,1.2)
    plot_limits[[7]] = c(0.05,0.35)
    plot_limits[[8]] = c(-0.35,0)
    plot_limits[[9]] = c(0,0.40)
    plot_limits[[10]] = c(0,2.5)
  }
  
  ## Compute mean parameter estimate? ##
  if(TRUE){
    mean_param = FALSE
  }
  
  ## Parameter estimation period ##
  if(TRUE){
    param_est_starts = list(5000,5000,10000,5000,10000,4000,4000,4000,3000,3000,3000,3000,3000)
    param_est_ends = rep(list(n_iterations*t),n_param)
  }
  
  ## RML Plot Point Frequency ##
  if(TRUE){
    RML_plot_point_freq = 25
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
    ROSP_plot = FALSE
    ROSP_plot2d = TRUE
  }
  
  ## ROSP Plot Point Frequency ##
  if(TRUE){
    ROSP_plot_point_freq = 50
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
  RML_ROSP_plot4 = RML_ROSP(spde_sim = spde_sim, n = n, n_obs = n_obs, t = t, 
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
saveRDS(RML_ROSP_plot4,"RML_ROSP_plot4")

##################

###############
#### PLOTS ####
###############

## Plots (RML) ##
if(TRUE){
  setwd(fig_wd)
  #pdf("RML_plot4.pdf",width=7,height=7)
  if(FALSE){
    setwd("/Users/ls616/Google Drive/MPE CDT/PhD/Year 2/Application_Paper/Current/Draft5/Figures/Section4/RML_ROSP_plot4")
    RML_ROSP_plot2 = readRDS("RML_ROSP_plot4")
  }
  param_est = RML_ROSP_plot4$param_est
  param_sig_indices = 1:8
  param_sig_true = param_sig_t_dep
  n_param_sig = 8
  param_obs_indices = 9:11
  param_obs_true = param_obs_t_dep
  n_param_obs = 3
  param_bias_indices = 12:13
  param_bias_true = param_bias_t_dep
  param_sig_names = list("$\\hat{\\rho}_0$","$\\hat{\\sigma}^2$",
                         "$\\hat{\\zeta}$","$\\hat{\\rho}_1$",
                         "$\\hat{\\gamma}$","$\\hat{\\alpha}$",
                         "$\\hat{\\mu}_x$","$\\hat{\\mu}_y$")
  param_obs_names = list("$\\hat{\\tau}^2$")
  param_bias_names = list("$\\hat{\\beta}$")
  
  
  ## Observation Parameters ##
  
  library(tikzDevice)
  tikz('RML_plot4_1.tex', standAlone=TRUE, height = 7, width = 9)
  par(mar=c(6,7,.6,2),mgp=c(4.0,1.7,0),family="LM Roman 10",cex=2)
  par(mfrow=c(1,1))
  plot_points = seq(1,t,10)

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
      
      plot(plot_points, param_est[plot_points,k], cex=0.8,
           xlim = c(1,n_iterations*t+1), ylab=y_lab, xlab="$t$",
           ylim = y_lim, main = "",cex.lab = 3.0, cex.axis = 3.0,pch=19)
    }
    
    if(k>9){
      points(plot_points,param_est[plot_points,k],pch=19,cex=0.8)
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
  tools::texi2dvi('RML_plot4_1.tex',pdf=T)
  
  
  ## Bias Parameters ##
  
  tikz('RML_plot4_2.tex', standAlone=TRUE, height = 7, width = 9)
  par(mar=c(6,7,.6,2),mgp=c(4.0,1.7,0),family="LM Roman 10",cex=2)
  par(mfrow=c(1,1))
  plot_points = seq(1,t,10)

  for (k in param_bias_indices){
    if(k==param_bias_indices[[1]]){
      
      y_lim = NULL
      y_lab = expression(theta)
      main = NULL
      
      if(!is.null(plot_limits)){
        y_lim = plot_limits[[10]]
      }
      
      if(!is.null(param_bias_names)){
        y_lab = param_bias_names[[k-(n_param_sig+n_param_obs)]]
        main = param_bias_names[[k-(n_param_sig+n_param_obs)]]
      }
      
      plot(plot_points, param_est[plot_points,k], cex=0.8,
           xlim = c(1,n_iterations*t+1) ,xlab="$t$", ylab=y_lab,
           ylim = y_lim, main = "",cex.lab = 3.0, cex.axis = 3.0,
           pch = 19)
    }
    
    if(k>param_bias_indices[[1]]){
      points(plot_points,param_est[plot_points,k],pch=19,cex=0.8)
    }
    
    if (!is.null(param_bias_true)){
      if(length(param_bias_true$t_vals)==1){
        abline(h=param_bias_true$param_vals[[1]]$bias_vals[[k-(n_param_sig+n_param_obs)]],
               col="red",lty=2)
      }
      else{
        for(l in 1:(length(param_bias_true$t_vals)-1)){
          segments(param_bias_true$t_vals[[l]],
                   param_bias_true$param_vals[[l]]$bias_vals[[k-(n_param_sig+n_param_obs)]],
                   param_bias_true$t_vals[[l+1]],
                   param_bias_true$param_vals[[l]]$bias_vals[[k-(n_param_sig+n_param_obs)]],
                   col="red",lty=2)
          segments(param_bias_true$t_vals[[l+1]],
                   param_bias_true$param_vals[[l]]$bias_vals[[k-(n_param_sig+n_param_obs)]],
                   param_bias_true$t_vals[[l+1]],
                   param_bias_true$param_vals[[l+1]]$bias_vals[[k-(n_param_sig+n_param_obs)]],
                   col="red",lty=2)
        }
        segments(param_bias_true$t_vals[[length(param_bias_true$t_vals)]],
                 param_bias_true$param_vals[[length(param_bias_true$t_vals)]]$bias_vals[[k-(n_param_sig+n_param_obs)]],
                 n_iterations*t+1,
                 param_bias_true$param_vals[[length(param_bias_true$t_vals)]]$bias_vals[[k-(n_param_sig+n_param_obs)]],
                 col="red",lty=2)
      }
    }
  }
  
  dev.off()
  tools::texi2dvi('RML_plot4_2.tex',pdf=T)
  
}

###############

