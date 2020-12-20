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
    n = 6; t = 40000; dt = 1; K = 21; nu = 1
  }
  
  ## RML & OSP Times ##
  if(TRUE){
    t_RML = 49; t_ROSP = 1
  }
  
  ## Initial Observation Locations ##
  if(TRUE){
    
    ## Coordinates ##
    obs0 = list(c(3.4,3.4),c(3.4,4.1),c(4.1,3.4),c(4.1,4.1),
                c(1.5,1.5),c(1.5,3),c(1.5,4.5),c(1.5,6),
                c(3,1.5),c(3,3),c(3,4.5),c(3,6),
                c(4.5,1.5),c(4.5,3),c(4.5,4.5),c(4.5,6),
                c(6,1.5),c(6,3),c(6,4.5),c(6,6))
    
    obs0 = lapply(obs0,function(x) x-c(0.75,0.75))
    
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
    
      rho0_t_dep = param_t_dep_func(list(1,16000),list(param_sig[1],0.5*param_sig[1]))
      sigma2_t_dep = param_t_dep_func(list(1,10000),list(param_sig[2],2*param_sig[2]))
      zeta_t_dep = param_t_dep_func(list(1,27000),list(param_sig[3],0.2))
      rho1_t_dep = param_t_dep_func(list(1,22000),list(param_sig[4],2*param_sig[4]))
      gamma_t_dep = param_t_dep_func(list(1,29000),list(param_sig[5],1.2))
      alpha_t_dep = param_t_dep_func(list(1,25000),list(param_sig[6],1.1))
      mu_x_t_dep = param_t_dep_func(list(1,13000),list(param_sig[7],.3*param_sig[7]))
      mu_y_t_dep = param_t_dep_func(list(1,19000),list(param_sig[8],param_sig[8]+0.2))
      
      param_sig_t_dep = param_sig_t_dep_func(rho0_t_dep,sigma2_t_dep,zeta_t_dep,
                                             rho1_t_dep,gamma_t_dep,alpha_t_dep,
                                             mu_x_t_dep,mu_y_t_dep)
      
    ## Observation parameters ##
    m = 1; m_indices = list(obs0_indices); tau2_vals = list(0.01)
    param_obs = param_obs_func(m,m_indices,tau2_vals)
    param_obs_t_dep =  param_t_dep_func(list(1,4000,7000,15000,19000),
                                        list(param_obs,
                                             param_obs_func(1,m_indices,list(0.05)),
                                             param_obs_func(1,m_indices,list(0.03)),
                                             param_obs_func(1,m_indices,list(0.1)),
                                             param_obs_func(1,m_indices,list(0.02))))
    
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
    param0 = c(rho0 = 0.25, sigma2 = 0.5, zeta = 0.3,rho1 = 0.2, gamma = 1.5, 
               alpha = pi/3, mu_x = 0.1, mu_y = -0.15, tau2 = 0.1)
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
  if(TRUE){
    spde_sim = spde_simulation(n,t,param_sig_t_dep,param_obs_t_dep,
                               param_bias_t_dep,seed = 2,n_obs = n_obs,
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
    grad_obs = 1:4
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
    step_size_vec = c(0.001,0.002,0.0015,0.0015,0.006,0.0015,0.0002,0.0002,0.0001)
    step_size_vec = c(0.001,0.002,0.0015,0.0015,0.008,0.0015,0.0004,0.0004,0.0001)
    
    step_size_RML = matrix(0,nrow=n_iterations*t,ncol=n_param)
    for (i in 1:dim(step_size_RML)[1]){
        step_size_RML[i,] = step_size_vec
    }
      
  }
  
  ## ROSP Step sizes ##
  if(TRUE){
    
    ## Initial step size vector ##
    step_size_vec = c(.01,.01)
    
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
    W_coords = list(c(1,1),c(2,5),c(5,3),c(4,1))
    
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
    diag(W)[unlist(W_indices)] = rep(1000,length(W_coords))
    
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
    plot_limits[[9]] = c(0,0.15)
    plot_limits[[10]] = c(0,2.5)
  }
  
  ## Compute mean parameter estimate? ##
  if(TRUE){
    mean_param = FALSE
  }
  
  ## Parameter estimation period ##
  if(TRUE){
    
    param_est_starts = list(c(4000,19000),c(4000,12000),c(5000,31000),c(4000,24000),
                            c(5000,33000),c(5000,27000),c(4000,15000),c(4000,21000),
                            c(500,4500,7500,15500,19500))
    param_est_ends = list(c(16000,t*n_iterations),c(10000,t*n_iterations),
                          c(27000,t*n_iterations),c(22000,t*n_iterations),
                          c(31000,t*n_iterations),c(25000,t*n_iterations),
                          c(13000,t*n_iterations),c(19000,t*n_iterations),
                          c(4000,7000,15000,19000,t*n_iterations))
  }
  
  ## RML Plot Point Frequency ##
  if(TRUE){
    RML_plot_point_freq = 25
  }
  
  ## RML Plot Frequency ##
  if(TRUE){
    RML_plot_freq = 2000
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
    ROSP_plot_freq = 20000
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
  RML_ROSP_plot3 = RML_ROSP(spde_sim = spde_sim, n = n, n_obs = n_obs, t = t, 
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
saveRDS(RML_ROSP_plot3,"RML_ROSP_plot3")

##################

###############
#### PLOTS ####
###############

## Plot (RML) ##
if(TRUE){
setwd(fig_wd)
#pdf("RML_plot3.pdf",width=7,height=7)
tikz('RML_plot3.tex', standAlone=TRUE, height = 7, width = 7)
par(mar=c(3,4.3,.6,2),mgp=c(2.3,1,0),family="LM Roman 10",cex=2)
par(mfrow=c(3,3))
plot_points = seq(1,t,25)
if(FALSE){
  setwd("/Users/ls616/Google Drive/MPE CDT/PhD/Year 2/Application_Paper/Current/Draft5/Figures/Section4/RML_ROSP_plot3")
  RML_ROSP_plot3 = readRDS("RML_ROSP_plot3")
}
param_est = RML_ROSP_plot3$param_est
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
tools::texi2dvi('RML_plot3.tex',pdf=T)
}

## Plot (ROSP) ##
if(TRUE){
setwd(fig_wd)
#pdf("ROSP_plot3_2d.pdf",width=7,height=7)
tikz('ROSP_plot3_2d.tex', standAlone=TRUE, height = 7, width = 7)
par(mar=c(3.2,3.2,.5,.5),mgp=c(1.7,.7,0),family="LM Roman 10")
plot_points = seq(1,t,ROSP_plot_point_freq)
obs = RML_ROSP_plot3$obs

## Initial Point ##
plot(obs[[1]][1,1],obs[[1]][1,2],xlim = c(0,n),ylim = c(0,n),
     cex = 2, pch = 21, col=brewer.pal(8,"Greens")[7],
     bg = brewer.pal(8,"Greens")[2], xlab="$x$",
     ylab="$y$",cex.lab = 1.5,cex.axis = 1.5,
     xaxt ="n", yaxt = "n", xaxs = "i", yaxs = "i")

axis(side=1,at=seq(0,n,length.out=6),labels = c(0,0.2,0.4,0.6,0.8,1.0),cex.lab = 1.5,cex.axis = 1.5)
axis(side=2,at=seq(0,n,length.out=6),labels = c(0,0.2,0.4,0.6,0.8,1.0),cex.lab = 1.5, cex.axis = 1.5)
      
for (k in grad_obs){
  points(obs[[k]][1,1],obs[[k]][1,2], cex = 2, pch = 21, 
        col=brewer.pal(8,"Greens")[7],
        bg = brewer.pal(8,"Greens")[2])
}
    
## Legend ##
legend(x=3.25, y = 5.15, pt.cex = c(2,2,2,2),
       pch = c(21,21,21),
       legend=c("Initial Sensor Locations",
                "Fixed Sensor Locations",
                "Target Sensor Locations", 
                "Final Sensor Locations"),
       pt.bg = c(brewer.pal(8,"Greens")[2],
                 brewer.pal(8,"Reds")[2],
                 brewer.pal(8,"Purples")[2],
                 brewer.pal(8,"Blues")[2]),
       col = c(brewer.pal(8,"Greens")[7],
               brewer.pal(8,"Reds")[7],
               brewer.pal(8,"Purples")[7],
               brewer.pal(8,"Blues")[7]),
       cex = 1.4)

abline(h=1:(n-1),v=1:(n-1),col="gray90",lty=2)
    
## Weighted coordinates ##
for(k in 1:length(W_coords)){
  points(W_coords[[k]][1],W_coords[[k]][2],cex=2,
         pch=21,col=brewer.pal(8,"Purples")[7],
         bg = brewer.pal(8,"Purples")[2])
}
    
## Plot Final Point ##
for (k in grad_obs){
  points(obs[[k]][n_iterations*t,1],obs[[k]][n_iterations*t,2],
         cex = 2, pch = 21, col = brewer.pal(8,"Blues")[7],
         bg = brewer.pal(8,"Blues")[2])
}
for (k in (1:n_obs)[-grad_obs]){
  points(obs[[k]][n_iterations*t,1],obs[[k]][n_iterations*t,2],
         cex = 2, pch = 21, col = brewer.pal(8,"Reds")[7],
         bg = brewer.pal(8,"Reds")[2])
}
    
## Plot remaining points ##
for (k in grad_obs){
  for (l in plot_points){
    points(obs[[k]][l,1],obs[[k]][l,2],
           col=colorRampPalette(c("grey10", "grey100"))(100)[5*k-4],cex=.6)
  }
}

dev.off()
tools::texi2dvi('ROSP_plot3_2d.tex',pdf=T)
}

###############