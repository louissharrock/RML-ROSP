########################
########################
#### SPDE FUNCTIONS ####
########################
########################

################
# COMMAND LINE #
################

## R CMD SHLIB FFBS_spectral.c -lfftw3
## R CMD SHLIB RML_spectral.c
## R CMD SHLIB propagate_spectral.c
## R CMD SHLIB tkf_spectral.c
## R CMD SHLIB tkf2_spectral.c

################

################
# .C FUNCTIONS #
################

## Define directories ##
main_wd = "/Users/ls616/Google Drive/MPE CDT/PhD/Year 1/Code/Current"
fig_wd = "/Users/ls616/Desktop"

## Set directory ##
setwd(main_wd)

## Load functions from C
dyn.load("FFBS_spectral.so")
dyn.load("RML_spectral.so")
dyn.load("propagate_spectral.so")
dyn.load("tkf_spectral.so")
dyn.load("tkf2_spectral.so")

################

#################
# .R LIBRARIES #
#################

## Additional libraries ##
library(mvtnorm)
library(coda)
library(microbenchmark)
library(lattice)
library(dlm)
library(RColorBrewer)
library(gtools)

#################

#################
# MISCALLENEOUS #
#################

## Remove non functions ##
rm_non_func <- function(){
  env <- parent.frame()
  rm(list = setdiff(ls(all.names=TRUE, env = env), 
                    lsf.str(all.names=TRUE, env = env)),
     envir = env)
}

## Colour Scheme ##
colors <- function(){
  return(c("#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF", "#0000DF", 
           "#0000EF", "#0000FF", "#0010FF", "#0020FF", "#0030FF", "#0040FF",
           "#0050FF", "#0060FF", "#0070FF", "#0080FF", "#008FFF", "#009FFF", 
           "#00AFFF", "#00BFFF", "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF",
           "#10FFEF", "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F", 
           "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50", "#BFFF40", 
           "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00", "#FFEF00", "#FFDF00", 
           "#FFCF00", "#FFBF00", "#FFAF00", "#FF9F00", "#FF8F00", "#FF8000",
           "#FF7000", "#FF6000", "#FF5000", "#FF4000", "#FF3000", "#FF2000", 
           "#FF1000", "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000",
           "#AF0000", "#9F0000", "#8F0000", "#800000"))
}

## Compute vector L2 norm ##
vec_norm <- function(x) return(sqrt(sum(x^2)))

## Compute linear predictor ##
linear_predictor <- function(x,beta){
  return(x%*%beta)
}

## Convert space-time matrix to vector (row = time, column = space) ##
STmatrix_to_vec <- function(mat){
  return(as.vector(t(mat)))
}

## Convert vector to space-time matrix (row = time, column = space) ##
vec_to_STmatrix <- function(vec,t=1){
  return(matrix(vec,nrow=t,byrow=TRUE))
}

## Map non-gridded observations to a unit n x n grid ##
map_obs_to_grid <- function(n,y_non_grid,coord,
                            length_x=NULL,length_y=NULL){
  
  ## Inputs ##
  # n: number of points per axis of the grid onto which the obs. are mapped
  # y_non_grid: obs. data in a T x N matrix; coords. of each point in 'coord'
  # coord: N x 2 matrix of the coords. of the N obs. data points
  # length_x: length of x-axis; must be at least as large as max(coord[,1])
  # length_y: length of y-axis; must be at least as large as max(coord[,2])
  
  
  ## Generate points on unit grid ##
  cell_x <- rep(1:n, n)/n - 1/(2*n)
  cell_y <- as.vector(apply(matrix(1:n), 1, rep, times = n))/n - 1/(2*n)
  
  ## Rescale coordinates ##
  
  ## If x-axis scale not provided, rescale coordinates using max value
  if (is.null(length_x)){
    coord[, 1] <- coord[, 1]/max(coord[, 1]) 
  }
  
  ## If x-axis scale provided, rescale accordingly
  else {
    coord[, 1] <- coord[, 1]/length_x
  }
  
  ## Similarly for y axis
  if (is.null(length_y)){
    coord[, 2] <- coord[, 2]/max(coord[, 2])
  }
  else{
    coord[, 2] <- coord[, 2]/length_y
  }
  
  ## Define 'obs_indices' variable
  obs_indices <- rep(0, dim(coord)[1])
  
  ## Populate w/ indices of the grid coordinates which are closest to each of 
  ## the observations
  for (i in 1:dim(coord)[1]) {
    
    ## Compute Euclidean distance 
    d = abs(cell_x-coord[i, 1])^2 + abs(cell_y-coord[i, 2])^2
    
    ## Compute 'obs_indices' as above
    obs_indices[i] <- which(d == min(d))
  }
  
  ## Compute gridded observations values
  y = array(NA,c(dim(y_non_grid)[1],n^2))
  
  for(t in 1:dim(y_non_grid)[1]){
    for(i in unique(obs_indices)){
      y[t,i]=mean(y_non_grid[t,which(i==obs_indices)])
    }
  }
  
  ## Output y
  return(y)
}

#################

## --- 

############################
# WAVE NUMBERS FOR REAL FT #
############################
## -> Generate all wavenumbers k_j used in real Fourier transform, for a given
##    spatial dimension n
## -> Inputs: n is the number of grid points on each (x,y) axis
## -> Output: list containing
##            -- 'wave', a 2 by n^2 matrix with wavenumbers
##            -- 'cosine_indices', a vector of indices indicating
##                column of wavenumbers for cosine terms (first
##                four terms not included)

wave_numbers <- function(n){
  
  ## Output an error message if n not even ##
  if(n%%2!=0){
    print("Error: n must be even")
    return()
  }
  
  ## Otherwise generate wavenumbers ##
  else{
    # Wavenumbers used in the real FT (unscaled by 2pi & unordered)
    x_index <- c(rep(0:(n/2),n/2+1),rep(1:(n/2-1),n/2-1))
    y_index <- c(as.vector(sapply(0:(n/2),function(i) rep(i,n/2+1))),
                 as.vector(sapply((-n/2+1):(-1),function(i) rep(i,n/2-1))))
    xy_index <- rbind(x_index,y_index)
    
    # Indices of wavenumbers which give a zero sine term
    zero_sine_indices <- c(which(x_index==0 & y_index==0),
                   which(x_index==(n/2) & y_index==0),
                   which(x_index==0 & y_index==(n/2)),
                   which(x_index==(n/2) & y_index==(n/2)))
    
    # Indices of cos basis functions (where cos & sin terms both non-zero) are
    # 5,7,9,11,...,n^2-1
    cosine_indices <- 2*(1:(dim(xy_index)[2]-4))-1+4
    
    ## Generate empty array to contain ordered & scaled wavenumbers ##
    wave_numbers <- array(0,c(2,n^2))
    
    # First 4 terms are wavenumbers for cosine-only terms
    wave_numbers[,1:4] <- xy_index[,zero_sine_indices]
    
    # Remaining n^2-4 terms are wavenumbers of (cosine, sine) terms
    wave_numbers[,cosine_indices] <- xy_index[,-zero_sine_indices]
    wave_numbers[,cosine_indices+1] <- xy_index[,-zero_sine_indices]
    
    # Scale wavenumbers by 2*pi
    wave_numbers <- wave_numbers*2*pi
    
    # Output
    return(list(wave=wave_numbers,cosine_indices=cosine_indices))
  }
}
############################

###############################
# REAL AND COMPLEX FT INDICES #
###############################
## -> Auxilary function for conversion between complex FFT and the real FT
## -> Inputs: n is the number of grid points on each (x,y) axis
## -> Output: list containing
##            -- 'cosine_indices', a vector of the indices of the wave-numbers 
##                of cosine terms for the real FT (not inc. first 4 cosine only
##                terms)
##            -- 'exp_indices', a vector of the indices of the wave-numbers in 
##                the complex FT, whic correspond to those used in the real FT 
##            -- 'exp_indices_conj', a vector of the indices of the wave-numbers
##                in the complex FT not used in the real FT

FFT_index_complex_and_real <- function(n){
  
  # Initialise matrices for complex and real wave numbers
  wave_complex <- array(0,c(2,n*n))
  wave_real <- array(0,c(2,n*n))
  
  # Wavenumbers used in the real FT (unscaled by 2pi & unordered)
  x_index <- c(rep(0:(n/2),n/2+1),rep(1:(n/2-1),n/2-1))
  y_index <- c(as.vector(sapply(0:(n/2),function(i) rep(i,n/2+1))),
               as.vector(sapply((-n/2+1):(-1),function(i) rep(i,n/2-1))))
  xy_index <- rbind(x_index,y_index)
  
  # Indices of wavenumbers which give a zero sine term
  zero_sine_indices <- c(which(x_index==0 & y_index==0),
                         which(x_index==(n/2) & y_index==0),
                         which(x_index==0 & y_index==(n/2)),
                         which(x_index==(n/2) & y_index==(n/2)))
  
  # Indices of cos basis functions (where cos & sin terms both non-zero) are
  # 5,7,9,11,...,n^2-1 
  cosine_indices <- seq(from=5,to=n^2,by=2)
  
  ## Real FT ##
  
  # First 4 terms are wavenumbers for cosine-only terms
  wave_real[,1:4] <- xy_index[,zero_sine_indices]
  
  # Remaining n^2-4 terms are wavenumbers of (cosine, sine) terms
  wave_real[,cosine_indices] <- xy_index[,-zero_sine_indices]
  wave_real[,cosine_indices+1] <- xy_index[,-zero_sine_indices]
  
  # Instead of wave-numbers in 0<=x<=n/2 and -n/2+1<=y<=n/2, use (equivalently)
  # wave-numbers in 0<=x<=n/2 and 0<=y<=n-1 
  wave_real[2,wave_real[2,]<0] <- wave_real[2,wave_real[2,]<0]+n
  
  
  ## Complex FT ##
  
  # Wavenumbers used in the complex FT (unscaled by 2pi & unordered), i.e. all 
  # grid points in 0<=x<=n^2-1, 0<=y<=n^2-1
  wave_complex <- rbind(rep(0:(n-1),n),
                        as.vector(apply(matrix(0:(n-1)),1,rep,times=n))) 
  
  # Indices of the wave-numbers in complex FT which correspond to wave-numbers
  # in real FT
  exp_indices <- match(wave_real[1,]+1i*wave_real[2,],
                       wave_complex[1,]+1i*wave_complex[2,])
  
  waveCon <- (n-wave_real[,cosine_indices])%%n
  
  # Indices of the wave-numbers in complex FT which don't correspond to wave-
  # numbers in real FT
  exp_indices_conj <- match(waveCon[1,]+1i*waveCon[2,],
                            wave_complex[1,]+1i*wave_complex[2,])
  
  ## Ouputs ##
  return(list(cosine_indices=cosine_indices,exp_indices=exp_indices,
              exp_indices_conj=exp_indices_conj))
}
###############################

## --- 

################################
# MATERN COV. SPECTRAL DENSITY #
################################
## -> Compute spectrum of the Matern covariance function
## -> Inputs:
##            -- 'wave', spatial wavenumbers (from wave_numbers$wave)
##            -- 'n', number of grid points on each ax_is (x,y)
##            -- 'n_cos', number of cosine only terms (default = 4)
##            -- 'rho0', range 
##            -- 'sigma2', marginal variance
##            -- 'nu', smoothness of Matern covariance
##            -- 'norm', if TRUE, spectrum is normalised
##            -- 'd', the Eucliden distance between two points (default=2)
## -> Output: 
##            -- 'spectrum', vector of Matern spectrum

matern_spectrum <- function(wave,n,n_cos=4,rho0,sigma2,nu=1,d=2,norm=TRUE){
  
  ## Error message if nu<=0 ##
  if(nu<=0){
    print("Error: nu must be positive")
    return()
  }
  
  ## Error message if sigma^2<=0 ##
  if(sigma2<=0){
    print("Error: sigma^2 must be positive")
    return()
  }
  
  ## Compute k_j^T k_j ##
  w2 <- apply(wave^2,2,sum)
  
  ## Compute spectrum (up to proportionality) ##
  numerator <- 1
  denominator <- ((1/rho0)^2 + w2)^(nu + 1)
  spectrum <- numerator/denominator

  ## Adjustment for cosine only terms ##
  spectrum[1:n_cos] <- spectrum[1:n_cos]/2
  
  ## Normalise ##
  if(norm){
    spectrum <- spectrum*(n^2)/sum(spectrum)
  }
  else{
    spectrum <- spectrum/sum(spectrum)
  }
  
  ## Rescale by sigma^2 ##
  spectrum <- sigma2 * spectrum
  return(spectrum)
}
################################

###############################
# INNOVATION SPECTRAL DENSITY #
###############################
## -> Compute spectrum of the innovation term (epsilon)
## -> Inputs:
##            -- 'wave', spatial wavenumbers (from wave_numbers$wave)
##            -- 'n', number of grid points on each ax_is (x,y)
##            -- 'n_cos', number of cosine only terms (default = 4)
##            -- 'rho0', range of Matern covariance function for innovation
##            -- 'sigma2', marginal variance of Mat. cov. fn. for innovation
##            -- 'zeta', damping parameter
##            -- 'rho1', range parameter of diffusion term
##            -- 'alpha', direction of anisotropy of diffusion term
##            -- 'gamma', amount of anisotropy of diffusion term
##            -- 'nu', smoothness of Matern covariance for innovations (default = 1)
##            -- 'dt', temporal lag between time points (default = 1)
##            -- 'norm', if TRUE, spectrum is normalised
## -> Output: 
##            -- 'Qhat', vector of spectrum of integrated innovation term
innovation_spectrum <- function(wave,n,n_cos=4,rho0,sigma2,
                                zeta,rho1,alpha,gamma,nu=1,
                                dt=1,norm=TRUE){
  
  ## Error messages if paramaters outside correct range ##
  if(nu<0){
    print("Error: nu must be positive")
    return()
  }
  if(sigma2<0){
    print("Error: sigma^2 must be positive")
    return()
  }
  if(zeta<0){
    print("Error: zeta must be positive")
    return()
  }
  if(gamma<0){
    print("Error: gamma must be positive")
    return()
  }
  if(alpha<0 & alpha>=(pi/2)){
    print("Error: alpha needs must be between 0 and pi/2")
    return()
  }
  
  ## Generate Matern spectrum using 'matern_spectrum' ##
  spec <- matern_spectrum(wave=wave,n=n,n_cos=n_cos,rho0=rho0,sigma2=1,
                          nu=nu,d=2,norm=norm)
  
  ## Number of Fourier terms ##
  K <- dim(wave)[2]
  
  ## Compute Sigma ## 
  
  # Special Case if rho1 = 0
  if(rho1==0){
    Sig <- cbind(c(0,0),c(0,0))
  }
  
  # General Case
  else{
    Sig_Chol <- cbind(c(cos(alpha),-gamma*sin(alpha)),
                      c(sin(alpha),gamma*cos(alpha)))/rho1
    Sig <- solve(t(Sig_Chol)%*%Sig_Chol)
  }
  
  ## Compute -(k_j^T * Sigma * k_j + zeta) ##
  ExpTerm <- -apply(wave*Sig%*%wave,2,sum)-rep(zeta,K)
  
  ## Compute Q hat ##
  Qhat <- sigma2*spec*(1-exp(2*dt*ExpTerm))/(-2*ExpTerm)
  
  ## Output ##
  return(Qhat)
}
###############################

## ---

##################
# REAL FT MATRIX #
##################
## -> Compute matrix which applies 2d real FT, for a set of n_obs spatial
##    locations
## -> Inputs:
##            -- 'wave', matrix of spatial wavenumbers (from wave_numbers$wave)
##            -- 'cosine_indices', vector of indices of columns of cos terms
##            -- 'n_cos', number of real Fourier terms with only cosine term
##            -- 'n', spatial dimensions
##            -- 'n_obs', the number of observations
##            -- 'coords', a list containing the coordinates of each of
##               the 'n_obs' observation locations, scaled by n [e.g. the input
##               (1,1) corresponds to actual coordinates (1/n,1/n)]

## -> Output: 
##            -- 'phi', a n_obs x K matrix which applies the 2D real FT

FFT_real_matrix <- function(wave,cosine_indices,n_cos=4,n,n_obs=NULL,
                            coords=NULL){
  
  ## If n_obs not provided, set n_obs = n^2 ##
  if(is.null(n_obs)){
    n_obs = n^2
  }
  
  ## If coords not provided, set coords = uniform n x n grid ##
  if(is.null(coords)){
    coords = lapply(1:(n^2),
                    function(i) as.numeric(expand.grid(0:(n-1),0:(n-1))[i,]))
  }
  
  ## Normalising Constants (Check these?)
  A <- sqrt(n^2/2)
  B <- sqrt(n^2)
  
  ## Number of Fourier basis functions
  K = dim(wave)[2]
  
  ## Initialise matrices for storage
  phi_row <- rep(0, K)
  phi <- matrix(0, nrow = n_obs, ncol = K)
  
  ## Loop over matrix indices
  for(i in 1:n_obs){
    
    # Cosine only terms
    phi_row[1:n_cos] <- cos(coords[[i]]%*%(wave[,1:n_cos])/n)/B
    
    # 'Other' cosine terms
    phi_row[cosine_indices] <- cos(coords[[i]]%*%(wave[,cosine_indices])/n)/A
    
    # 'Other' sine terms
    phi_row[cosine_indices+1] <- sin(coords[[i]]%*%(wave[,cosine_indices+1])/n)/A
    
    # Add to matrix
    phi[i,] <- phi_row
  }
  
  # Output
  return(phi)
}

##################

###############################
# REAL FT MATRIX SPATIAL GRAD #
###############################
## -> Compute gradient of matrix which applies 2d real FT, for a set of n_obs 
##    spatial locations
## -> Inputs:
##            -- 'wave', matrix of spatial wavenumbers (from wave_numbers$wave)
##            -- 'cosine_indices', vector of indices of columns of cos terms
##            -- 'n_cos', number of real Fourier terms with only cosine term
##            -- 'n', spatial dimensions
##            -- 'n_obs', the number of observations
##            -- 'coords', a list containing the coordinates of each of
##               the 'n_obs' observation locations, scaled by n [e.g. the input
##               (1,1) corresponds to actual coordinates (1/n,1/n)]
##            -- 'grad_coord_index', the indices of the coordinates with respect 
##                to which to compute the gradient

## -> Output: 
##            -- 'phi', a list of length 2*length(grad_coord_index), containing
##                the gradients w.r.t x,y of the locations whose index in 
##                coords is given by grad_coord_index

FFT_real_matrix_grad <- function(wave,cosine_indices,n_cos=4,n,n_obs,coords,
                                 grad_coord_index=1:n_obs){
  
  ## Normalising Constants (Check these?)
  A <- sqrt(n^2/2)
  B <- sqrt(n^2)
  
  ## Number of Fourier basis functions
  K = dim(wave)[2]
  
  ## Initialise matrices
  phi_grad_x <- matrix(0,nrow=n_obs,ncol=K)
  phi_grad_y <- matrix(0,nrow=n_obs,ncol=K)
  
  ## Initalise rows
  phi_row_grad_x <- rep(0,K)
  phi_row_grad_y <- rep(0,K)
  
  ## Loop over matrix indices
  for(i in 1:n_obs){
    
    # Grad of cosine only terms
    phi_row_grad_x[1:n_cos] <- -wave[1,1:n_cos]*sin(coords[[i]]%*%(wave[,1:n_cos])/n)/B
    phi_row_grad_y[1:n_cos] <- -wave[2,1:n_cos]*sin(coords[[i]]%*%(wave[,1:n_cos])/n)/B
    
    # Grad of 'other' cosine terms
    phi_row_grad_x[cosine_indices] <- -wave[1,cosine_indices]*sin(coords[[i]]%*%(wave[,cosine_indices])/n)/A
    phi_row_grad_y[cosine_indices] <- -wave[2,cosine_indices]*sin(coords[[i]]%*%(wave[,cosine_indices])/n)/A
    
    # Grad of 'other' sine terms
    phi_row_grad_x[cosine_indices+1] <- wave[1,cosine_indices+1]*cos(coords[[i]]%*%(wave[,cosine_indices+1])/n)/A
    phi_row_grad_y[cosine_indices+1] <- wave[2,cosine_indices+1]*cos(coords[[i]]%*%(wave[,cosine_indices+1])/n)/A
    
    # Add to grad. matrix
    phi_grad_x[i,] <- phi_row_grad_x
    phi_grad_y[i,] <- phi_row_grad_y
  }
  
  ## Initialise output ##
  output = list()
  
  ## Set non grad-coords to zero
  if(!is.null(grad_coord_index)){
    for (i in grad_coord_index){
      phi_grad_x_tmp = phi_grad_x
      phi_grad_y_tmp = phi_grad_y
      phi_grad_x_tmp[-i,] = rep(0,K)
      phi_grad_y_tmp[-i,] = rep(0,K)
      output = c(output,list(phi_grad_x_tmp,phi_grad_y_tmp))
    }
  }
  else{
    output = list(phi_grad_x,phi_grad_y)
  }
  
  # Output
  return(output)
}

###############################

## ---

#######################
# WEIGHTING FUNCTIONS #
#######################

## Bivariate Normal ##
bivar_norm = function(x,y,centre_x,centre_y,sigma_x,sigma_y,rho){
  
  z = (x-centre_x)^2/sigma_x^2 + (y-centre_y)^2/sigma_y^2 -
    2*rho*(x-centre_x)*(y-centre_y)/(sigma_x*sigma_y) 
  
  weights = 1/(2*pi*sigma_x*sigma_y*sqrt(1-rho^2))*exp(-z/(2*(1-rho^2)))
  
  return(weights)
}
bivar_norm_weights = function(n,centre_x,centre_y,sigma_x,sigma_y,rho){
  
  ## Generate grid points ##
  grid = expand.grid(x = 0:(n-1),y = 0:(n-1))/n
  
  ## Compute weights using 'bivar_norm' ##
  weights = with(grid,bivar_norm(x = x,y = y, centre_x = centre_x, 
                                 centre_y = centre_y, sigma_x = sigma_x, 
                                 sigma_y=sigma_y, rho = rho))
  
  ## Ouput ##
  return(weights)
  
}

## Sech ##
sech2d = function(x,y,centre_x,centre_y,sigma_x,sigma_y,rho){
  
  z = (x-centre_x)^2/sigma_x^2 + (y-centre_y)^2/sigma_y^2 -
    2*rho*(x-centre_x)*(y-centre_y)/(sigma_x*sigma_y) 
  
  weights = pracma::sech(sqrt(z))
  
  return(weights)
}
sech2d_weights = function(n,centre_x,centre_y,sigma_x,sigma_y,rho){
  
  ## Generate grid points ##
  grid = expand.grid(x = 0:(n-1),y = 0:(n-1))/n
  
  ## Compute weights using 'bivar_norm' ##
  weights = with(grid,sech2d(x = x,y = y, centre_x = centre_x, 
                             centre_y = centre_y, sigma_x = sigma_x, 
                             sigma_y=sigma_y))
  
  ## Ouput ##
  return(weights)
  
}

#######################

####################
# WEIGHTING MATRIX #
####################

weights_mat_func = function(weights){
  weights_mat = diag(weights)
  return(weights_mat)
}

####################

#########################
# SPEC WEIGHTING MATRIX #
#########################

spec_weights_mat_func = function(n, K, weights = NULL, 
                                 weights_mat = NULL){
  
  ## SPDE_FT Object ##
  SPDE_FT = spde_initialise(n,t,K)
  wave = SPDE_FT$wave
  cosine_indices = SPDE_FT$cosine_indices
  n_cos = SPDE_FT$n_cos
  n_cos_sin = length(cosine_indices)
  K = SPDE_FT$K
  
  ## Weighting Matrix ##
  if(is.null(weights_mat)){
    weights_mat = weights_mat_func(weights)
  }
  
  ## Phi ##
  Phi = FFT_real_matrix(wave,cosine_indices,n_cos,n)
  
  ## Spectral Weighting Matrix ##
  spec_weights_mat = t(Phi)%*%weights_mat%*%Phi

  ## Output ##
  return(spec_weights_mat)
}

#########################

## ---

##############################
# SPECTRAL PROPAGATOR MATRIX #
##############################
## -> Compute spectral propagator matrix G of vector autoregressive model
## -> Inputs:
##            -- 'wave', matrix of spatial wavenumbers (from wave_numbers$wave)
##            -- 'cosine_indices', vector of indices of columns of cos terms
##            -- 'zeta', damping parameter
##            -- 'rho1', range parameter of diffusion term
##            -- 'gamma', amount of anisotropy in diffusion
##            -- 'alpha', direction of anisotropy in diffusion
##            -- 'mu_x', x component of drift
##            -- 'mu_y', y component of drift
##            -- 'dt', temporal lag between time points (default = 1)
##            -- 'n_cos', number of real Fourier terms with only cosine term

## -> Output: 
##            -- 'G', the propagator matrix

propagator_matrix <- function(wave,cosine_indices,zeta,rho1,gamma,alpha,
                                mu_x,mu_y,dt=1,n_cos=4){
  
  ## Number of terms in Fourier expansion
  K <- dim(wave)[2]
  
  ## Special case when rho = 0
  if(rho1==0){
    Sig <- cbind(c(0,0),c(0,0))
  }
  
  ## Compute Sigma (see pg 4 of r_spate documentation)
  else{
    Sig_Chol <- cbind(c(cos(alpha),-gamma*sin(alpha)),
                      c(sin(alpha),gamma*cos(alpha)))/rho1
    Sig <- solve(t(Sig_Chol)%*%Sig_Chol)
  }
  
  ## Compute -Delta*(k_j^T * Sigma * k_j + zeta)
  H_diag <- -dt*apply(wave*Sig%*%wave,2,sum)-dt*rep(zeta,K)
  
  ## Compute -Delta * mu * k_j
  H_off_diag <- dt*c(mu_x,mu_y)%*%wave
  
  ## Compute G
  
  # Initialise
  G <- matrix(0,ncol=K,nrow=K)
  
  # Diagonal elements
  diag(G) <- c(exp(H_diag[1:n_cos]),
               exp(H_diag[(n_cos+1):K])*cos(H_off_diag[(n_cos+1):K]))
  
  # Off diagonal elements
  diag(G[cosine_indices,cosine_indices+1]) <-
    -exp(H_diag[cosine_indices])*sin(H_off_diag[cosine_indices])
  
  diag(G[cosine_indices+1,cosine_indices]) <- 
    exp(H_diag[cosine_indices])*sin(H_off_diag[cosine_indices])
  
  
  ## Alternative: straightforward but MUCH slower 
  
  # Compute Delta*H matrix
  # H = matrix(0,ncol=K,nrow=K)
  # diag(H) = H_diag
  # for (i in cosine_indices){
  #  H[i,i+1]=-H_off_diag[i]
  #  H[i+1,i]=H_off_diag[i]
  #}
  
  # Compute G = exp(Delta*H)
  # library(expm)
  # G <- expm(H)
  
  ## Output
  return(G)
}
##############################

####################################
# SPECTRAL PROPAGATOR MATRIX (VEC) #
####################################
## -> Compute spectral propagator matrix G of vector autoregressive model, in 
##    vector form
## -> Inputs:
##            -- 'wave', matrix of spatial wavenumbers (from wave_numbers$wave)
##            -- 'cosine_indices', vector of indices of columns of cos terms
##            -- 'zeta', damping parameter
##            -- 'rho1', range parameter of diffusion term
##            -- 'gamma', amount of anisotropy in diffusion
##            -- 'alpha', direction of anisotropy in diffusion
##            -- 'mu_x', x component of drift
##            -- 'mu_y', y component of drift
##            -- 'dt', temporal lag between time points (default = 1)
##            -- 'n_cos', number of real Fourier terms with only cosine term

## -> Output: list containing
##            -- 'G_cos', first n_cos diagonal entries of G, corresponding to 
##                cosine only terms
##            -- 'G_1', remaining diagonal entries of G, corresponding to other 
##                cosine terms 
##            -- 'G_2', off (upper) diagonal entries of G, corresponding to 
##                G_1 entries

propagator_vec <- function(wave,cosine_indices,zeta,rho1,gamma,alpha,mu_x,mu_y,
                           dt=1,n_cos=4){
  
  ## Number of terms in Fourier expansion
  K <- dim(wave)[2]
  
  ## Special case when rho = 0
  if(rho1==0){
    Sig <- cbind(c(0,0),c(0,0))
  }
  
  ## Compute Sigma (see pg 4 of r_spate documentation)
  else{
    Sig_Chol <- cbind(c(cos(alpha),-gamma*sin(alpha)),
                      c(sin(alpha),gamma*cos(alpha)))/rho1
    Sig <- solve(t(Sig_Chol)%*%Sig_Chol)
  }
  
  ## Compute -Delta*(k_j^T * Sigma * k_j + zeta)
  H_diag <- -dt*apply(wave*Sig%*%wave,2,sum)-dt*rep(zeta,K)
  
  ## Compute -Delta * mu * k_j
  H_off_diag <- dt*c(mu_x,mu_y)%*%wave
  
  ## First n_cos diagonal entries, corresponding to first n_cos
  ## cosine terms
  G_cos <- exp(H_diag[1:n_cos])
  
  ## Other diagonal entries, corresponding to other cosine terms
  G_1 <- exp(H_diag[cosine_indices]) * cos(H_off_diag[cosine_indices])
  
  ## Off (upper) diagonal entries, corresponding to diagonal entries above
  G_2 <- -exp(H_diag[cosine_indices]) * sin(H_off_diag[cosine_indices])
  
  ## Alternative: more straightforward, but SLIGHTLY slower
  
  ## Compute propagator matrix
  # G_matrix = propagator_matrix(wave,cosine_indices,zeta,rho1,gamma,alpha,mu_x,mu_y,dt=1,n_cos=4)
  
  ## First n_cos diagonal entries, corresponding to first n_cos cosine terms
  # G_cos = diag(G_matrix)[1:n_cos]
  
  ## Other diagonal entries, corresponding to other cosine terms
  # G_1 = diag(G_matrix)[cosine_indices]
  
  ## Off (upper) diagonal entries, corresponding to diagonal entries above
  # G_2 = diag(G_matrix[cosine_indices,cosine_indices+1])
  
  ## Output list of vectors
  G <- list(G_cos=G_cos,G_1=G_1,G_2=G_2)
  return(G)
}
####################################

## ---

######################
# SPDE OBJECT FOR FT #
######################
## -> Generate 'SPDE_FT' objects used for two-dim Fourier transform
## -> Inputs:
##            -- 'n', number of points on each ax_is (x,y)
##            -- 't', number of time points
##            -- 'K', number of Fourier basis functions (default = n^2)

## -> Output: list containing
##            -- 'n', number of points on each ax_is (x,y)
##            -- 't', number of time points 
##            -- 'wave', 2 x n^2 matrix of wave numbers
##            -- 'cosine_indices', vector of indices of cosine terms
##            -- 'n_cos', number of cosine-only terms
##            -- 'FFT_indices', indices used for conversion between complex FFT 
##                and real FT

spde_initialise <- function(n,t,K=n^2){
  
  ## Generate wave_numbers object (contains all wave numbers, and 
  ## cosine indices)
  wave_list <- wave_numbers(n) 
  
  ## Add element to list: number of cosine only terms equals 4
  wave_list$n_cos <- 4
  
  ## Indices of wavenumbers of basis functions: 1,...,n^2
  basis_indices <- 1:(n^2)
  
  ## If number of Fourier basis functions is less than n^2 (the full basis):
  if(K<(n^2)){
    
    ## Sort wave-numbers in ascending order of distance from the origin
    
    ## Compute 'cut' as the [K/n^2]^th sample quantile of the 
    ## set of distances of all wave-numbers from the origin
    cut <- quantile(apply(wave_list$wave,2,vec_norm)
                    [order(apply(wave_list$wave,2,vec_norm))],
                    probs=K/(n^2))
    
    ## Indices of wave-numbers of new set of basis functions = wave-numbers
    ## whose distance to origin is less than 'cut'
    basis_indices <- which(apply(wave_list$wave,2,vec_norm)<=cut)
    
    ## Print warning message if number of basis functions to be used (i.e. 
    ## number of wave-numbers whose distance to the origin is less than 'cut') 
    ## isn't equal to the number specified by the user
    if(length(basis_indices)!=K){
      print(paste("Warning: ",length(basis_indices)," Fourier functions", 
                  " (instead of ",K,") are used since when", 
                  " using only ",K," there is anisotropy in", 
                  " the reduced dimensional basis.",sep=""))
    } 
    
    ## Recompute 'wave_list' using only the required wave-numbers
    wave_list$wave <- wave_list$wave[,basis_indices]
    wave_list$n_cos <- 4-sum(is.na(match(1:4,basis_indices)))
    wave_list$cosine_indices <- 
      match(basis_indices[match(wave_list$cosine_indices,basis_indices)
                        [!is.na(match(wave_list$cosine_indices,
                                      basis_indices))]],basis_indices)
    
    K <- length(basis_indices)
    }
  
  ## Compute FFT indices
  FFT_indices <- FFT_index_complex_and_real(n)
  
  ## Generate 'SPDE_FT' object 
  SPDE_FT <- list(n=n,t=t,wave=wave_list$wave,
                  cosine_indices=wave_list$cosine_indices,
                  n_cos=wave_list$n_cos,basis_indices=basis_indices,
                  FFT_indices=FFT_indices, K = length(basis_indices))
  class(SPDE_FT) <- "SPDE_FT"
  
  ## Output
  return(SPDE_FT)
}

######################

## ---

########################
# SPECTRAL PROPAGATION #
########################
## -> Propagate `alpha_t’ vector in time.
## -> Inputs:
##            -- 'alpha_t', vector of spectral coefficients.
##            -- 'SPDE_FT', SPDE_FT object from `spde_initialise' 
##                (require this or 'n')
##            -- 'n', number of points on each axis (x,y) (require
##                this or 'SPDE_FT')
##            -- 'G_vec', spectral propagator matrix in vector form 
##                (require this or 'param')
##            -- 'param', parameters from the SPDE (require this or 
##                'G_vec')

## -> Output:
##            -- 'alpha_t_plus_1', vector of propagated spectral 
##                coefficients (G*alpha_t)


## NB:
## -- 'spectral_propagate_matrix' provides same functionality, 
##    but with G in matrix form. This fn. is considerably slower.
## -- 'spectral_propagate' is a wrapper function for 'propagate_function'

propagate_function = function(xtp1,xt,G_cos,G_1,G_2,n_cos_sin,n_cos){
  for(i in 1:n_cos){
    xtp1[i] = G_cos[i]*xt[i]
  }
  for(i in 1:n_cos_sin){
    xtp1[n_cos + 2*i - 1] = G_1[i] * xt[n_cos + 2*i - 1] + G_2[i] * xt[n_cos + 2*i]
    xtp1[n_cos + 2*i] = G_1[i] * xt[n_cos + 2*i] - G_2[i] * xt[n_cos + 2*i - 1]
  }
  return(xtp1)
}

spectral_propagate <- function(alpha_t,SPDE_FT=NULL,n=NULL,G_vec=NULL,
                               param=NULL){
  
  ## If SPDE_FT not given: create SPDE_FT object
  if(is.null(SPDE_FT)){
    SPDE_FT <- spde_initialise(n=n,t=1)
  }
  
  ## If G_vec not given:
  if(is.null(G_vec)){
    
    ## If also param. not given: report error message
    if(is.null(param)){
      print("Either 'G_vec' or 'param' is required.")
      return()
    }

    ## Else create G (vector form) using input param.
    else{
      G_vec <- propagator_vec(wave=SPDE_FT$wave,
                              cosine_indices=SPDE_FT$cosine_indices,
                              zeta=param[3],rho1=param[4],gamma=param[5],
                              alpha=param[6],mu_x=param[7],mu_y=param[8],
                              dt=1,n_cos=SPDE_FT$n_cos)
    }
  }
  
  ## Compute alpha_t_plus_1 = G*alpha_t
  #alpha_t_plus_1 <- propagate_function(xtp1=(rep(0,SPDE_FT$n*SPDE_FT$n)),
  #                                     xt = alpha_t, 
  #                                     G_cos = G_vec$G_cos,
  #                                     G_1 = G_vec$G_1,
  #                                     G_2 = G_vec$G_2,
  #                                     n_cos_sin = length(SPDE_FT$cosine_indices),
  #                                     n_cos = SPDE_FT$n_cos)
  
  ## Alternative: compute alpha_t_plus_1 = G*alpha_t
  alpha_t_plus_1 <- .C("propagate_spectral", 
                       xtp1 = as.double(rep(0,SPDE_FT$n*SPDE_FT$n)),
                       as.double(alpha_t),as.double(G_vec$G_cos),
                       as.double(G_vec$G_1),as.double(G_vec$G_2),
                       as.integer(length(SPDE_FT$cosine_indices)),
                       as.integer(SPDE_FT$n_cos))$xtp1
                       
  return(alpha_t_plus_1)
}

spectral_propagate_matrix <- function(alpha_t,SPDE_FT=NULL,
                                      n=NULL,G=NULL,param=NULL){
  
  ## If SPDE_FT not given: create SPDE_FT object
  if(is.null(SPDE_FT)){
    SPDE_FT <- spde_initialise(n=n,t=1)
  }
  
  ## If G not given:
  if(is.null(G)){
    
    ## If also param. not given: report error message
    if(is.null(param)){
      print("Either 'G' or 'param' is required.")
      return()
    }
    
    ## Else create G (matrix form) using input param. 
    else{
      G <- propagator_matrix(wave=SPDE_FT$wave,
                             cosine_indices=SPDE_FT$cosine_indices,
                             zeta=param[3],rho1=param[4],gamma=param[5],
                             alpha=param[6],mu_x=param[7],mu_y=param[8],
                             dt=1,n_cos=SPDE_FT$n_cos)
    }
  }
  
  ## Compute alpha_t_plus_1 = G*alpha_t
  alpha_t_plus_1 <- as.vector(G %*% alpha_t)
  return(alpha_t_plus_1)
}

########################

## ---

###########
# REAL FT #
###########

## -> Compute Fast Fourier Transform of spatial field (in vector form) 
## -> Inputs:
##            -- 'w', a spatial field as a stacked vector of length N=n^2.
##            -- 'n', no. of grid  points on each ax_is (x,y)
##            -- 'inv', indicates whether to compute inverse FFT
##            -- 'FFT_indices', list  containing vectors of natural numbers
##                representing indices to transform between Re and Complex FT

## -> Output:
##            -- 'w_FFT_real', vector of real (inverse) FFT of 'w'

FFT_real <- function(w,n,inv=TRUE,FFT_indices=NULL){
  
  ## If not provided, generate FFT_indices
  if(is.null(FFT_indices)){
    FFT_indices <- FFT_index_complex_and_real(n)
  }
  
  ## Compute Fourier Transform
  w_FFT_real <- .C("real_fft",n=as.integer(n),yh=as.double(w),
                inverse=as.integer(sum(inv)),
                indCos=as.integer(FFT_indices$cosine_indices),
                indW=as.integer(FFT_indices$exp_indices),
                indWCon=as.integer(FFT_indices$exp_indices_conj),
                n_cos_sin=as.integer(length(FFT_indices$cosine_indices)))$yh
  return(w_FFT_real)
}

###########

##############
# ST REAL FT #
##############
## -> Compute 2D Fast Fourier Transform of spatio-temporal field
## -> Inputs:
##            -- 'w', a spatio-temporal field as a stacked vector of length T*n^2, 
##                or a matrix of dimensions T x n^2.
##            -- 'n', no. of grid  points on each spatial axis (x,y)
##            -- 'inv', indicates whether to compute inverse FFT (default = TRUE)
##            -- 'FFT_indices', list  containing vectors of natural numbers
##                representing indices to transform between Re and Complex FT

## -> Output:
##            -- 'w_FFT_real', matrix of real (inverse) FFT of 'w'

FFT_real_ST <- function(w,n,t,inv=TRUE,FFT_indices=NULL){
  if(class(w)=="matrix"){
    w <- STmatrix_to_vec(w)
    mat <- TRUE
  }
  else{
    mat <- FALSE
  }
  if(length(w)!=(n*n*t)){
    print("Error: 'w' needs to be a vector of length n*n*T, 
          or a matrix of dimension T x n*n.")
    return()
  }
  if(is.null(FFT_indices)){
    FFT_indices <- FFT_index_complex_and_real(n)
  }
  w_FFT_real <- .C("TSreal_fft",n=as.integer(n),T=as.integer(t),
                   yh=as.double(w),inverse=as.integer(sum(inv)), 
                   indCos=as.integer(FFT_indices$cosine_indices), 
                   indW=as.integer(FFT_indices$exp_indices), 
                   as.integer(FFT_indices$exp_indices_conj), 
                   as.integer(length(FFT_indices$cosine_indices)))$yh
  if(mat){
    w_FFT_real <- vec_to_STmatrix(w_FFT_real,t=t)
  }
  
  return(w_FFT_real)
}
##############

## ---

#################################################
# CONVERT COORDS INDICES TO OBS. MATRIX INDICES #
#################################################
## -> Convert a list of indices of coordinates on a n x n grid, to the 
##    corresponding list of matrix indices for use in (e.g.) the Kalman Filter.
## -> This is particularly useful when n_obs<n^2: in this case, indices of 
##    coordinates can range from 1:n^2, but matrix indices are from 1:n_obs.

## -> Specifically, this function converts a list of n_obs indices, from 1:n^2,
##    to the list of the ranks of these indices, from 1:n_obs

## -> E.g... ind = list(c(1,7),c(4,8)) -> obs_mat_ind = list(c(1,3),c(2,4))

## -> Inputs:
##            -- 'indices', a list (of any length) of vectors of indices in the
##                range 1:n^2

## -> Output:
##            -- 'obs_mat_ind', a list (of the same length as 'indices') of 
##                vectors of indices in the range 1:n_obs

obs_mat_indices = function(indices){
  
  ## Convert indices to their corresponding rank from 1:n_obs ##
  if(!is.null(indices)){
    obs_mat_ind = 
      lapply(1:length(indices),function(i) match(indices[[i]], 
                                                 sort(unlist(indices))))
  }
  else{
    obs_mat_ind = NULL
  }
  
  ## Output ##
  return(obs_mat_ind)
}

#################################################

##########################
# OBSERVATION PARAMETERS #
##########################
## -> Generate a list containing the number of sensor classes (m), the variance
##    of each of the sensor classes (tau2_vals), and the indices of sensors in
##    each of the sensor classes
## -> Inputs:
##            -- 'm', the number of sensor classes (default = 1)
##            -- 'm_indices', a list of length m, containing the indices which
##                correspond to the spatial locations of each of the m sensor
##                classes 
##            -- 'tau2_vals', a list of length m, containing the variances of
##               of each of the m sensor classes
##            -- 'n_obs_indices', a boolean which indicates whether or not
##               to output the matrix conversion of m_indices, as computed by
##               the function obs_mat_indices, or the standard version

## -> Output:
##            -- 'param_obs', a list containing m, tau2_vals, m_indices

param_obs_func <- function(m,m_indices,tau2_vals,n_obs_indices=TRUE){
  if(n_obs_indices){
    param_obs = list(m=m,m_indices=obs_mat_indices(m_indices),tau2_vals=tau2_vals)
  }
  else{
    param_obs = list(m=m,m_indices=m_indices,tau2_vals=tau2_vals)
  }
  
  return(param_obs)
}

##########################

##############################
# OBSERVATION NOISE (MATRIX) #
##############################
## -> Generate observation noise matrix for n x n grid of observations
## -> Inputs:
##            -- 'n', no. of grid  points on each spatial axis (x,y)
##            -- 'param_obs', list containing elements which uniquely determine
##                the observation noise, as generated by 'param_obs_func' 

## -> Output:
##            -- 'obs_noise', the observation noise matrix

obs_noise_mat <- function(n_obs=n^2,param_obs){
  obs_noise <- rep(0,n_obs)
  for (i in 1:(param_obs$m)){
    obs_noise[param_obs$m_indices[[i]]]=param_obs$tau2_vals[[i]]
  }
  return(diag(obs_noise,nrow=n_obs,ncol=n_obs))
}

##############################

##############################
# OBSERVATION NOISE (VECTOR) #
##############################
## -> Generate observation noise matrix for n x n grid of observations
## -> Inputs:
##            -- 'n', no. of grid  points on each spatial axis (x,y)
##            -- 'param_obs', list containing elements which uniquely determine
##                the observation noise, as generated by 'param_obs_func' 

## -> Output:
##            -- 'obs_noise', the diagonal of the observation noise matrix

obs_noise_vec <- function(n_obs=n^2,param_obs){
  obs_noise <- rep(0,n_obs)
  for (i in 1:(param_obs$m)){
    obs_noise[param_obs$m_indices[[i]]]=param_obs$tau2_vals[[i]]
  }
  return(obs_noise)
}

##############################

#################################
# OBSERVATION PARAMETERS (BIAS) #
#################################
## -> Generate a list containing the number of sensor classes (m), the bias
##    of each of the sensor classes (bias_vals), and the indices of sensors in
##    each of the sensor classes
## -> Inputs:
##            -- 'm', the number of sensor classes (default = 1)
##            -- 'm_indices', a list of length m, containing the indices which
##                correspond to the spatial locations of each of the m sensor
##                classes 
##            -- 'bias_vals', a list of length m, containing the biases of
##               of each of the m sensor classes
##            -- 'n_obs_indices', a boolean which indicates whether or not
##               to output the matrix conversion of m_indices, as computed by
##               the function obs_mat_indices, or the standard version

## -> Output:
##            -- 'param_bias', a list containing m, bias_vals, m_indices

param_bias_func <- function(m=0,m_indices=NULL,bias_vals=0,n_obs_indices=TRUE){
  if(n_obs_indices){
    param_bias = list(m=m,m_indices=obs_mat_indices(m_indices),bias_vals=bias_vals)
  }
  else{
    param_bias = list(m=m,m_indices=m_indices,bias_vals=bias_vals)
  }
  return(param_bias)
}

#################################

#############################
# OBSERVATION BIAS (MATRIX) #
#############################
## -> Generate t x n^2 observation bias matrix 
## -> Inputs:
##            -- 'n', no. of grid  points on each spatial axis (x,y)
##            -- 't', no. of time points
##            -- 'param_bias', list containing elements which uniquely determine
##                the observation bias, as generated by 'param_bias_func' 

## -> Output:
##            -- 'obs_bias', the t x n^2 observation bias matrix

obs_bias_mat <- function(n_obs=n^2,t,param_bias){
  
  if(param_bias$m==0){
    obs_bias = matrix(0,nrow=t,ncol=n_obs)
  }
  else{
    obs_bias_row <- rep(0,n_obs)
    for (i in 1:(param_bias$m)){
      obs_bias_row[param_bias$m_indices[[i]]]=param_bias$bias_vals[[i]]
    }
    obs_bias <- matrix(rep(obs_bias_row,t),nrow=t,byrow=TRUE)
  }
  
  return(obs_bias)
}

#############################

############################
# TIME DEPENDENT PARAMETER #
############################
## -> Generate a list containing the initial times - 1=t_0,t_1,t_2,... -  for 
##    which a parameter takes new values - θ_0,θ_1,θ_2,... - and those values
## -> That is to say:
##    - θ = θ_0 for t∈[t_0,t_1-1]
##    - θ = θ_1 for t∈[t_1,t_2-1]
##    - θ = θ_2 for t∈[t_2,t_3-1] 
##    - ...

## -> Inputs:
##            -- 't_vals', a list of the initial times (as above)
##            -- 'param_vals', a list of the corresponding parameter values (or
##                a list of parameter lists generated by param_obs_func or
##                param_bias_func)

## -> Output:
##            -- 'param_t_dep', a list of length two containing the list of
##                times, and the list of param values

param_t_dep_func = function(t_vals,param_vals){
  if(length(t_vals)!=length(param_vals)){
    cat("Error: the number of t-values must be the same of the number of param-values")
  }
  param_t_dep = list(t_vals=t_vals,param_vals=param_vals)
  return(param_t_dep)
}

############################

#############################
# TIME DEPENDENT PARAMETERS #
#############################

## -> Generate a list containing the initial times - 1=t_0,t_1,t_2,... -  for 
##    which the signal, noise, or bias parameters takes new values - θ_0,θ_1,θ_2,...

## -> Inputs:
##            -- 'rho0_t_dep', a list of initial times & parameter values for 
##                rho0, as generated by param_t_dep_func
##            -- 'sigma2_t_dep', a list of initial times & parameter values 
##                for sigma2, as generated by param_t_dep_func
##            -- ...

##            -- 'param_obs_t_dep', a list of initial times & parameter lists
##                for the observation noise, as generated by param_t_dep_func

##            -- 'param_bias_t_dep', a list of initial times & parameter lists
##                for the observation bias, as generated by param_t_dep_func

## -> Output:
##            -- 'param_sig_t_dep', a list of length two containing the list of
##                times, and the list of all signal parameter values

##            -- 'param_obs_t_dep', a list of length two containing the list of
##                times, and the list of observation parameter lists

##            -- 'param_bias_t_dep', a list of length two containing the list of
##                times, and the list of bias parameter lists
 

param_sig_t_dep_func = function(rho0_t_dep,sigma2_t_dep,zeta_t_dep,rho1_t_dep,
                                gamma_t_dep,alpha_t_dep,mu_x_t_dep,mu_y_t_dep){
  
  all_param_vals = list(unlist(rho0_t_dep$param_vals),unlist(sigma2_t_dep$param_vals),
                        unlist(zeta_t_dep$param_vals),unlist(rho1_t_dep$param_vals),
                        unlist(gamma_t_dep$param_vals),unlist(alpha_t_dep$param_vals),
                        unlist(mu_x_t_dep$param_vals),unlist(mu_y_t_dep$param_vals))
  n_param = length(all_param_vals)
  
  all_t_vals = list(unlist(rho0_t_dep$t_vals), unlist(sigma2_t_dep$t_vals),
                    unlist(zeta_t_dep$t_vals), unlist(rho1_t_dep$t_vals),
                    unlist(gamma_t_dep$t_vals), unlist(alpha_t_dep$t_vals),
                    unlist(mu_x_t_dep$t_vals), unlist(mu_y_t_dep$t_vals))
  
  all_t_vals_vec = unlist(all_t_vals)
  
  t_vals = all_t_vals_vec[!duplicated(all_t_vals_vec)]
  t_vals = sort(t_vals)
  t_vals = as.list(t_vals)
  n_t = length(t_vals)
  
  param_vals = rep(list(rep(0,n_param)),n_t)
  for (i in 1:n_param){
    for (j in 1:n_t){
      if(t_vals[[j]]%in%all_t_vals[[i]]){
        t_index = which(all_t_vals[[i]]==t_vals[[j]])
      }
      else{
        t_index = max(which(all_t_vals[[i]]<t_vals[[j]]))
      }
      param_vals[[j]][i] = all_param_vals[[i]][t_index]
    }
  }
  
  param_sig_t_dep = list(t_vals=t_vals,param_vals=param_vals)
  
  return(param_sig_t_dep)
}

param_obs_t_dep_func = function(param_obs_t_dep){
  return(param_obs_t_dep)
}

param_bias_t_dep_func = function(param_bias_t_dep){
  return(param_bias_t_dep)
}

#############################


## ---

###################
# SPDE SIMULATION #
###################
## -> Generate a single sample from the GP specified through the SPDE
## -> Inputs:
##            -- 'n', no. of grid  points on each spatial axis (x,y)
##            -- 't', no. of time points
##            -- 'param_sig_t_dep', possibly time-dependent signal parameters,
##                as generated by param_sig_t_dep_func 
##            -- 'param_obs_t_dep', possibly time-dependent observation noise
##                parameters, as generated by param_obs_t_dep_func 
##            -- 'param_bias_t_dep', possibly time-dependent observation bias
##                parameters, as generated by param_bias_t_dep_func 
##            -- 'nu', smoothness parameter of Matern cov. fn. for 
##                innovations (default = 1 <-> Whittle cov. fn.)
##            -- 'initial_value', starting value (field) for SPDE
##            -- 'seed', seed for random number generation
##            -- 'dt', time step
##            -- 'n_obs', number of observations
##            -- 'obs_loc', a list of length n_obs, whose j^th element is a 
##                t x 2 matrix containing the (x,y) coords of the j^th obs. 
##                location for each time.
##            -- 'alpha_sim', a boolean indicating whether to simulate the
##                Fourier coefficients alpha
##            -- 'x_sim', a boolean indicating whether to simulate the latent
##                field x
##            -- 'obs_noise_sim', a boolean indicating whether to simulate the
##                observation noise
##            -- 'obs_bias_sim', a boolean indicating whether to simulate the
##                observation bias
##            -- 'y_sim', a boolean indicating whether to simulate the
##                observations y

## -> Output:
##            -- 'simulation', depending on the inputs, a list containing
##                α, x, η, β, y, param_sig, param_obs, param_bias, n, t, dt

## We also provide auxiliary functions to format the print / summary of this
## function: 'print.SPDEsimulation' and 'summary.SPDEsimulation'; and to plot
## the output: 'plot.SPDEsimulation'.

spde_simulation <- function(n, t, param_sig_t_dep, param_obs_t_dep,
                            param_bias_t_dep, nu = 1, initial_value = NULL,  
                            seed = NULL, dt = 1, n_obs = NULL, obs_loc = NULL,
                            alpha_sim = TRUE, x_sim = TRUE, 
                            obs_noise_sim = TRUE, obs_bias_sim = TRUE,
                            y_sim = TRUE, weights = NULL){
  
  ## ------------- ##
  ## Preliminaries ##
  
  ## If n_obs not provided, set n_obs = n^2 ##
  if(is.null(n_obs)){
    n_obs = n^2
  }
  
  ## SPDE_FT object ##
  SPDE_FT <- spde_initialise(n=n,t=t)
  wave = SPDE_FT$wave
  cosine_indices = SPDE_FT$cosine_indices
  n_cos = SPDE_FT$n_cos
  
  ## Output ##
  simulation = list()
  
  ## ------------- ##
  
  
  ## ---------------------------- ##
  ## Compute α_t [T x n^2 matrix] ##
  if(x_sim || y_sim){
    alpha_sim = TRUE
  }
  if(alpha_sim){
    
    ## Initialise alpha ##
    alpha = matrix(0,nrow=t,ncol=n^2)
    
    ## Number of time epochs (for signal parameters)
    n_t = length(param_sig_t_dep$t_vals)
    
    ## Set seed
    if(!is.null(seed)){set.seed(seed)} 
    
    ## Simulate α_t in each epoch
    for (i in 1:n_t){
      
      ## Initial time of this epoch
      t_init = param_sig_t_dep$t_vals[[i]]
      
      ## Final time of this epoch
      if(i<n_t){t_fin = param_sig_t_dep$t_vals[[i+1]]-1}
      if(i==n_t){t_fin = t}
      
      ## Length of this epoch
      t_tmp = t_fin + 1 - t_init
      
      ## Signal parameters
      param_sig = param_sig_t_dep$param_vals[[i]]
      
      ## Spectral density f(k_j)
      spectrum <- innovation_spectrum(wave=wave,n=n,n_cos=n_cos,
                                      rho0=param_sig[1],sigma2=param_sig[2],
                                      zeta=param_sig[3],rho1=param_sig[4],
                                      gamma=param_sig[5],alpha=param_sig[6],
                                      nu=nu,dt=dt,norm=TRUE)
      
      ## Propagator matrix G (vector form)
      G_vec <- propagator_vec(wave=wave,cosine_indices=cosine_indices,
                              zeta=param_sig[3],rho1=param_sig[4],
                              gamma=param_sig[5],alpha=param_sig[6],
                              mu_x=param_sig[7],mu_y=param_sig[8],
                              dt=dt,n_cos=n_cos)
      
      ## t_tmp x n^2 matrix of innovations ε_t
      innovation <- vec_to_STmatrix(rep(sqrt(spectrum),t_tmp)*rnorm(n*n*t_tmp),
                                    t = t_tmp)
      
      ## t_tmp x n^2 matrix of weighted innovations 
      if(!is.null(weights)){
        innovation <- FFT_real_ST(t(t(FFT_real_ST(innovation,n=n,t=t_tmp,inv=FALSE))*weights),n=n,t=t_tmp,inv=TRUE)
        #innovation <- FFT_real_ST(sweep(FFT_real_ST(innovation,n=n,t=t_tmp,inv=FALSE),MARGIN=2,weights,'*'),n=n,t=t_tmp,inv=TRUE)
      }
      
      ## Initial value
      if(i==1){
        if(!is.null(initial_value)){
          alpha[t_init,] <- FFT_real(initial_value,n,inv=TRUE)
        }
        else{
          alpha[t_init,] <- innovation[t_init,]
        }
      }
      
      ## Compute α_t = G*α_{t-1} + epsilon_t
      if (i>1){
        alpha[t_init,] <- spectral_propagate(alpha_t = alpha[t_init-1,], 
                                             SPDE_FT = SPDE_FT,
                                             G_vec = G_vec) + innovation[1,]
      }
      
      if(t_tmp>1){
        for(j in (t_init+1):(t_fin)){
          alpha[j,] <- spectral_propagate(alpha_t = alpha[j-1,], SPDE_FT = SPDE_FT,
                                          G_vec = G_vec) + innovation[j+1-t_init,]
        }
      }
    }
    
    ## Output
    simulation = c(simulation,list(alpha = alpha))
    
  }
  ## ---------------------------- ##
  
  
  ## ------------------------------ ##
  ## Compute x_t [T x n_obs matrix] ##
  if(y_sim){
    x_sim = TRUE
  }
  if(x_sim){
    
    ## If observation locations not provided [so observation locations assumed
    ## to be a unit n x n grid], compute x using FFT
    if(is.null(obs_loc)){
      x <- FFT_real_ST(alpha,n=n,t=t,inv=FALSE)
    }
    
    ## If observation locations provided, compute x = Phi*α
    else{
      
      ## Initialise x_t ##
      x = matrix(0,nrow=t,ncol=n_obs)
      
      ## Compute ##
      for (i in 1:t){
        obs_loc_i = lapply(obs_loc,function(x) as.numeric(x[i,]))
        x[i,] = FFT_real_matrix(wave,cosine_indices,n_cos,n,n_obs,
                                obs_loc_i)%*%alpha[i,]
      }
    }
    
    ## Output ##
    simulation = c(simulation,list(x = x))
    
  }
  ## ------------------------------ ##
  
  
  ## Compute η_t [T x n_obs matrix] ##
  if(y_sim){
    obs_noise_sim = TRUE
  }
  if(obs_noise_sim){
    
    ## Initialise η_t ##
    obs_error = matrix(0,nrow=t,ncol=n_obs)
    
    ## Number of time epochs (obs. noise parameters)
    n_t = length(param_obs_t_dep$t_vals)
    
    ## Compute η_t in each epoch
    for (i in 1:n_t){
      
      ## Initial time of this epoch
      t_init = param_obs_t_dep$t_vals[[i]]
      
      ## Final time of this epoch
      if(i<n_t){t_fin = param_obs_t_dep$t_vals[[i+1]]-1}
      else{t_fin = t}
      
      ## Length of this epoch
      t_tmp = t_fin + 1 - t_init
      
      ## Observation noise
      obs_noise = obs_noise_vec(n_obs,param_obs_t_dep$param_vals[[i]])
      
      ## Observation error
      obs_error[t_init:t_fin,] = 
        vec_to_STmatrix(rep(sqrt(obs_noise),t_tmp)*rnorm(n_obs*t_tmp),t=t_tmp)
    }
    
    ## Output ##
    simulation = c(simulation,list(obs_error = obs_error))
    
  }
  ## ------------------------------ ##
  
  
  ## ------------------------------ ##
  ## Compute β_t [T x n_obs matrix] ##
  if(y_sim){
    obs_bias_sim = TRUE
  }
  if(obs_bias_sim){
    ## Initialise β_t ##
    obs_bias = matrix(0,nrow=t,ncol=n_obs)
    
    ## Number of time epochs (obs. bias parameters)
    n_t = length(param_bias_t_dep$t_vals)
    
    for (i in 1:n_t){
      
      ## Initial time of this epoch
      t_init = param_bias_t_dep$t_vals[[i]]
      
      ## Final time of this epoch
      if(i<n_t){t_fin = param_bias_t_dep$t_vals[[i+1]]-1}
      else{t_fin = t}
      
      ## Length of this epoch
      t_tmp = t_fin + 1 - t_init
      
      ## Observation noise
      obs_bias[t_init:t_fin,] = 
        obs_bias_mat(n_obs,t_tmp,param_bias_t_dep$param_vals[[i]])
    }
    
    ## Output ##
    simulation = c(simulation,list(obs_bias = obs_bias))
    
  }
  ## ------------------------------ ##
  
  
  ## ------------------------------ ##
  ## Compute y_t [T x n_obs matrix] ##
  if(y_sim){
    y <- x + obs_bias + obs_error
    
    ## Output ##
    simulation = c(simulation, list(y = y))
    
  }
  ## ------------------------------ ##
  
  ## Output ##
  simulation <- c(simulation,list(param_sig_t_dep = param_sig_t_dep,
                                  param_obs_t_dep = param_obs_t_dep,
                                  param_bias_t_dep = param_bias_t_dep, 
                                  n = n, t = t, dt = dt))
  
  class(simulation) <- "SPDEsimulation"
  
  return(simulation)
}

print.SPDEsimulation  <- function(x,...){
  print(paste("SPDEsimulation object with t=",x$t,", dt=",x$dt,", n=",x$n,
              sep=""))
}

summary.SPDEsimulation  <- function(x,...){
  print(paste("SPDEsimulation object with t=",x$t,", dt=",x$dt,", n=",x$n,
              sep=""))
}

plot.SPDEsimulation <- function(x,...,plot_x=TRUE,plot_y=FALSE){
  if(plot_x) SPDE_plot(x_i=x$x,dt=x$dt,...)
  if(plot_y) SPDE_plot(x_i=x$y,dt=x$dt,...)
}

###################

#############
# SPDE PLOT #
#############
## -> Generate plots of ST field
## -> Inputs:
##            -- 'x_i', ST field represented as a TxN matrix
##            -- 'n_x', number (integer) of points on the x axis
##            -- 'which_t', vector of time points to be plotted
##            -- 'format', string specifying the format of the final
##                plot (="images_together" or "images_separate")
##            -- 'save_file', whether or not to save file
##            -- 'path', string specifying where to save file
##            -- 'file', file name
##            -- 'abs_scale', whether colour scale for all plots
##                is the same
##            -- 'main', title for plots (must be of length 1, or the
##                same length as 'which_t')
##            -- 'mfrow', see ?par 
##            -- 'image_size', size of .jpeg is 'save_file'=TRUE
##            -- 'z_lim', colour scale limits of plots, see ?image
##            -- 'breaks', colour scale breaks of plots, see ?image
##            -- '...', any other graphical parameters to be passed to 
##               'image' and 'plot'

##            -- NB: this fn should be used for SPDEsimulations with 
##               n_obs = n^2 only

## -> Output:
##            -- plot(s) of the spatio-temporal field

SPDE_plot <- function(x_i, n_x = NULL, which_t = NULL, dt = 1,
                      format = "images_together", save_file = FALSE, 
                      path = NULL, file = NULL, abs_scale = FALSE , 
                      main = NULL, mfrow = NULL, image_size = c(1000, 1000), 
                      zlim = NULL, breaks = NULL,...){
  
  ## If 'which_t' not specified, plot all time points
  if(is.null(which_t)){
    which_t <- 1:dim(x_i)[1]
  }
  
  ## If 'zlim' not specified, define as min. & max. of values of x_i
  if(is.null(zlim)){
    zlim <- c(min(x_i[which_t,]), max(x_i[which_t,]))
  }
  
  ## If 'breaks' not specified, define as an equally partitioned sequence 
  ## from zlim[1] to zlim[2], with length  equal to length of cols() vector
  if(is.null(breaks)){
    breaks = seq(from = zlim[1],to = zlim[2], 
                 length.out = length(colors())+1)
  }
  
  ## If 'n_x' not specified, define n_x = n_y = dim(x_i)[2]..
  if(is.null(n_x)){
    n_x <- n_y <- sqrt(dim(x_i)[2])
  }
  ## If 'n_x' is specified, define n_y = dim(x_i)[2]/n_x. 
  else{
    n_y <- dim(x_i)[2]/n_x
  }
  
  ## Define x,y axes
  x_axis = 1/n_x*(0:(n_x-1))
  y_axis = 1/n_y*(0:(n_y-1))
  
  ## Define grid
  grid_xy = expand.grid(x_axis,y_axis)
  
  ## Define plotting function
  plot.x_i <- function(){
    
    level_plots = list()
    
    ## Loop over time points
    for(i in which_t){
      
      ## Plot title(s) ##
      
      ## If 'main' not specified, specify title as "t=...'
      if(is.null(main)){
        maint <- paste("t = ",i*dt,sep="")
      }
      
      ## If 'main' specified, and same length as 'which_t', set title of  
      ## plot i equal to i^th element of 'main'
      else if(length(main)==length(which_t)){
        maint <- main[i]
      }
      
      ## If 'main' specified, and of length 1, set title of all plots to
      ## that title
      else if(length(main)==1){
        maint <- main
      }
      
      ## Otherwise, write error message
      else print("'main' must either be NULL or a 
                 character vector of length equal 
                 to the number of time points or 1.")
      
      ## Other options ##
      
      ## If format="images_separate" and save_file = TRUE, open jpeg graphics 
      ## device at each iteration of the loop
      if(format=="images_separate" & save_file){
        jpeg(paste(path,file,t,".jpeg",sep=""), width = image_size[1], 
             height = image_size[2])
      }
      
      ## If abs_scale = TRUE, plot each image with its own colour scale
      if(abs_scale){
        image(x_axis,y_axis,matrix(x_i[i,],nrow=n_x),col=colors(),
              main=maint,xlim=c(0,1),ylim=c(0,1),xlab="x",ylab="y",
              axes=F,xaxt="n",yaxt="n",...)
        axis(1,at=seq(0,1,0.5),labels=as.character(seq(0,1,0.5)))
        axis(2,at=seq(0,1,0.5),labels=as.character(seq(0,1,0.5)))
        axis(3, lwd.tick=0, labels=FALSE); axis(4, lwd.tick=0, labels=FALSE);
      } 
      
      ## If abs_scale = FALSE, plot each image on the same colour scale
      else{
        image(x_axis,y_axis,matrix(x_i[i,],nrow=n_x),col = colors(),
                 zlim=zlim,breaks=breaks,main=maint,xlim=c(0,1),ylim=c(0,1),
                 xlab="x",ylab="y",xaxt="n",yaxt="n",...)
        axis(1,at=seq(0,1,0.5),labels=as.character(seq(0,1,0.5)))
        axis(2,at=seq(0,1,0.5),labels=as.character(seq(0,1,0.5)))
        axis(3, lwd.tick=0, labels=FALSE); axis(4, lwd.tick=0, labels=FALSE);
      }
      
      ## If format = "images_separate", and save_file = TRUE, turn graphics 
      ## device off before next iteration
      if(format=="images_separate" & save_file){
        dev.off()
      }
    }
  }
  
  ## Now use this function to plot:
  
  ## If format = 'images_together': 
  if(format=="images_together"){
    
    ## If save_file = TRUE, open jpeg graphics device
    if(save_file) jpeg(paste(path,file,".jpeg",sep=""), width = imagesize[1], 
                       height = imagesize[2])
    
    ## If 'mfrow' not specified, define according to length of 'which_t'
    if(is.null(mfrow)){
      if(length(which_t)==1) mfrow=c(1,1)
      else if(length(which_t)==2) mfrow=c(1,2)
      else if(length(which_t)<=4) mfrow=c(2,2)
      else if(length(which_t)<=6) mfrow=c(2,3)
      else if(length(which_t)<=9) mfrow=c(3,3)
      else if(length(which_t)<=12) mfrow=c(3,4)
      else if(length(which_t)<=16) mfrow=c(4,4)
      else mfrow=c(5,5)
    }
    
    ## Set number of plot rows
    par(mfrow=mfrow,...)
    
    ## Plot!
    plot.x_i()
    
    ## If save_file = TRUE, close graphics device 
    if(save_file) dev.off()
  }
  
  ## If format = 'images_separate': 
  else if(format=="images_separate"){
    
    ## Plot!
    plot.x_i()
  }
}

#############

################
# Y SIMULATION #
################
## -> Generate observations from the latent field at 'n_obs' observation 
##    locations 'obs_loc', and at times 't_vals', given an existing simulation 
##    of the Fourier coefficients, the observation noise, and the observation 
##    bias, over a time period including 't_vals'
##    
## -> Inputs:
##            -- 'spde_sim', a simulation of the latent field, generated by
##                spde_simulation (for efficiency, set x_sim = FALSE and 
##                y_sim = FALSE)
##            -- 't_vals', the time points at which to generate the observations
##            -- 'n_obs', the number of observation locations.
##            -- 'obs_loc', a list of length n_obs, whose j^th element is a 
##                matrix of dimensions length(t_vals) x 2, containing the (x,y) 
##                coords of the j^th obs location at times t_vals; ALTERNATIVELY,
##                if 'obs_loc_t_dep' = FALSE, a list of length n_obs, whose
##                j^th element is a vector of length 2, containing the constant
##                (x,y) corods of the j^th obs location at times t-vals
##            -- 'obs_loc_t_dep', a boolean input indicating whether or not the
##                observation locations are time dependent 

y_simulation <- function(spde_sim, t_vals, n_obs, obs_loc,
                         obs_loc_t_dep = FALSE){
  
  ## Compute n ##
  n = sqrt(dim(spde_sim$alpha)[2])
  
  ## SPDE_FT object ##
  SPDE_FT <- spde_initialise(n=n,t=t)
  wave = SPDE_FT$wave
  cosine_indices = SPDE_FT$cosine_indices
  n_cos = SPDE_FT$n_cos
  
  ## Extract alpha,obs_bias,obs_error at time t ##
  alpha_t = matrix(spde_sim$alpha[t_vals,],nrow=length(t_vals))
  obs_bias_t = matrix(spde_sim$obs_bias[t_vals,],nrow=length(t_vals))
  obs_error_t = matrix(spde_sim$obs_error[t_vals,],nrow=length(t_vals))
  
  ## If observation locations are CONSTANT ##
  if(!obs_loc_t_dep){
    
    ## Compute x_t ##
    x_t = t(FFT_real_matrix(wave,cosine_indices,n_cos,n,n_obs,
                            obs_loc)%*%t(alpha_t))
  }
  
  ## If observation locations are TIME VARYING ##
  if(obs_loc_t_dep){
    
    ## Initialise x_t ##
    x_t = matrix(0,nrow=length(t_vals),ncol=n_obs)
    
    ## Compute x_t ##
    for (i in 1:length(t_vals)){
      obs_loc_i = lapply(obs_loc,function(x) x[i,])
      x_t[i,] = FFT_real_matrix(wave,cosine_indices,n_cos,n,n_obs,
                                obs_loc_i)%*%alpha_t[i,]
    }
  }
  
  ## Compute y_t ##
  y_t = x_t + obs_bias_t + obs_error_t
  
  ## Output ##
  return(y_t = y_t)
}

################


## ---

#################
# KALMAN FILTER #
#################
## -> Kalman Filter for the Gaussian state space model
## -> Inputs:
##            -- 'y', observed data in t x N matrix
##            -- 'linear_pred', linear predictor in t x N matrix
##            -- 'Phi', observation matrix relation y to alpha
##            -- 'R', covariance of observation error nu
##            -- 'G', propagation matrix
##            -- 'Q', innovation covariance matrix of alpha
##            -- 'N', number of spatial points
##            -- 't', number of time points
##            -- 'K', dimension of latent process alpha
##            -- 'log_lik', whether to also return log-likelihood
##            -- 'm0', initial mean
##            -- 'RO', initial covariance

## -> Output: list containing
##            -- 'm_i_i', (t+1) x N matrix of filtered values of alpha (first row
##                is equal to m_0)
##            -- 'R_i_i', (t+1) x N x N array of filtered covariances of alpha
##                (first row is equal to R0)
##            -- 'm_i_i1', t x N matrix of forecast values of alpha
##            -- 'R_i_i1', t x N x N array of forecast covariances of alpha
##            -- 'innov_resid', t x N matrix of innovation residuals
##            -- 'innov_cov', t x N X N array of innovation residual covariances
##            -- 'K_i', N x N array containing Kalman Gain at time T
##            -- 'll', log-likelihood

## The model is in the following form:
# --- y_t = linpred_t + Phi*alpha_t + nu_t; nu_t ~ N(0,R)
# --- alpha_t = G*alpha_{t-1} + epsilon_t; epsilon_t ~ N(0,Q)

KalmanFilter <- function(y, linear_pred = matrix(0, dim(y)[1], dim(y)[2]), Phi,
                         R, G, Q, N = dim(y)[2], t = dim(y)[1], K = dim(G)[1],
                         log_lik = TRUE, m0 = rep(0,K), R0 = as.matrix(Q)){
  
  ## Preliminaries
  Phi_T <- t(Phi)
  
  ## Initialise various arrays
  
  # Predictions
  m_i_i1 <- array(0, c(t, K)) # Mean
  R_i_i1 <- array(0, c(t, K, K)) # Covariance
  
  # Updates
  m_i_i <- array(0, c(t+1, K)) # Mean
  R_i_i <- array(0, c(t+1, K, K)) # Covariance
  
  # Innovations
  innov_resid <- array(0, c(t, dim(y)[2])) # Residual 
  innov_cov <- array(0, c(t, dim(y)[2], dim(y)[2])) # Covariance
  
  # Initial covariance 
  R_i_i[1, , ] <- R0
  
  # Initial mean
  m_i_i[1,] <- m0
  
  ## Kalman Filter
  for (i in 1:t) {
    
    ## Prediction
    ## -- Mean
    m_i_i1[i, ] <- G %*% m_i_i[i, ] 
    
    ## -- Covariance
    R_i_i1[i, , ] <- Q + G %*% R_i_i[i, , ] %*% t(G)
    R_i_i1[i, , ] <- (R_i_i1[i, , ] + t(R_i_i1[i, , ]))/2 
    
    ## Update 
    ## -- Innovation Covariance
    S_i <- R + Phi %*% R_i_i1[i, , ] %*% Phi_T
    S_i_inv <- solve(S_i)
    innov_cov[i,,] <- S_i_inv
    
    ## -- Optimal Kalman Gain
    K_i <- R_i_i1[i, , ] %*% Phi_T %*% S_i_inv
    
    ## -- Innovation Residual
    innov_resid[i, ] <- y[i, ] - linear_pred[i, ] - Phi %*% m_i_i1[i, ]
    
    ## -- Mean
    m_i_i[i+1, ] <- m_i_i1[i, ] + K_i %*% (innov_resid[i, ]) 
    
    ## -- Covariance
    R_i_i[i+1, , ] <- R_i_i1[i, , ] - K_i %*% Phi %*% R_i_i1[i, , ] 
    R_i_i[i+1, , ] <- (R_i_i[i+1, , ] + t(R_i_i[i+1, , ]))/2
  }
  
  
  ## If 'log_lik'=TRUE: compute likelihood 
  if (log_lik) {
    
    # Initial value
    ll <- 0
    
    # Sequentially compute log-likelihood
    for (i in 1:t) {
      ll <- ll + determinant(as.matrix(innov_cov[i, , ]), logarithm = TRUE)$modulus[[1]] - 
        t(innov_resid[i, ]) %*% innov_cov[i, , ] %*% innov_resid[i,]
    }
    
    # Final value
    ll <- (ll - t*K*log(2*pi))/2
  }
  
  ## Output
  output = list(m_i_i = m_i_i, R_i_i = R_i_i, m_i_i1 = m_i_i1, R_i_i1 = R_i_i1,
                innov_resid = innov_resid, innov_cov = innov_cov, 
                K_i = K_i)
  
  ## Amend output according to if log_lik = TRUE
  if (log_lik){
    output = c(output, list(ll=ll))
  }
  
  ## Final output
  return(output)
}

#################

##########################
# SPECTRAL KALMAN FILTER #
##########################

KF_spectral <- function(w = NULL, w_FFT_real = NULL, innov_spectrum = NULL,
                        G_vec = NULL, tau2 = NULL, param_obs = NULL, 
                        param_sig = NULL, n, t, K=n^2, log_lik = TRUE, 
                        cosine_indices=(1:((n*n-4)/2)*2+3),
                        n_cos=4, nu=1, dt=1, m0 = rep(0,K), R0 = rep(0,K)){
  
  ## If 'w_FFT_real' not given, compute as inverse-FFT of w
  if(is.null(w_FFT_real)){
    w_FFT_real <- FFT_real_ST(STmatrix_to_vec(w),n=n,t=t,inv=TRUE)
  }
  
  ## If 'w_FFT_real' is given:
  else{
    # If 'w_FFT_real' in matrix form, convert into vector
    if(class(w_FFT_real)=="matrix"){
      w_FFT_real <- STmatrix_to_vec(w_FFT_real)
    }
    # If length of 'w_FFT_real' != t x n x n, output error message
    if(length(w_FFT_real)!=(n^2*t)){
      print(paste("Error: 'w_FFT_real' must be a vector of length t*n^2 or a ",
                  "matrix of dimension t x n^2.",sep=""))
      return()
    }
  }
  
  ## Select required terms from w_FFT_real
  if(K<n^2){
    basis_indices <- spde_initialise(n,t,K)$basis_indices
    w_FFT_indices <- rep(seq(0,n^2*(t-1),n^2),each=K)+basis_indices
    w_FFT_real = w_FFT_real[w_FFT_indices]
  }
  
  ## If required inputs not specified, output error message
  if((is.null(G_vec) & is.null(param_sig)) | (is.null(innov_spectrum) & 
                                          is.null(param_sig))){
    print("Either {'G_vec', 'spectrum', 'tau2'} or 'param' must be supplied.")
    return()
  }
  
  ## Compute propagator matrix G (vector form)
  if(is.null(G_vec)){
    SPDE_FT <- spde_initialise(n,t,K)
    wave <- SPDE_FT$wave
    cosine_indices <- SPDE_FT$cosine_indices
    n_cos <- SPDE_FT$n_cos
    G_vec <- propagator_vec(wave=wave,cosine_indices=cosine_indices,
                            zeta=param_sig[3],rho1=param_sig[4],
                            gamma=param_sig[5],alpha=param_sig[6],
                            mu_x=param_sig[7],mu_y=param_sig[8],dt=dt,
                            n_cos=n_cos)
  }
  
  ## Compute innovation spectrum f(k_j) 
  if(is.null(innov_spectrum)){
    if(is.null(SPDE_FT)){
      SPDE_FT <- spde_initialise(n,t,K)
      wave <- SPDE_FT$wave
      cosine_indices <- SPDE_FT$cosine_indices
      n_cos <- SPDE_FT$n_cos
    }
    innov_spectrum <- innovation_spectrum(wave=wave,n=n,n_cos=n_cos,
                                          rho0=param_sig[1],sigma2=param_sig[2],
                                          zeta=param_sig[3],rho1=param_sig[4],
                                          gamma=param_sig[5],alpha=param_sig[6],
                                          nu=nu,dt=dt,norm=TRUE)
  }
  
  ## Compute tau^2 
  if(is.null(tau2)){
    tau2=param_obs$tau2_vals[[1]]
  }
  
  
  ## Run 'kf_spectral'
  KF <- kf_spectral(w_FFT_real = w_FFT_real,m_i_i1 = rep(0,K*t),
                    m_i_i=c(m0,rep(0,K*t)),R_i_i1=rep(0,K*t),
                    R_i_i=c(R0,rep(0,K*t)),G_cos=G_vec$G_cos,G_1=G_vec$G_1, 
                    G_2=G_vec$G_2,spectrum=innov_spectrum,t=t,n_cos=n_cos,
                    n_cos_sin=length(cosine_indices),tau2=tau2)
  
  m_i_i = vec_to_STmatrix(KF$m_i_i,t=t+1)
  m_i_i1 = vec_to_STmatrix(KF$m_i_i1,t=t)
  R_i_i = vec_to_STmatrix(KF$R_i_i,t=t+1)
  R_i_i1 = vec_to_STmatrix(KF$R_i_i1,t=t)
  
  ## Output(s)
  output <- list(m_i_i = m_i_i, R_i_i = R_i_i, m_i_i1 = m_i_i1, R_i_i1 = R_i_i1)
  
  ## If 'log_lik'=TRUE: compute likelihood  & amend output
  if(log_lik){
    ll = ll_spectral(log_lik=0,w_FFT_real, KF$m_i_i1, KF$R_i_i1, t, K, 
                     tau2)$ll
    
    if(K<n^2){
      ll = ll_non_spectral(log_lik=0,w,KF$m_i_i1,KF$R_i_i1,t,n,K,tau2)$ll
    }
    
    output = c(output, list(ll=ll))
  }
  return(output)
}

kf_spectral = function(w_FFT_real,m_i_i1,m_i_i,R_i_i1,R_i_i,G_cos,G_1,G_2,
                       spectrum,t,n_cos,n_cos_sin,tau2){
  
  kf = .C("kf_spectral_2",y_FFT=as.double(w_FFT_real),
          mu_i_i1=as.double(m_i_i1),mu_i_i=as.double(m_i_i),
          Sigma_i_i1=as.double(R_i_i1),Sigma_i_i=as.double(R_i_i), 
          G_cos = as.double(G_cos),G_1 = as.double(G_1),
          G_2 = as.double(G_2),Q = as.double(spectrum), 
          T=as.integer(t),n_cos=as.integer(n_cos),
          n_cos_sin = as.integer(n_cos_sin),tau2=as.double(tau2))
  
  output = list(m_i_i1 = kf$mu_i_i1, m_i_i = kf$mu_i_i, R_i_i1 = kf$Sigma_i_i1,
                R_i_i = kf$Sigma_i_i)
  
  return(output)
}

##########################

## ---

###################
# KALMAN SMOOTHER #
###################
## -> Kalman Smoother for the Gaussian state space model
## -> Inputs:
##            -- 'y', observed data in t x N matrix
##            -- 'linear_pred', linear predictor in t x N matrix
##            -- 'Phi', observation matrix relation y to alpha
##            -- 'R', covariance of observation error nu
##            -- 'G', propagation matrix
##            -- 'Q', innovation covariance matrix of alpha
##            -- 'N', number of spatial points
##            -- 't', number of time points
##            -- 'K', dimension of latent process alpha
##            -- 'log_lik', whether to also return log-likelihood
##            -- 'm0', initial mean
##            -- 'RO', initial covariance

## -> Output: list containing
##            -- 'm_i_t', (t+1) x N matrix of smoothed values of alpha
##            -- 'R_i_t', (t+1) x N x N array of smoothed covariances of alpha
##            -- 'J_i', t x N x N array of 'smoothing matrices'
##            -- 'm_i_i', (t+1) x N matrix of filtered values of alpha
##                (first row = m_0)
##            -- 'R_i_i', (t+1) x N x N array of filtered covariances of alpha
##                (first row is equal to R0)
##            -- 'm_i_i1', t x N matrix of forecast values of alpha
##            -- 'R_i_i1', t x N x N array of forecast covariances of alpha
##            -- 'innov_resid', t x N matrix of innovation residuals
##            -- 'innov_cov', t x N X N array of innovation residual covariances
##            -- 'K_i', N x N array containing Kalman Gain at time T
##            -- 'll', log-likelihood

## The model is in the following form:
# --- y_t = linpred_t + Phi*alpha_t + nu_t; nu_t ~ N(0,R)
# --- alpha_t = G*alpha_{t-1} + epsilon_t; epsilon_t ~ N(0,Q)

KalmanSmoother <- function(y, linear_pred = matrix(0, dim(y)[1], dim(y)[2]), 
                           Phi, R, G, Q, N = dim(y)[2], t = dim(y)[1], 
                           K = dim(G)[1],log_lik = TRUE, m0 = rep(0,K), 
                           R0 = as.matrix(Q)){
  
  ## Kalman Filter
  KF = KalmanFilter(y = y, linear_pred = linear_pred, Phi = Phi, R = R, G = G, 
                    Q = Q, N = N, t = t, K = K, log_lik = log_lik, m0 = m0, 
                    R0 = R0)
  
  m_i_i = KF$m_i_i
  m_i_i1 = KF$m_i_i1
  R_i_i = KF$R_i_i
  R_i_i1 = KF$R_i_i1
  innov_resid = KF$innov_resid
  innov_cov = KF$innov_cov
  K_i = KF$K_i
  
  ## Initialise arrays
  m_i_t = array(0, c(t+1, K))
  R_i_t = array(0, c(t+1, K, K))
  J_i = array(0, c(t, K, K)) 
  
  ## Starting values (output from KalmanFilter)
  m_i_t[t+1,]=m_i_i[t+1,] 
  R_i_t[t+1,,]=R_i_i[t+1,,]
  
  ## Smoothing
  for(i in t:1)  {
    J_i[i,,]=(R_i_i[i,,]%*%t(G))%*%solve(R_i_i1[i,,])
    m_i_t[i,]=m_i_i[i,]+J_i[i,,]%*%(m_i_t[i+1,]-m_i_i1[i,])
    R_i_t[i,,]=R_i_i[i,,]+J_i[i,,]%*%(R_i_t[i+1,,]-R_i_i1[i,,])%*%t(J_i[i,,])
  }
  
  ## Output
  output = list(m_i_t = m_i_t, R_i_t = R_i_t, J_i = J_i, m_i_i = m_i_i, 
                R_i_i = R_i_i, m_i_i1 = m_i_i1, R_i_i1 = R_i_i1, 
                innov_resid = innov_resid, innov_cov = innov_cov, K_i = K_i)
  
  ## Amend output according to if log_lik = TRUE
  if (log_lik){
    output = c(output, list(ll=KF$ll))
  }
  
  ## Final output
  return(output)
}


###################

## ---

####################
# BACKWARD SAMPLER #
####################

bs_spectral = function(alpha_sim,m_i_i,m_i_i1,R_i_i, R_i_i1,spectrum,
                       G_cos,G_1,G_2,t,n_cos_sin,n_cos){
  
  bs = .C("bs_spectral",simAlpha=as.double(alpha_sim),
          mtt=as.double(m_i_i),mtt1=as.double(m_i_i1),
          rtt=as.double(R_i_i),rtt1=as.double(R_i_i1),
          spec=as.double(spectrum),G11C=as.double(G_cos),
          G11=as.double(G_1),G12=as.double(G_2),T=as.integer(t),
          NFc=as.integer(n_cos_sin),ns=as.integer(n_cos))
  
  output = list(alpha_sim = bs$simAlpha)
  return(output)
}

####################

## ---

########
# FFBS #
########
## -> Forward Filtering Backward Sampling to sample from joint full conditional 
##    of the hidden state of a linear, Gaussian state space model
## -> Inputs:
##            -- 'y', observed data in t x N matrix
##            -- 'linpred', linear predictor in t x N matrix
##            -- 'Phi', observation matrix relation y to alpha
##            -- 'R', covariance of observation error nu
##            -- 'G', propagation matrix
##            -- 'Q', innovation covariance matrix of alpha
##            -- 'N', number of spatial points
##            -- 't', number of time points
##            -- 'K', dimension of latent process alpha
##            -- 'log_lik', whether to also return log-likelihood
##            -- 'backward_sample', whether to return sample of full 
##                conditional of latent process alpha
##            -- 'filtered', whether to return filtered values for alpha

## -> Output: list containing
##            -- 'alpha_simulation', t x N matrix containing sample from full 
##                conditional of latent process alpha
##            -- 'll', the log-likelihood
##            -- 'm_i_i', t x N matrix containing mean of filtered values of 
##                alpha

## The model is in the following form:
# --- y_t = linpred_t + Phi*alpha_t + nu_t; nu_t ~ N(0,R)
# --- alpha_t = G*alpha_{t-1} + epsilon_t; epsilon_t ~ N(0,Q)

FFBS <- function (y, linear_pred=matrix(0,ncol=dim(y)[2],nrow=dim(y)[1]), Phi, 
                  R, G, Q, N = dim(y)[2], t = dim(y)[1], K = dim(G)[1], 
                  log_lik = FALSE, backward_sample = TRUE, filtered = FALSE){
  
  ## Kalman Filter
  KF = KalmanFilter(y = y, linear_pred = linear_pred, Phi = Phi, R = R, G = G, 
                    Q = Q, N = N, t = t, K = K, log_lik = log_lik)
  
  m_i_i = KF$m_i_i
  m_i_i1 = KF$m_i_i1
  R_i_i = KF$R_i_i
  R_i_i1 = KF$R_i_i1
  ll = KF$ll
  alpha_simulation <- array(0, c(t+1, K))
  
  ## Backward Sample
  if (backward_sample) {
    R_i_i[t+1, , ] = (R_i_i[t+1, , ] + t(R_i_i[t+1, , ]))/2
    
    
    alpha_simulation[t+1, ] <- rmvnorm(1, mean = m_i_i[t+1, ],
                                       sigma = R_i_i[t+1, , ], 
                                       method = "chol")
    
    for (i in t:1) {
      L_i <- R_i_i[i, , ] %*% t(G) %*% solve(R_i_i1[i, , ])
      
      # Covariance
      R_i <- R_i_i[i, , ] - L_i %*% G %*% R_i_i[i, , ]
      R_i <- (R_i + t(R_i))/2
      
      # Mean
      m_i <- m_i_i[i, ] + L_i %*% (alpha_simulation[i+1, ] - m_i_i1[i,])
      
      # Simulate 
      alpha_simulation[i, ] <- rmvnorm(1, mean = m_i, sigma = R_i, 
                                       method = "chol")
    }
  }
  
  ## Initialise output
  output <- list()
  
  ## If 'backward_sample'=TRUE: return sample from full conditional of alpha
  if (backward_sample){
    
    ## Required block of 'alpha_simulation'
    alpha_simulation = alpha_simulation[-1, ]
    
    ## Format
    if(class(alpha_simulation)=="numeric"){
      alpha_simulation=t(matrix(alpha_simulation))
    }
    
    ## Amend output
    output <- c(output, list(alpha_simulation = alpha_simulation))
  } 
  
  ## If 'log_lik'=TRUE, return log-likelihood 
  if (log_lik) 
    output <- c(output, list(ll = ll))
  
  ## If 'filtered'=TRUE, return filtered values of alpha
  if (filtered) 
    output <- c(output, list(m_i_i = m_i_i[-1, ]))
  
  return(output)
}

########

#################
# SPECTRAL FFBS #
#################
## -> Forward Filtering Backward Sampling in the spectral space of the SPDE: 
##    i.e. to sample from joint full conditional of the coefficients alpha, and 
##    to evaluate the corresponding log-likelihood.
## -> Inputs:
##            -- 'w', observed data or latent process w in t x N matrix 
##                (default = NULL)
##            -- 'w_FFT_real', real Fourier transform of w in vector of length 
##                tN
##            -- 'spectrum', spectral density of the innovations in vector of 
##                length n^2. If not specified, generated from 'param' 
##            -- 'G_vec', propagation matrix in vector form. If not specified, 
##                generated from 'param'
##            -- 'tau2', variance of the measurement error. If not specified, 
##                set equal to param[9].
##            -- 'param', vector of parameters for the SPDE, in standard order.
##            -- 'n', number of spatial points on each axis (x,y)
##            -- 't', number of time points
##            -- 'K', number of Fourier functions (default = n^2)
##            -- 'log_lik', whether to also return log-likelihood (default = 
##                FALSE)
##            -- 'backward_sample', whether to return sample of full 
##                conditional of latent process alpha (default = TRUE)
##            -- 'cosine_indices', vector of indices of columns of cos terms 
##                in the 1:K real Fourier functions
##            -- 'n_cos', number of real Fourier functions containing only a 
##                cosine term (default = 4)
##            -- 'nu', smoothness parameter of Matern cov. fn. for innovations 
##                (default = 1 <-> Whittle cov. fn.)
##            -- 'dt', temporal lag between time points

## -> Output: list containing
##            -- 'alpha_simulation', t x N matrix containing sample from full 
##                conditional of latent process alpha
##            -- 'll', the log-likelihood

## The principal function 'FFBS_spectral' is a wrapper for the corefunction 
## 'ffbs_spectral', which in turn calls 'kf_spectral', 'bs_spectral', and 
## 'll_spectral'. 

FFBS_spectral <- function(w=NULL, w_FFT_real=NULL, innov_spectrum=NULL,
                          G_vec=NULL, tau2=NULL, param=NULL, n, t,
                          K=n^2, log_lik=FALSE, backward_sample=TRUE,
                          cosine_indices=(1:((n*n-4)/2)*2+3),
                          n_cos=4, nu=1, dt=1,m0 = rep(0,K), R0 = rep(0,K)){
  
  ## If 'w_FFT_real' not given, compute as inverse-FFT of w
  if(is.null(w_FFT_real)){
    w_FFT_real <- FFT_real_ST(STmatrix_to_vec(w),n=n,t=t,inv=TRUE)
  }
  
  ## If 'w_FFT_real' is given:
  else{
    # If 'w_FFT_real' in matrix form, convert into vector
    if(class(w_FFT_real)=="matrix"){
      w_FFT_real <- STmatrix_to_vec(w_FFT_real)
    }
    # If length of 'w_FFT_real' != t x n x n, output error message
    if(length(w_FFT_real)!=(n^2*t)){
      print(paste("Error: 'w_FFT_real' must be a vector of length t*n^2 or a ",
                  "matrix of dimension t x n^2.",sep=""))
      return()
    }
  }
  
  ## Select required terms from w_FFT_real
  if(K<n^2){
    basis_indices <- spde_initialise(n,t,K)$basis_indices
    w_FFT_indices <- rep(seq(0,n^2*(t-1),n^2),each=K)+basis_indices
    w_FFT_real = w_FFT_real[w_FFT_indices]
  }
  
  ## If required inputs not specified, output error message
  if((is.null(G_vec) & is.null(param)) | (is.null(innov_spectrum) & 
                                          is.null(param))){
    print("Either {'G_vec', 'spectrum', 'tau2'} or 'param' must be supplied.")
    return()
  }
  
  ## Compute propagator matrix G (vector form)
  if(is.null(G_vec)){
    SPDE_FT <- spde_initialise(n,t,K)
    wave <- SPDE_FT$wave
    cosine_indices <- SPDE_FT$cosine_indices
    n_cos <- SPDE_FT$n_cos
    G_vec <- propagator_vec(wave=wave,cosine_indices=cosine_indices,
                            zeta=param[3],rho1=param[4],gamma=param[5],
                            alpha=param[6],mu_x=param[7],mu_y=param[8],dt=dt,
                            n_cos=n_cos)
  }
  
  ## Compute innovation spectrum f(k_j) 
  if(is.null(innov_spectrum)){
    if(is.null(SPDE_FT)){
      SPDE_FT <- spde_initialise(n,t,K)
      wave <- SPDE_FT$wave
      cosine_indices <- SPDE_FT$cosine_indices
      n_cos <- SPDE_FT$n_cos
    }
    innov_spectrum <- innovation_spectrum(wave=wave,n=n,n_cos=n_cos,
                                          rho0=param[1],sigma2=param[2],
                                          zeta=param[3],rho1=param[4],
                                          gamma=param[5],alpha=param[6],
                                          nu=nu,dt=dt,norm=TRUE)
  }
  
  ## Compute tau^2 
  if(is.null(tau2)){
    tau2=param[9]
  }
  
  ## Run 'kf_spectral'
  KF <- kf_spectral(w_FFT_real = w_FFT_real,m_i_i1 = rep(0,K*t),
                    m_i_i=c(m0,rep(0,K*t)),R_i_i1=rep(0,K*t),
                    R_i_i=c(R0,rep(0,K*t)),G_cos=G_vec$G_cos,G_1=G_vec$G_1, 
                    G_2=G_vec$G_2,spectrum=innov_spectrum,t=t,n_cos=n_cos,
                    n_cos_sin=length(cosine_indices),tau2=tau2)
  
  ## Output
  output = list()
  
  ## Run 'bs_spectral'
  if(backward_sample==TRUE){
    alpha_sim = bs_spectral(alpha_sim = w_FFT_real, 
                            m_i_i = KF$m_i_i, m_i_i1 = KF$m_i_i1,
                            R_i_i = KF$R_i_i, R_i_i1 = KF$R_i_i1,
                            spectrum = innov_spectrum,G_cos = G_vec$G_cos, 
                            G_1 = G_vec$G_1, G_2 = G_vec$G_2,
                            t = t, n_cos_sin = length(cosine_indices),
                            n_cos = n_cos)$alpha_sim
    alpha_sim = vec_to_STmatrix(alpha_sim,t=t)
    output = c(output,list(alpha_sim=alpha_sim))
  }
  
  ## Run 'll_spectral'
  if(log_lik){
    ll = ll_spectral(log_lik=0,w_FFT_real, KF$m_i_i1, KF$R_i_i1, t, K, 
                     tau2)$ll
    
    if(K<n^2){
      ll = ll_non_spectral(log_lik=0,w,KF$m_i_i1,KF$R_i_i1,t,n,K,tau2)$ll
    }
    
    output = c(output,list(ll=ll))
  }
  
  return(output)
}

#################

## ---

##################
# LOG-LIKELIHOOD #
##################
## -> Evaluate log-likelihood of the parameters; this is a wrapper function for
##    FFBS_spectral, which itself is a wrapper for ll_non_spectral (K<n^2) or
##    ll_spectral (K=n^2)
## -> Inputs:
##            -- 'w', observed data or latent process w in t x N matrix 
##                (default = NULL)
##            -- 'w_FFT_real', real Fourier transform of w in vector of length 
##                tN
##            -- 'spectrum', spectral density of the innovations in vector of 
##                length n^2. If not specified, generated from 'param' 
##            -- 'G_vec', propagation matrix in vector form. If not specified, 
##                generated from 'param'
##            -- 'tau2', variance of the measurement error. If not specified, 
##                set equal to param[9].
##            -- 'param', vector of parameters for the SPDE, in  standard order.
##            -- 'n', number of spatial points on each axis (x,y)
##            -- 't', number of time points
##            -- 'K', number of Fourier functions (default = n^2)
##            -- 'cosine_indices', vector of indices of columns of cos terms 
##                in the 1:K real Fourier functions
##            -- 'n_cos', number of real Fourier functions containing only a 
##                cosine term (default = 4)
##            -- 'nu', smoothness parameter of Matern cov. fn. for innovations 
##                (default = 1 <-> Whittle cov. fn.)
##            -- 'dt', temporal lag between time points
##            -- 'x', covariates in array of dimension p x t x N
##            -- 'log_scale', indicates whether any parameters are provided on 
##                a logarithmic scale
##            -- 'log_indices', vector of integers provided on logarithmic 
##                scale
##            -- 'negative', whether to output the negative log-likelihood

## -> Output:
##            -- 'll', the log-likelihood

log_likelihood <- function(param_obs=NULL, param_sig = NULL, w=NULL, 
                           w_FFT_real=NULL, innov_spectrum=NULL, G_vec=NULL,
                           tau2=NULL, n, t, K=n^2, 
                           cosine_indices=(1:((n*n-4)/2)*2+3), n_cos=4,
                           nu=1, dt=1, x=NULL, log_scale=FALSE,
                           log_indices=c(1,2,3,4,5,9), negative=FALSE,
                           m0 = rep(0,K), R0 = rep(0,K)){
  
  ## Parameters ##
  param = c(param_sig,as.numeric(unlist(param_obs$tau2_vals)))
  
  ## If parameters provided on log-scale, rescale
  if(log_scale){
    param[log_indices] <- exp(param[log_indices])
  }
  
  ## Signal and observation parameters ##
  param_sig = param[1:8]
  param_obs = param_obs_func(m=1,m_indices=list(1:n^2),
                             tau2_vals=as.numeric(list(param[9])))
  
  ## Compute linear predictor X*beta
  if(!is.null(x)){
    linear_pred <- apply(x,c(2,3),linear_predictor,beta = param[-c(1:9)])
  }
  else{
    linear_pred <- 0
  }
  
  ## If 'w_FFT_real' not given, compute as inverse-FFT of (w-x*beta)
  if(is.null(w_FFT_real)){
    w_FFT_real <- FFT_real_ST(STmatrix_to_vec(w-linear_pred),n=n,t=t,inv=TRUE)
  }
  
  ## Compute log-likelihood using 'FFBS_spectral'
  ll <- KF_spectral(w=w,w_FFT_real=w_FFT_real,innov_spectrum=innov_spectrum,
                      G_vec=G_vec,tau2=tau2,param=param,n=n,t=t,log_lik=TRUE,
                      backward_sample=FALSE,K=K,cosine_indices=cosine_indices,
                      n_cos=n_cos,nu=nu,dt=dt,m0=m0,R0=R0)$ll
  
  ## Output
  if(negative){
    return(-ll)
  }
  else{
    return(ll)
  }
}

ll_non_spectral = function(log_lik=0,w,m_i_i1,R_i_i1,t,n,K,tau2){
  
  ## SPDE_FT object
  SPDE_FT = spde_initialise(n,t,K)
  wave = SPDE_FT$wave
  cosine_indices = SPDE_FT$cosine_indices
  n_cos = SPDE_FT$n_cos
  
  ## FFT matrix
  Phi = FFT_real_matrix(wave,cosine_indices,n_cos,n)
  Phi_T = t(Phi)
  
  ## R matrix
  R = diag(rep(tau2,n^2))
  
  ## Initialise
  ll = log_lik
  
  for (i in 1:t){
    S_i = R + Phi%*%diag(R_i_i1[(i-1)*K + 1:K])%*%Phi_T
    S_i_inv = solve(S_i)
    resid = w[i,] - Phi%*%m_i_i1[(i-1)*K + 1:K]
    
    ll = ll + determinant(S_i_inv, logarithm = TRUE)$modulus[[1]] - 
      t(resid) %*% S_i_inv %*% resid
  }
  
  # Normalise
  ll <- (ll - t*K*log(2*pi))/2
  
  return(list(ll=ll))
}

ll_spectral = function(log_lik=0,w_FFT_real,m_i_i1,R_i_i1,t,K,tau2){
  
  log_lik = .C("ll_spectral",ll=as.double(log_lik),wFT=as.double(w_FFT_real),
               mtt1=as.double(m_i_i1),rtt1=as.double(R_i_i1),
               T=as.integer(t),NF=as.integer(K),tau2=as.double(tau2))
  
  output = list(ll = log_lik$ll)
}

##################

## -- 

##########################
# SAMPLE FT COEFFICIENTS #
##########################
## -> Sample from the full conditional of the Fourier coefficients
##    (wrapper function for FFBS_spectral)
## -> Inputs:
##            -- 'w', observed data or latent process w in t x N matrix 
##                (default = NULL)
##            -- 'w_FFT_real', real Fourier transform of w in vector of length 
##                tN
##            -- 'innov_spectrum', spectral density of the innovations in 
##                vector of length n^2. If not specified, generated from 'param' 
##            -- 'G_vec', propagation matrix in vector form. If not specified, 
##                generated from 'param'
##            -- 'tau2', variance of the measurement error. If not specified, 
##                set equal to param[9].
##            -- 'param', vector of parameters for the SPDE, in standard order.
##            -- 'n', number of spatial points on each axis (x,y)
##            -- 't', number of time points
##            -- 'K', number of Fourier functions (default = n^2)
##            -- 'cosine_indices', vector of indices of columns of cos terms in
##                the 1:K real Fourier functions
##            -- 'n_cos', number of real Fourier functions containing only a 
##                cosine term (default = 4)
##            -- 'nu', smoothness parameter of Matern cov. fn. for innovations 
##                (default = 1 <-> Whittle cov. fn.)
##            -- 'dt', temporal lag between time points

## -> Output:
##            -- 'alpha', a t x N matrix containing a sample from the full 
##                conditional of latent process alpha

alpha_sample <- function(w = NULL,w_FFT_real = NULL,innov_spectrum = NULL, 
                         G_vec = NULL, tau2 = NULL, param = NULL, n, t, K = n^2,
                         cosine_indices = (1:((n*n-4)/2)*2+3), n_cos = 4, 
                         nu = 1, dt = 1){
  
  # Sample using FFBS_spectral
  alpha <- FFBS_spectral(w=w,w_FFT_real=w_FFT_real,
                         innov_spectrum=innov_spectrum,G_vec=G_vec,tau2=tau2,
                         param=param,n=n,t=t,K=K,log_lik=FALSE,
                         backward_sample=TRUE,cosine_indices=cosine_indices,
                         n_cos=n_cos,nu=nu,dt=dt)$alpha_sim
  
  return(alpha)
}

##########################

## --

##################
# LGSSM MATRICES #
##################

Phi_f <- function(wave,cosine_indices,n_cos=4,n,n_obs,coords){
  Phi_obs <- FFT_real_matrix(wave,cosine_indices,n_cos,n,n_obs,coords)
  return(Phi_obs)
}

XBeta_f <- function(param_bias,n,t,n_obs){
  
  ## Compute XBeta ##
  XBeta_obs = obs_bias_mat(n_obs,t,param_bias)
  
  ## Output ##
  return(XBeta_obs)
}

R_f <- function(param_obs,n,n_obs){
  
  ## Compute R matrix ##
  R_obs = obs_noise_mat(n_obs,param_obs)
  
  ## Output ##
  return(R_obs)
}

R_f_vec <- function(param_obs,n,n_obs){
  R = obs_noise_vec(n_obs,param_obs)
  return(R)
}

G_f <- function(param_sig,wave,cosine_indices,dt=1,n_cos=4){
  
  G = propagator_matrix(wave=wave,cosine_indices=cosine_indices,zeta=param_sig[3],
                        rho1=param_sig[4],gamma=param_sig[5],alpha=param_sig[6],
                        mu_x=param_sig[7],mu_y=param_sig[8],dt=dt,n_cos=n_cos)
  
  return(G)
}

G_f_vec <- function(param_sig,wave,cosine_indices,dt=1,n_cos=4){
  
  G_vec = propagator_vec(wave=wave,cosine_indices=cosine_indices,zeta=param_sig[3],
                         rho1=param_sig[4],gamma=param_sig[5],alpha=param_sig[6],
                         mu_x=param_sig[7],mu_y=param_sig[8],dt=dt,n_cos=n_cos)
  
  return(list(G_cos=G_vec$G_cos,G_1=G_vec$G_1,G_2=G_vec$G_2))
}

Q_f <- function(param_sig,wave,n,n_cos=4,nu=1,dt=1,norm=TRUE,
                spec_weights_mat = NULL){
  
  ## Compute Q matrix ##
  Q = diag(innovation_spectrum(wave=wave,n=n,n_cos=n_cos,rho0=param_sig[1],
                               sigma2=param_sig[2],zeta=param_sig[3],rho1=param_sig[4],
                               gamma=param_sig[5],alpha=param_sig[6],nu=nu,dt=dt,
                               norm = norm))
  
  ## Compute weighted Q matrix ##
  if(!is.null(spec_weights_mat)){
    Q = spec_weights_mat%*%Q%*%t(spec_weights_mat)
  }
  
  ## Output ##
  return(Q)
}

Q_f_vec <- function(param_sig,wave,n,n_cos=4,nu=1,dt=1,norm=TRUE){
  
  Q_vec = innovation_spectrum(wave=wave,n=n,n_cos=n_cos,rho0=param_sig[1],
                              sigma2=param_sig[2],zeta=param_sig[3],rho1=param_sig[4],
                              gamma=param_sig[5],alpha=param_sig[6],nu=nu,dt=dt,
                              norm=norm)
  
  return(Q_vec)
}

##################

######################
# LGSSM MATRIX GRADS #
######################

Phi_grad_f <- function(param_obs,param_sig,n,K,
                           param_bias = param_bias_func(),n_obs){
  
  ## Number of parameters
  n_param = param_obs$m + length(param_sig) + param_bias$m
  
  ## Compute Phi_grad
  Phi_obs_grad =  rep(list(matrix(0,nrow=n_obs,ncol=K)),n_param)
  
  return(Phi_obs_grad)
}

XBeta_grad_f <- function(param_obs, param_sig, n, t,
                             param_bias = param_bias_func(), n_obs){
  
  ## Number of parameters ##
  n_param = param_obs$m + length(param_sig) + param_bias$m
  
  ## Initalise ##
  XBeta_obs_grad = rep(list(matrix(0,nrow=t,ncol=n_obs)),n_param) 
  
  ## Compute XBeta_grad w.r.t observation bias parameters ##
  if(param_bias$m>0){
    for (j in 1:t){
      for (i in 1:param_bias$m){
        XBeta_obs_grad[[param_obs$m+length(param_sig)+i]][j,param_bias$m_indices[[i]]]=1
      }
    }
  }
  
  return(XBeta_obs_grad)
}

R_grad_f <- function(param_obs, param_sig, n,
                         param_bias = param_bias_func(), n_obs){
  
  ## Number of parameters ##
  n_param = param_obs$m + length(param_sig) + param_bias$m
  
  ## Initialise ##
  R_obs_grad = rep(list(matrix(0,nrow=n_obs,ncol=n_obs)),n_param)
  
  ## Compute R_grad w.r.t observation parameters ##
  for (i in 1:param_obs$m){
    diag(R_obs_grad[[length(param_sig)+i]])[param_obs$m_indices[[i]]]=1
  }
  
  ## Output ##
  return(R_obs_grad)
}

R_grad_f_vec <- function(param_obs,param_sig,n,
                             param_bias = param_bias_func(),n_obs){
  
  ## Number of parameters
  n_param = param_obs$m + length(param_sig) + param_bias$m
  
  ## Initialise
  R_grad_vec = rep(list(rep(0,n_obs)),n_param)
  
  ## Compute R_grad w.r.t observation parameters
  for (i in 1:(param_obs$m)){
    R_grad_vec[[length(param_sig)+i]][param_obs$m_indices[[i]]]=1
  }

  return(R_grad_vec)
}

G_grad_f <- function(param_obs,param_sig,n,K,wave,cosine_indices,dt=1,nu=1,
                         n_cos=4,param_bias = param_bias_func()){
  
  ## Number of parameters
  n_param = param_obs$m + length(param_sig) + param_bias$m
  
  ## Required parameters
  rho0 = param_sig[1]; sigma2 = param_sig[2]; zeta = param_sig[3]; rho1 = param_sig[4]
  gamma = param_sig[5]; alpha = param_sig[6]; mu_x = param_sig[7]; mu_y = param_sig[8]
  
  ## Initialise list
  G_grad = rep(list(matrix(0,nrow=K,ncol=K)),n_param)
  
  ## rho0, sigma2, tau2 ##
  G_grad[[1]] = G_grad[[2]] = G_grad[[9]] = matrix(0,nrow=K,ncol=K)
  
  
  ## Compute G ##
  G = propagator_matrix(wave=wave,cosine_indices=cosine_indices,zeta=zeta,
                        rho1=rho1,gamma=gamma,alpha=alpha,mu_x=mu_x,mu_y=mu_y,
                        dt=dt,n_cos=n_cos)
  ## zeta ##
  G_grad[[3]] = -dt*G
  
  ## rho1, gamma, psi, mu_x, mu_y ##
  
  ## Compute Sigma 
  if(rho1==0){
    Sig <- cbind(c(0,0),c(0,0))
  }
  else{
    Sig_Chol <- cbind(c(cos(alpha),-gamma*sin(alpha)),
                      c(sin(alpha),gamma*cos(alpha)))/rho1
    Sig <- solve(t(Sig_Chol)%*%Sig_Chol)
  }
  
  ## Compute Grad. Sigma
  Sigma_grad_rho1 <- 2 / rho1 * Sig
  
  Sigma_grad_gamma <- 
    2 * rho1^2 / gamma^3 * cbind(c(-sin(alpha)^2,sin(alpha)*cos(alpha)),
                                 c(sin(alpha)*cos(alpha),-cos(alpha)^2))
  
  Sigma_grad_alpha <- 
    rho1^2*(gamma^2-1) / gamma^2 * cbind(c(-2*sin(alpha)*cos(alpha),cos(alpha)^2-sin(alpha)^2),
                                         c(cos(alpha)^2-sin(alpha)^2,2*sin(alpha)*cos(alpha))) 
  
  ## Compute k_j^T * Grad Sigma * k_j
  rho1_multiplier <- apply(wave*Sigma_grad_rho1%*%wave,2,sum)
  gamma_multiplier <- apply(wave*Sigma_grad_gamma%*%wave,2,sum)
  alpha_multiplier <- apply(wave*Sigma_grad_alpha%*%wave,2,sum)
  
  ## Compute -Delta*(k_j^T * Sigma * k_j + zeta)
  ExpTerm <- -dt*apply(wave*Sig%*%wave,2,sum)-dt*rep(zeta,dim(wave)[2])
  
  ## Compute Delta * mu * k_j
  TrigTerm <- dt*c(mu_x,mu_y)%*%wave
  
  
  ## Remaining Gradients ##
  
  ## rho1 ##
  
  # Diagonal
  diag(G_grad[[4]]) <- 
    -dt*rho1_multiplier*c(exp(ExpTerm[1:n_cos]),exp(ExpTerm[(n_cos+1):K])*cos(TrigTerm[(n_cos+1):K]))
  
  # Off Diagonal
  diag(G_grad[[4]][cosine_indices,cosine_indices+1]) <-
    dt*rho1_multiplier[cosine_indices]*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices])
  
  diag(G_grad[[4]][cosine_indices+1,cosine_indices]) <- 
    -dt*rho1_multiplier[cosine_indices]*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices])
  
  
  ## gamma ##
  
  # Diagonal
  diag(G_grad[[5]]) <- 
    -dt*gamma_multiplier*c(exp(ExpTerm[1:n_cos]),
                           exp(ExpTerm[(n_cos+1):K])*cos(TrigTerm[(n_cos+1):K]))
  
  # Off Diagonal
  diag(G_grad[[5]][cosine_indices,cosine_indices+1]) <-
    dt*(gamma_multiplier[cosine_indices]*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices]))
  
  diag(G_grad[[5]][cosine_indices+1,cosine_indices]) <- 
    -dt*(gamma_multiplier[cosine_indices]*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices]))
  
  
  ## alpha ##
  
  # Diagonal
  diag(G_grad[[6]]) <- 
    -dt*alpha_multiplier*c(exp(ExpTerm[1:n_cos]),
                           exp(ExpTerm[(n_cos+1):K])*cos(TrigTerm[(n_cos+1):K]))
  
  # Off Diagonal
  diag(G_grad[[6]][cosine_indices,cosine_indices+1]) <-
    dt*alpha_multiplier[cosine_indices]*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices])
  
  diag(G_grad[[6]][cosine_indices+1,cosine_indices]) <- 
    -dt*alpha_multiplier[cosine_indices]*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices])
  
  
  ## mu_x ##
  
  # Diagonal
  diag(G_grad[[7]]) <- 
    -dt*wave[1,]*c(rep(0,n_cos),
                   exp(ExpTerm[(n_cos+1):K])*sin(TrigTerm[(n_cos+1):K]))
  
  # Off Diagonal
  diag(G_grad[[7]][cosine_indices,cosine_indices+1]) <-
    -dt*wave[1,cosine_indices]*exp(ExpTerm[cosine_indices])*cos(TrigTerm[cosine_indices])
  
  diag(G_grad[[7]][cosine_indices+1,cosine_indices]) <- 
    dt*wave[1,cosine_indices]*exp(ExpTerm[cosine_indices])*cos(TrigTerm[cosine_indices])
  
  
  ## mu_y ##
  
  # Diagonal
  diag(G_grad[[8]]) <- 
    -dt*wave[2,]*c(rep(0,n_cos),
                   exp(ExpTerm[(n_cos+1):K])*sin(TrigTerm[(n_cos+1):K]))
  
  # Off Diagonal
  diag(G_grad[[8]][cosine_indices,cosine_indices+1]) <- 
    -dt*wave[2,cosine_indices]*exp(ExpTerm[cosine_indices])*cos(TrigTerm[cosine_indices])
  
  diag(G_grad[[8]][cosine_indices+1,cosine_indices]) <- 
    dt*wave[2,cosine_indices]*exp(ExpTerm[cosine_indices])*cos(TrigTerm[cosine_indices])
  
  
  ## Output
  return(G_grad)
}

G_grad_f_vec <- function(param_obs,param_sig,n,K,wave,cosine_indices,dt=1,
                             nu=1,n_cos=4,param_bias = param_bias_func()){
  
  ## Number of parameters
  n_param = param_obs$m + length(param_sig) + param_bias$m
  
  ## Compute auxiliary values
  n_cos_sin = length(cosine_indices)
  
  ## Required parameters
  rho0 = param_sig[1]; sigma2 = param_sig[2]; zeta = param_sig[3]; rho1 = param_sig[4]
  gamma = param_sig[5]; alpha = param_sig[6]; mu_x = param_sig[7]; mu_y = param_sig[8]
  
  ## Initialise vectors
  G_grad_aux = list(G_grad_cos = rep(0,n_cos),G_grad_1 = rep(0,n_cos_sin),
                    G_grad_2 = rep(0,n_cos_sin))
  G_grad = rep(list(G_grad_aux),n_param)
  
  ## rho0, sigma2, tau2: leave as zero ##
  
  ## Compute G ##
  G = propagator_vec(wave=wave,cosine_indices=cosine_indices,zeta=zeta,
                     rho1=rho1,gamma=gamma,alpha=alpha,mu_x=mu_x,mu_y=mu_y,
                     dt=dt,n_cos=n_cos)
  
  ## zeta ##
  G_grad[[3]]$G_grad_cos = -dt*G$G_cos
  G_grad[[3]]$G_grad_1 = -dt*G$G_1
  G_grad[[3]]$G_grad_2 = -dt*G$G_2
  
  ## rho1, gamma, psi, mu_x, mu_y ##
  
  # Sigma 
  if(rho1==0){
    Sig <- cbind(c(0,0),c(0,0))
  }
  else{
    Sig_Chol <- cbind(c(cos(alpha),-gamma*sin(alpha)),
                      c(sin(alpha),gamma*cos(alpha)))/rho1
    Sig <- solve(t(Sig_Chol)%*%Sig_Chol)
  }
  
  # Grad. Sigma
  Sigma_grad_rho1 <- 2 / rho1 * Sig
  
  Sigma_grad_gamma <- 
    2 * rho1^2 / gamma^3 * cbind(c(-sin(alpha)^2,sin(alpha)*cos(alpha)),
                                 c(sin(alpha)*cos(alpha),-cos(alpha)^2))
  
  Sigma_grad_alpha <- 
    rho1^2*(gamma^2-1) / gamma^2 * cbind(c(-2*sin(alpha)*cos(alpha),cos(alpha)^2-sin(alpha)^2),
                                         c(cos(alpha)^2-sin(alpha)^2,2*sin(alpha)*cos(alpha))) 
  
  # Compute k_j^T * Grad.Sigma * k_j
  rho1_multiplier <- apply(wave*Sigma_grad_rho1%*%wave,2,sum)
  gamma_multiplier <- apply(wave*Sigma_grad_gamma%*%wave,2,sum)
  alpha_multiplier <- apply(wave*Sigma_grad_alpha%*%wave,2,sum)
  
  # Compute -Delta*(k_j^T * Sigma * k_j + zeta)
  ExpTerm <- -dt*apply(wave*Sig%*%wave,2,sum)-dt*rep(zeta,dim(wave)[2])
  
  # Compute Delta * mu * k_j
  TrigTerm <- dt*c(mu_x,mu_y)%*%wave
  
  ## Compute remaining gradients ##
  
  ## rho1 ##
  G_grad[[4]]$G_grad_cos <- -dt*rho1_multiplier[1:n_cos]*exp(ExpTerm[1:n_cos])
  G_grad[[4]]$G_grad_1 <- -dt*rho1_multiplier[cosine_indices]*exp(ExpTerm[cosine_indices])*cos(TrigTerm[cosine_indices])
  G_grad[[4]]$G_grad_2 <- dt*rho1_multiplier[cosine_indices]*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices])
  
  ## gamma ##
  G_grad[[5]]$G_grad_cos <- -dt*gamma_multiplier[1:n_cos]*exp(ExpTerm[1:n_cos])
  G_grad[[5]]$G_grad_1 <- -dt*gamma_multiplier[cosine_indices]*exp(ExpTerm[cosine_indices])*cos(TrigTerm[cosine_indices])
  G_grad[[5]]$G_grad_2 <- dt*gamma_multiplier[cosine_indices]*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices])
  
  ## alpha ##
  G_grad[[6]]$G_grad_cos <- -dt*alpha_multiplier[1:n_cos]*exp(ExpTerm[1:n_cos])
  G_grad[[6]]$G_grad_1 <- -dt*alpha_multiplier[cosine_indices]*exp(ExpTerm[cosine_indices])*cos(TrigTerm[cosine_indices])
  G_grad[[6]]$G_grad_2 <- dt*alpha_multiplier[cosine_indices]*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices])
  
  ## mu_x ##
  G_grad[[7]]$G_grad_1 <- -dt*wave[1,cosine_indices]*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices])
  G_grad[[7]]$G_grad_2 <- -dt*wave[1,cosine_indices]*exp(ExpTerm[cosine_indices])*cos(TrigTerm[cosine_indices])
  
  ## mu_y ##
  G_grad[[8]]$G_grad_1 <- -dt*wave[2,cosine_indices]*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices])
  G_grad[[8]]$G_grad_2 <- -dt*wave[2,cosine_indices]*exp(ExpTerm[cosine_indices])*cos(TrigTerm[cosine_indices])
  
  return(G_grad)
}

Q_grad_f <- function(param_obs,param_sig,n,K,wave,n_cos=4,nu=1,dt=1,norm=TRUE,
                     param_bias = param_bias_func(),spec_weights_mat = NULL){
  
  ## Number of parameters
  n_param = param_obs$m + length(param_sig) + param_bias$m
  
  ## Required parameters
  rho0 = param_sig[1]; sigma2 = param_sig[2]; zeta = param_sig[3]; rho1 = param_sig[4]
  gamma = param_sig[5]; alpha = param_sig[6]; mu_x = param_sig[7]; mu_y = param_sig[8]
  
  ## Initialise list
  Q_grad = rep(list(matrix(0,nrow=K,ncol=K)),n_param)
  
  ## mu_x,mu_y,tau2 ##
  Q_grad[[7]] = Q_grad[[8]] = Q_grad[[9]] = diag(0,nrow=K,ncol=K)
  
  ## rho0,sigma2,zeta,rho1,gamma,alpha ##
  
  ## Compute spectral density f(k_j) ##
  
  ## Compute k_j^T k_j
  w2 <- apply(wave^2,2,sum)
  
  ## Compute f(k_j)
  numerator <- 1
  denominator <- ((1/rho0)^2 + w2)^(nu + 1)
  spec <- numerator/denominator
  
  ## Adjustment for cosine only terms
  spec[1:n_cos] <- spec[1:n_cos]/2
  
  ## Normalise
  if(norm){
    spec_normalised <- spec*(n^2)/sum(spec)
  }
  else{
    spec_normalised <- spec/sum(spec)
  }
  
  ## Rescale by sigma^2
  spec_normalised <- sigma2 * spec_normalised
  
  ## Compute df(k_j)/d(theta)
  spec_normalised_grad_sigma2 <- 1/sigma2*spec_normalised
  spec_normalised_grad_rho0 <- 
    4/rho0^3*((w2+(1/rho0^2))^(-1)-
                sum((w2+(1/rho0^2))^(-1)*spec)/sum(spec))*spec_normalised
  
  ## Compute Sigma
  if(rho1==0){
    Sig <- cbind(c(0,0),c(0,0))
  }
  else{
    Sig_Chol <- cbind(c(cos(alpha),-gamma*sin(alpha)),
                      c(sin(alpha),gamma*cos(alpha)))/rho1
    Sig <- solve(t(Sig_Chol)%*%Sig_Chol)
  }
  
  ## Compute Grad. Sigma
  Sigma_grad_rho1 <- 2 / rho1 * Sig
  
  Sigma_grad_gamma <- 
    2 * rho1^2 / gamma^3 * cbind(c(-sin(alpha)^2,sin(alpha)*cos(alpha)),
                                 c(sin(alpha)*cos(alpha),-cos(alpha)^2))
  
  Sigma_grad_alpha <- 
    rho1^2*(gamma^2-1) / gamma^2 * cbind(c(-2*sin(alpha)*cos(alpha),cos(alpha)^2-sin(alpha)^2),
                                         c(cos(alpha)^2-sin(alpha)^2,2*sin(alpha)*cos(alpha))) 
  
  ## Compute k_j^T * Grad. Sigma * k_j
  rho1_multiplier <- apply(wave*Sigma_grad_rho1%*%wave,2,sum)
  gamma_multiplier <- apply(wave*Sigma_grad_gamma%*%wave,2,sum)
  alpha_multiplier <- apply(wave*Sigma_grad_alpha%*%wave,2,sum)
  
  ## Compute -(k_j^T * Sigma * k_j + zeta)
  ExpTerm <- -apply(wave*Sig%*%wave,2,sum)-rep(zeta,dim(wave)[2])
  
  
  ## Remaining Gradients ##
  
  ## rho0 ##
  Q_grad[[1]] = diag(spec_normalised_grad_rho0*(1-exp(2*dt*ExpTerm))/(-2*ExpTerm))
  
  ## sigma2 ##
  Q_grad[[2]] = diag(spec_normalised_grad_sigma2*(1-exp(2*dt*ExpTerm))/(-2*ExpTerm))
  
  ## zeta ##
  numerator_zeta = (-2*dt*ExpTerm+1)*exp(2*dt*ExpTerm)-rep(1,dim(wave)[2])
  denominator_zeta = 2*(-ExpTerm)^2
  Q_grad[[3]] = diag(spec_normalised*numerator_zeta/denominator_zeta)
  
  ## rho1 ##
  numerator_rho1 = ((-2*dt*ExpTerm+1)*exp(2*dt*ExpTerm)-
                      rep(1,dim(wave)[2]))*rho1_multiplier
  denominator_rho1 = 2*(-ExpTerm)^2
  Q_grad[[4]] = diag(spec_normalised*numerator_rho1/denominator_rho1)
  
  ## gamma ##
  numerator_gamma = ((-2*dt*ExpTerm+1)*exp(2*dt*ExpTerm)-
                       rep(1,dim(wave)[2]))*gamma_multiplier
  denominator_gamma = 2*(-ExpTerm)^2
  Q_grad[[5]] = diag(spec_normalised*numerator_gamma/denominator_gamma)
  
  ## alpha ##
  numerator_alpha = ((-2*dt*ExpTerm+1)*exp(2*dt*ExpTerm)-
                       rep(1,dim(wave)[2]))*alpha_multiplier
  denominator_alpha = 2*(-ExpTerm)^2
  Q_grad[[6]] = diag(spec_normalised*numerator_alpha/denominator_alpha)
  
  
  ## Compute weighted matrix ##
  if(!is.null(spec_weights_mat)){
    Q_grad = lapply(1:n_param, function(i) spec_weights_mat%*%Q_grad[[i]]%*%t(spec_weights_mat))
  }
  
  ## Output ##
  return(Q_grad)
}

Q_grad_f_vec <- function(param_obs,param_sig,n,K,wave,n_cos=4,nu=1,dt=1,
                             norm=TRUE,param_bias = param_bias_func()){
  
  ## Number of parameters
  n_param = param_obs$m + length(param_sig) + param_bias$m
  
  ## Required parameters
  rho0 = param_sig[1]; sigma2 = param_sig[2]; zeta = param_sig[3]; rho1 = param_sig[4]
  gamma = param_sig[5]; alpha = param_sig[6]; mu_x = param_sig[7]; mu_y = param_sig[8]
  
  ## Initialise array
  Q_grad_vec = rep(list(rep(0,K)),n_param)
  
  ## mu_x,mu_y,tau2 : leave as zero ##
  
  ## rho0,sigma2,zeta,rho1,gamma,alpha ##
  
  ## Compute spectral density f(k_j) ##
  
  ## Compute k_j^T k_j
  w2 <- apply(wave^2,2,sum)
  
  ## Compute f(k_j)
  numerator <- 1
  denominator <- ((1/rho0)^2 + w2)^(nu + 1)
  spec <- numerator/denominator
  
  ## Adjustment for cosine only terms
  spec[1:n_cos] <- spec[1:n_cos]/2
  
  ## Normalise
  if(norm){
    spec_normalised <- spec*(n^2)/sum(spec)
  }
  else{
    spec_normalised <- spec/sum(spec)
  }
  
  ## Rescale by sigma^2
  spec_normalised <- sigma2 * spec_normalised
  
  ## Compute df(k_j)/d(theta)
  spec_normalised_grad_sigma2 <- 1/sigma2*spec_normalised
  spec_normalised_grad_rho0 <- 
    4/rho0^3*((w2+(1/rho0^2))^(-1)-
                sum((w2+(1/rho0^2))^(-1)*spec)/sum(spec))*spec_normalised
  
  ## Compute Sigma
  if(rho1==0){
    Sig <- cbind(c(0,0),c(0,0))
  }
  else{
    Sig_Chol <- cbind(c(cos(alpha),-gamma*sin(alpha)),
                      c(sin(alpha),gamma*cos(alpha)))/rho1
    Sig <- solve(t(Sig_Chol)%*%Sig_Chol)
  }
  
  ## Compute Grad. Sigma
  Sigma_grad_rho1 <- 2 / rho1 * Sig
  
  Sigma_grad_gamma <- 
    2 * rho1^2 / gamma^3 * cbind(c(-sin(alpha)^2,sin(alpha)*cos(alpha)),
                                 c(sin(alpha)*cos(alpha),-cos(alpha)^2))
  
  Sigma_grad_alpha <- 
    rho1^2*(gamma^2-1) / gamma^2 * cbind(c(-2*sin(alpha)*cos(alpha),cos(alpha)^2-sin(alpha)^2),
                                         c(cos(alpha)^2-sin(alpha)^2,2*sin(alpha)*cos(alpha))) 
  
  ## Compute k_j^T * Grad. Sigma * k_j
  rho1_multiplier <- apply(wave*Sigma_grad_rho1%*%wave,2,sum)
  gamma_multiplier <- apply(wave*Sigma_grad_gamma%*%wave,2,sum)
  alpha_multiplier <- apply(wave*Sigma_grad_alpha%*%wave,2,sum)
  
  ## Compute -(k_j^T * Sigma * k_j + zeta)
  ExpTerm <- -apply(wave*Sig%*%wave,2,sum)-rep(zeta,dim(wave)[2])
  
  
  ## Compute Remaining Gradients ##
  
  ## rho0 ##
  Q_grad_vec[[1]] = spec_normalised_grad_rho0*(1-exp(2*dt*ExpTerm))/(-2*ExpTerm)
  
  ## sigma2 ##
  Q_grad_vec[[2]] = spec_normalised_grad_sigma2*(1-exp(2*dt*ExpTerm))/(-2*ExpTerm)
  
  ## zeta ##
  numerator_zeta = (-2*dt*ExpTerm+1)*exp(2*dt*ExpTerm)-rep(1,dim(wave)[2])
  denominator_zeta = 2*(-ExpTerm)^2
  Q_grad_vec[[3]] = spec_normalised*numerator_zeta/denominator_zeta
  
  ## rho1 ##
  numerator_rho1 = ((-2*dt*ExpTerm+1)*exp(2*dt*ExpTerm)-
                      rep(1,dim(wave)[2]))*rho1_multiplier
  denominator_rho1 = 2*(-ExpTerm)^2
  Q_grad_vec[[4]] = spec_normalised*numerator_rho1/denominator_rho1
  
  ## gamma ##
  numerator_gamma = ((-2*dt*ExpTerm+1)*exp(2*dt*ExpTerm)-
                       rep(1,dim(wave)[2]))*gamma_multiplier
  denominator_gamma = 2*(-ExpTerm)^2
  Q_grad_vec[[5]] = spec_normalised*numerator_gamma/denominator_gamma
  
  ## alpha ##
  numerator_alpha = ((-2*dt*ExpTerm+1)*exp(2*dt*ExpTerm)-
                       rep(1,dim(wave)[2]))*alpha_multiplier
  denominator_alpha = 2*(-ExpTerm)^2
  Q_grad_vec[[6]] = spec_normalised*numerator_alpha/denominator_alpha
  
  ## Output ##
  return(Q_grad_vec)
}

######################

##############################
# LGSSM MATRIX SPATIAL GRADS #
##############################

Phi_grad_space_f <- function(wave,cosine_indices,n_cos,n,n_obs,coords,
                                     grad_coord_index){
  Phi_obs_grad_space <- FFT_real_matrix_grad(wave,cosine_indices,n_cos,n,
                                             n_obs,coords,grad_coord_index)
  return(Phi_obs_grad_space)
}

XBeta_grad_space_f <- function(n_obs,t){
  XBeta_obs_grad_space <- rep(list(matrix(0,nrow=t,ncol=n_obs)),2*n_obs)
  return(XBeta_obs_grad_space)
}

R_grad_space_f <- function(n_obs){
  R_obs_grad_space = rep(list(matrix(0,ncol=n_obs,nrow=n_obs)),2*n_obs)
  return(R_obs_grad_space)
}

G_grad_space_f <- function(n_obs,K){
  G_grad_space <- rep(list(matrix(0,nrow=K,ncol=K)),2*n_obs)
  return(G_grad_space)
}

Q_grad_space_f <- function(n_obs,K){
  Q_grad_space <- rep(list(matrix(0,nrow=K,ncol=K)),2*n_obs)
  return(Q_grad_space)
}

##############################

#######################
# LGSSM MATRIX GRAD2s #
#######################

Phi_grad2_f <- function(param_obs,param_sig,n,K,
                            param_bias = param_bias_func(),n_obs){
  
  ## Number of parameters
  n_param = param_obs$m + length(param_sig) + param_bias$m
  
  ## Compute Phi_grad
  Phi_obs_grad =  rep(list(matrix(0,nrow=n_obs,ncol=K)),n_param)
  
  return(Phi_obs_grad)
}

XBeta_grad2_f <- function(param_obs,param_sig,n,t,
                              param_bias=param_bias_func(),n_obs){
  
  ## Number of parameters ##
  n_param = param_obs$m + length(param_sig) + param_bias$m
  
  ## Compute ##
  XBeta_obs_grad = rep(list(matrix(0,nrow=t,ncol=n_obs)),n_param) 
  
  ## Output ##
  return(XBeta_obs_grad)
}

R_grad2_f <- function(param_obs,param_sig,n,
                          param_bias = param_bias_func(),n_obs){
  
  ## Number of parameters ##
  n_param = param_obs$m + length(param_sig) + param_bias$m
  
  ## Compute ##
  R_obs_grad = rep(list(matrix(0,nrow=n_obs,ncol=n_obs)),n_param)
  
  ## Output ##
  return(R_obs_grad)
}

G_grad2_f <- function(param_obs,param_sig,n,K,wave,cosine_indices,dt=1,
                          nu=1,n_cos=4,param_bias = param_bias_func()){
  
  ## Number of parameters
  n_param = param_obs$m + length(param_sig) + param_bias$m
  
  ## Required parameters
  rho0 = param_sig[1]; sigma2 = param_sig[2]; zeta = param_sig[3]; rho1 = param_sig[4]
  gamma = param_sig[5]; alpha = param_sig[6]; mu_x = param_sig[7]; mu_y = param_sig[8]
  
  ## Initialise list
  G_grad2 = rep(list(matrix(0,nrow=K,ncol=K)),n_param)
  
  ## rho0, sigma2, tau2 ##
  G_grad2[[1]] = G_grad2[[2]] = G_grad2[[9]] = matrix(0,nrow=K,ncol=K)
  
  
  ## Compute G ##
  G = propagator_matrix(wave=wave,cosine_indices=cosine_indices,zeta=zeta,
                        rho1=rho1,gamma=gamma,alpha=alpha,mu_x=mu_x,mu_y=mu_y,
                        dt=dt,n_cos=n_cos)
  ## zeta ##
  G_grad2[[3]] = dt^2*G
  
  ## rho1, gamma, psi, mu_x, mu_y ##
  
  ## Compute Sigma 
  if(rho1==0){
    Sig <- cbind(c(0,0),c(0,0))
  }
  else{
    Sig_Chol <- cbind(c(cos(alpha),-gamma*sin(alpha)),
                      c(sin(alpha),gamma*cos(alpha)))/rho1
    Sig <- solve(t(Sig_Chol)%*%Sig_Chol)
  }
  
  ## Compute Grad. Sigma
  Sigma_grad_rho1 <- 2 / rho1 * Sig
  
  Sigma_grad_gamma <- 
    2 * rho1^2 / gamma^3 * cbind(c(-sin(alpha)^2,sin(alpha)*cos(alpha)),
                                 c(sin(alpha)*cos(alpha),-cos(alpha)^2))
  
  Sigma_grad_alpha <- 
    rho1^2*(gamma^2-1) / gamma^2 * cbind(c(-2*sin(alpha)*cos(alpha),cos(alpha)^2-sin(alpha)^2),
                                         c(cos(alpha)^2-sin(alpha)^2,2*sin(alpha)*cos(alpha))) 
  
  ## Compute k_j^T * Grad Sigma * k_j
  rho1_multiplier <- apply(wave*Sigma_grad_rho1%*%wave,2,sum)
  gamma_multiplier <- apply(wave*Sigma_grad_gamma%*%wave,2,sum)
  alpha_multiplier <- apply(wave*Sigma_grad_alpha%*%wave,2,sum)
  
  
  ## Compute Grad^2 Sigma
  Sigma_grad2_rho1 <- 2 / rho1^2 * Sig
  
  Sigma_grad2_gamma <- 
    -6 * rho1^2 / gamma^4 * cbind(c(-sin(alpha)^2,sin(alpha)*cos(alpha)),
                                  c(sin(alpha)*cos(alpha),-cos(alpha)^2))
  
  Sigma_grad2_alpha <- 
    2*rho1^2*(gamma^2-1) / gamma^2 * cbind(c(sin(alpha)^2-cos(alpha)^2,-2*sin(alpha)*cos(alpha)),
                                           c(-2*sin(alpha)*cos(alpha),cos(alpha)^2-sin(alpha)^2))
  
  ## Compute k_j^T * Grad^2 Sigma * k_j
  rho1_multiplier_2 <- apply(wave*Sigma_grad2_rho1%*%wave,2,sum)
  gamma_multiplier_2 <- apply(wave*Sigma_grad2_gamma%*%wave,2,sum)
  alpha_multiplier_2 <- apply(wave*Sigma_grad2_alpha%*%wave,2,sum)
  
  
  ## Compute -Delta*(k_j^T * Sigma * k_j + zeta)
  ExpTerm <- -dt*apply(wave*Sig%*%wave,2,sum)-dt*rep(zeta,dim(wave)[2])
  
  ## Compute Delta * mu * k_j
  TrigTerm <- dt*c(mu_x,mu_y)%*%wave
  
  
  ## Remaining Gradients ##
  
  ## rho1 ##
  
  # Diagonal
  diag(G_grad2[[4]]) <- 
    (-dt*rho1_multiplier_2+dt^2*(rho1_multiplier)^2)*c(exp(ExpTerm[1:n_cos]),exp(ExpTerm[(n_cos+1):K])*cos(TrigTerm[(n_cos+1):K]))
  
  # Off Diagonal
  diag(G_grad2[[4]][cosine_indices,cosine_indices+1]) <-
    -(-dt*rho1_multiplier_2[cosine_indices]+dt^2*(rho1_multiplier[cosine_indices])^2)*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices])
  
  diag(G_grad2[[4]][cosine_indices+1,cosine_indices]) <- 
    (-dt*rho1_multiplier_2[cosine_indices]+dt^2*(rho1_multiplier[cosine_indices])^2)*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices])
  
  
  ## gamma ##
  
  # Diagonal
  diag(G_grad2[[5]]) <- 
    (-dt*gamma_multiplier_2+dt^2*(gamma_multiplier)^2)*c(exp(ExpTerm[1:n_cos]),exp(ExpTerm[(n_cos+1):K])*cos(TrigTerm[(n_cos+1):K]))
  
  # Off Diagonal
  diag(G_grad2[[5]][cosine_indices,cosine_indices+1]) <-
    -(-dt*gamma_multiplier_2[cosine_indices]+dt^2*(gamma_multiplier[cosine_indices])^2)*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices])
  
  diag(G_grad2[[5]][cosine_indices+1,cosine_indices]) <- 
    (-dt*gamma_multiplier_2[cosine_indices]+dt^2*(gamma_multiplier[cosine_indices])^2)*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices])
  
  
  ## alpha ##
  
  # Diagonal
  diag(G_grad2[[6]]) <- 
    (-dt*alpha_multiplier_2+dt^2*(alpha_multiplier)^2)*c(exp(ExpTerm[1:n_cos]),exp(ExpTerm[(n_cos+1):K])*cos(TrigTerm[(n_cos+1):K]))
  
  # Off Diagonal
  diag(G_grad2[[6]][cosine_indices,cosine_indices+1]) <-
    -(-dt*alpha_multiplier_2[cosine_indices]+dt^2*(alpha_multiplier[cosine_indices])^2)*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices])
  
  diag(G_grad2[[6]][cosine_indices+1,cosine_indices]) <- 
    (-dt*alpha_multiplier_2[cosine_indices]+dt^2*(alpha_multiplier[cosine_indices])^2)*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices])
  
  
  ## mu_x ##
  
  # Diagonal
  diag(G_grad2[[7]]) <- 
    -dt^2*wave[1,]^2*c(rep(0,n_cos),exp(ExpTerm[(n_cos+1):K])*cos(TrigTerm[(n_cos+1):K]))
  
  # Off Diagonal
  diag(G_grad2[[7]][cosine_indices,cosine_indices+1]) <-
    dt^2*wave[1,cosine_indices]^2*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices])
  
  diag(G_grad2[[7]][cosine_indices+1,cosine_indices]) <- 
    -dt^2*wave[1,cosine_indices]^2*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices])
  
  
  ## mu_y ##
  
  # Diagonal
  diag(G_grad2[[8]]) <- 
    -dt^2*wave[2,]^2*c(rep(0,n_cos),exp(ExpTerm[(n_cos+1):K])*cos(TrigTerm[(n_cos+1):K]))
  
  # Off Diagonal
  diag(G_grad2[[8]][cosine_indices,cosine_indices+1]) <- 
    dt^2*wave[2,cosine_indices]^2*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices])
  
  diag(G_grad2[[8]][cosine_indices+1,cosine_indices]) <- 
    -dt^2*wave[2,cosine_indices]^2*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices])
  
  
  ## Output
  return(G_grad2)
}

G_grad2_f_vec <- function(param_obs,param_sig,n,K,wave,cosine_indices,dt,
                              nu=1,n_cos=4,param_bias = param_bias_func()){
  
  ## Number of parameters
  n_param = param_obs$m + length(param_sig) + param_bias$m
  
  ## Compute auxiliary values
  n_cos_sin = length(cosine_indices)
  
  ## Required parameters
  rho0 = param_sig[1]; sigma2 = param_sig[2]; zeta = param_sig[3]; rho1 = param_sig[4]
  gamma = param_sig[5]; alpha = param_sig[6]; mu_x = param_sig[7]; mu_y = param_sig[8]
  
  ## Initialise vectors
  G_grad2_aux = list(G_grad2_cos = rep(0,n_cos),G_grad2_1 = rep(0,n_cos_sin),
                     G_grad2_2 = rep(0,n_cos_sin))
  G_grad2 = rep(list(G_grad2_aux),n_param)
  
  ## rho0, sigma2, tau2: leave as zero ##
  
  ## Compute G ##
  G = propagator_vec(wave=wave,cosine_indices=cosine_indices,zeta=zeta,
                     rho1=rho1,gamma=gamma,alpha=alpha,mu_x=mu_x,mu_y=mu_y,
                     dt=dt,n_cos=n_cos)
  
  ## zeta ##
  G_grad2[[3]]$G_grad2_cos = dt^2*G$G_cos
  G_grad2[[3]]$G_grad2_1 = dt^2*G$G_1
  G_grad2[[3]]$G_grad2_2 = dt^2*G$G_2
  
  ## rho1, gamma, psi, mu_x, mu_y ##
  
  ## Compute Sigma 
  if(rho1==0){
    Sig <- cbind(c(0,0),c(0,0))
  }
  else{
    Sig_Chol <- cbind(c(cos(alpha),-gamma*sin(alpha)),
                      c(sin(alpha),gamma*cos(alpha)))/rho1
    Sig <- solve(t(Sig_Chol)%*%Sig_Chol)
  }
  
  ## Compute Grad. Sigma
  Sigma_grad_rho1 <- 2 / rho1 * Sig
  
  Sigma_grad_gamma <- 
    2 * rho1^2 / gamma^3 * cbind(c(-sin(alpha)^2,sin(alpha)*cos(alpha)),
                                 c(sin(alpha)*cos(alpha),-cos(alpha)^2))
  
  Sigma_grad_alpha <- 
    rho1^2*(gamma^2-1) / gamma^2 * cbind(c(-2*sin(alpha)*cos(alpha),cos(alpha)^2-sin(alpha)^2),
                                         c(cos(alpha)^2-sin(alpha)^2,2*sin(alpha)*cos(alpha))) 
  
  ## Compute k_j^T * Grad Sigma * k_j
  rho1_multiplier <- apply(wave*Sigma_grad_rho1%*%wave,2,sum)
  gamma_multiplier <- apply(wave*Sigma_grad_gamma%*%wave,2,sum)
  alpha_multiplier <- apply(wave*Sigma_grad_alpha%*%wave,2,sum)
  
  
  ## Compute Grad^2 Sigma
  Sigma_grad2_rho1 <- 2 / rho1^2 * Sig
  
  Sigma_grad2_gamma <- 
    -6 * rho1^2 / gamma^4 * cbind(c(-sin(alpha)^2,sin(alpha)*cos(alpha)),
                                  c(sin(alpha)*cos(alpha),-cos(alpha)^2))
  
  Sigma_grad2_alpha <- 
    2*rho1^2*(gamma^2-1) / gamma^2 * cbind(c(sin(alpha)^2-cos(alpha)^2,-2*sin(alpha)*cos(alpha)),
                                           c(-2*sin(alpha)*cos(alpha),cos(alpha)^2-sin(alpha)^2))
  
  ## Compute k_j^T * Grad^2 Sigma * k_j
  rho1_multiplier_2 <- apply(wave*Sigma_grad2_rho1%*%wave,2,sum)
  gamma_multiplier_2 <- apply(wave*Sigma_grad2_gamma%*%wave,2,sum)
  alpha_multiplier_2 <- apply(wave*Sigma_grad2_alpha%*%wave,2,sum)
  
  
  ## Compute -Delta*(k_j^T * Sigma * k_j + zeta)
  ExpTerm <- -dt*apply(wave*Sig%*%wave,2,sum)-dt*rep(zeta,dim(wave)[2])
  
  ## Compute Delta * mu * k_j
  TrigTerm <- dt*c(mu_x,mu_y)%*%wave
  
  ## Compute remaining gradients ##
  
  ## rho1 ##
  G_grad2[[4]]$G_grad2_cos <- (-dt*rho1_multiplier_2[1:n_cos]+dt^2*(rho1_multiplier[1:n_cos])^2)*exp(ExpTerm[1:n_cos])
  G_grad2[[4]]$G_grad2_1 <- (-dt*rho1_multiplier_2[cosine_indices]+dt^2*(rho1_multiplier[cosine_indices])^2)*exp(ExpTerm[cosine_indices])*cos(TrigTerm[cosine_indices])
  G_grad2[[4]]$G_grad2_2 <- -(-dt*rho1_multiplier_2[cosine_indices]+dt^2*(rho1_multiplier[cosine_indices])^2)*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices])
  
  ## gamma ##
  G_grad2[[5]]$G_grad2_cos <- (-dt*gamma_multiplier_2[1:n_cos]+dt^2*(gamma_multiplier[1:n_cos])^2)*exp(ExpTerm[1:n_cos])
  G_grad2[[5]]$G_grad2_1 <- (-dt*gamma_multiplier_2[cosine_indices]+dt^2*(gamma_multiplier[cosine_indices])^2)*exp(ExpTerm[cosine_indices])*cos(TrigTerm[cosine_indices])
  G_grad2[[5]]$G_grad2_2 <- -(-dt*gamma_multiplier_2[cosine_indices]+dt^2*(gamma_multiplier[cosine_indices])^2)*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices])
  
  ## alpha ##
  G_grad2[[6]]$G_grad2_cos <- (-dt*alpha_multiplier_2[1:n_cos]+dt^2*(alpha_multiplier[1:n_cos])^2)*exp(ExpTerm[1:n_cos])
  G_grad2[[6]]$G_grad2_1 <- (-dt*alpha_multiplier_2[cosine_indices]+dt^2*(alpha_multiplier[cosine_indices])^2)*exp(ExpTerm[cosine_indices])*cos(TrigTerm[cosine_indices])
  G_grad2[[6]]$G_grad2_2 <- -(-dt*alpha_multiplier_2[cosine_indices]+dt^2*(alpha_multiplier[cosine_indices])^2)*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices])
  
  ## mu_x ##
  G_grad2[[7]]$G_grad2_1 <- -dt^2*wave[1,cosine_indices]^2*exp(ExpTerm[cosine_indices])*cos(TrigTerm[cosine_indices])
  G_grad2[[7]]$G_grad2_2 <- dt^2*wave[1,cosine_indices]^2*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices])
  
  ## mu_y ##
  G_grad2[[8]]$G_grad2_1 <- -dt^2*wave[2,cosine_indices]^2*exp(ExpTerm[cosine_indices])*cos(TrigTerm[cosine_indices])
  G_grad2[[8]]$G_grad2_2 <- dt^2*wave[2,cosine_indices]^2*exp(ExpTerm[cosine_indices])*sin(TrigTerm[cosine_indices])
  
  return(G_grad2)
}

Q_grad2_f <- function(param_obs,param_sig,n,K,wave,n_cos=4,nu=1,dt=1,norm=TRUE,
                      param_bias = param_bias_func(),spec_weights_mat = NULL){
  
  ## Number of parameters
  n_param = param_obs$m + length(param_sig) + param_bias$m
  
  ## Required parameters
  rho0 = param_sig[1]; sigma2 = param_sig[2]; zeta = param_sig[3]; rho1 = param_sig[4]
  gamma = param_sig[5]; alpha = param_sig[6]; mu_x = param_sig[7]; mu_y = param_sig[8]
  
  ## Initialise list
  Q_grad2 = rep(list(matrix(0,nrow=K,ncol=K)),n_param)
  
  ## mu_x,mu_y,tau2 ##
  Q_grad2[[7]] = Q_grad2[[8]] = Q_grad2[[9]] = diag(0,nrow=K,ncol=K)
  
  ## rho0,sigma2,zeta,rho1,gamma,alpha ##
  
  ## Compute spectral density f(k_j) ##
  
  ## Compute k_j^T k_j
  w2 <- apply(wave^2,2,sum)
  
  ## Compute f(k_j)
  numerator <- 1
  denominator <- ((1/rho0)^2 + w2)^(nu + 1)
  spec <- numerator/denominator
  
  ## Adjustment for cosine only terms
  spec[1:n_cos] <- spec[1:n_cos]/2
  
  ## Normalise
  if(norm){
    spec_normalised <- spec*(n^2)/sum(spec)
  }
  else{
    spec_normalised <- spec/sum(spec)
  }
  
  ## Rescale by sigma^2
  spec_normalised <- sigma2 * spec_normalised
  
  ## Compute df(k_j)/d(sigma^2) & df(k_j)/d(rho_0)
  spec_normalised_grad_sigma2 <- 1/sigma2*spec_normalised
  spec_normalised_grad_rho0 <- 
    4/rho0^3*((w2+(1/rho0^2))^(-1)-sum((w2+(1/rho0^2))^(-1)*spec)/sum(spec))*spec_normalised
  
  ## Compute d^2f(k_j)/d(sigma^2)^2 & df^2(k_j)/d(rho_0)^2
  spec_normalised_grad2_sigma2 <- 0
  spec_normalised_grad2_rho0 <- 
    -12/rho0^4*((w2+(1/rho0^2))^(-1)-sum((w2+(1/rho0^2))^(-1)*spec)/sum(spec))*spec_normalised +
    24/rho0^6*((w2+(1/rho0^2))^(-2)-sum((w2+(1/rho0^2))^(-2)*spec)/sum(spec))*spec_normalised -
    32/rho0^6*((w2+(1/rho0^2))^(-1)*(sum((w2+(1/rho0^2))^(-1)*spec))/sum(spec)-(sum((w2+(1/rho0^2))^(-1)*spec)/sum(spec))^2)*spec_normalised
  
  
  ## Compute Sigma 
  if(rho1==0){
    Sig <- cbind(c(0,0),c(0,0))
  }
  else{
    Sig_Chol <- cbind(c(cos(alpha),-gamma*sin(alpha)),
                      c(sin(alpha),gamma*cos(alpha)))/rho1
    Sig <- solve(t(Sig_Chol)%*%Sig_Chol)
  }
  
  ## Compute Grad. Sigma
  Sigma_grad_rho1 <- 2 / rho1 * Sig
  
  Sigma_grad_gamma <- 
    2 * rho1^2 / gamma^3 * cbind(c(-sin(alpha)^2,sin(alpha)*cos(alpha)),
                                 c(sin(alpha)*cos(alpha),-cos(alpha)^2))
  
  Sigma_grad_alpha <- 
    rho1^2*(gamma^2-1) / gamma^2 * cbind(c(-2*sin(alpha)*cos(alpha),cos(alpha)^2-sin(alpha)^2),
                                         c(cos(alpha)^2-sin(alpha)^2,2*sin(alpha)*cos(alpha))) 
  
  ## Compute k_j^T * Grad Sigma * k_j
  rho1_multiplier <- apply(wave*Sigma_grad_rho1%*%wave,2,sum)
  gamma_multiplier <- apply(wave*Sigma_grad_gamma%*%wave,2,sum)
  alpha_multiplier <- apply(wave*Sigma_grad_alpha%*%wave,2,sum)
  
  
  ## Compute Grad^2 Sigma
  Sigma_grad2_rho1 <- 2 / rho1^2 * Sig
  
  Sigma_grad2_gamma <- 
    -6 * rho1^2 / gamma^4 * cbind(c(-sin(alpha)^2,sin(alpha)*cos(alpha)),
                                  c(sin(alpha)*cos(alpha),-cos(alpha)^2))
  
  Sigma_grad2_alpha <- 
    2*rho1^2*(gamma^2-1) / gamma^2 * cbind(c(sin(alpha)^2-cos(alpha)^2,-2*sin(alpha)*cos(alpha)),
                                           c(-2*sin(alpha)*cos(alpha),cos(alpha)^2-sin(alpha)^2))
  
  ## Compute k_j^T * Grad^2 Sigma * k_j
  rho1_multiplier_2 <- apply(wave*Sigma_grad2_rho1%*%wave,2,sum)
  gamma_multiplier_2 <- apply(wave*Sigma_grad2_gamma%*%wave,2,sum)
  alpha_multiplier_2 <- apply(wave*Sigma_grad2_alpha%*%wave,2,sum)
  
  
  ## Compute -Delta*(k_j^T * Sigma * k_j + zeta)
  ExpTerm <- -dt*apply(wave*Sig%*%wave,2,sum)-dt*rep(zeta,dim(wave)[2])
  
  
  ## Remaining Gradients ##
  
  ## rho0 ##
  Q_grad2[[1]] = diag(spec_normalised_grad2_rho0*(1-exp(2*dt*ExpTerm))/(-2*ExpTerm))
  
  ## sigma2: leave as zero ##
  
  ## zeta ##
  numerator_zeta = rep(1,dim(wave)[2])-(2*dt^2*(-ExpTerm)^2+2*dt*(-ExpTerm)+1)*exp(2*dt*ExpTerm)
  denominator_zeta = (-ExpTerm)^3
  Q_grad2[[3]] = diag(spec_normalised*numerator_zeta/denominator_zeta)
  
  ## rho1 ##
  numerator1_rho1 = ((2*dt*(-ExpTerm)+1)*exp(2*dt*ExpTerm)-rep(1,dim(wave)[2]))*rho1_multiplier_2
  denominator1_rho1 = 2*(-ExpTerm)^2
  numerator2_rho1 = (rep(1,dim(wave)[2])-(2*dt^2*(-ExpTerm)^2+2*dt*(-ExpTerm)+1)*exp(2*dt*ExpTerm))*rho1_multiplier^2
  denominator2_rho1 = (-ExpTerm)^3
  Q_grad2[[4]] = diag(spec_normalised*(numerator1_rho1/denominator1_rho1+numerator2_rho1/denominator2_rho1))
  
  ## gamma ##
  numerator1_gamma = ((2*dt*(-ExpTerm)+1)*exp(2*dt*ExpTerm)-rep(1,dim(wave)[2]))*gamma_multiplier_2
  denominator1_gamma = 2*(-ExpTerm)^2
  numerator2_gamma = (rep(1,dim(wave)[2])-(2*dt^2*(-ExpTerm)^2+2*dt*(-ExpTerm)+1)*exp(2*dt*ExpTerm))*gamma_multiplier^2
  denominator2_gamma = (-ExpTerm)^3
  Q_grad2[[5]] = diag(spec_normalised*(numerator1_gamma/denominator1_gamma+numerator2_gamma/denominator2_gamma))
  
  ## alpha ##
  numerator1_alpha = ((2*dt*(-ExpTerm)+1)*exp(2*dt*ExpTerm)-rep(1,dim(wave)[2]))*alpha_multiplier_2
  denominator1_alpha = 2*(-ExpTerm)^2
  numerator2_alpha = (rep(1,dim(wave)[2])-(2*dt^2*(-ExpTerm)^2+2*dt*(-ExpTerm)+1)*exp(2*dt*ExpTerm))*alpha_multiplier^2
  denominator2_alpha = (-ExpTerm)^3
  Q_grad2[[6]] = diag(spec_normalised*(numerator1_alpha/denominator1_alpha+numerator2_alpha/denominator2_alpha))
  
  
  ## Compute weighted matrix ##
  if(!is.null(spec_weights_mat)){
    Q_grad2 = lapply(1:n_param, function(i) spec_weights_mat%*%Q_grad2[[i]]%*%t(spec_weights_mat))
  }
  
  ## Output ##
  return(Q_grad2)
}

Q_grad2_f_vec <- function(param_obs,param_sig,n,K,wave,n_cos=4,nu=1,dt=1,
                              norm=TRUE,param_bias = param_bias_func()){
  
  ## Number of parameters
  n_param = param_obs$m + length(param_sig) + param_bias$m
  
  ## Required parameters
  rho0 = param_sig[1]; sigma2 = param_sig[2]; zeta = param_sig[3]; rho1 = param_sig[4]
  gamma = param_sig[5]; alpha = param_sig[6]; mu_x = param_sig[7]; mu_y = param_sig[8]
  
  ## Initialise array
  Q_grad2_vec = rep(list(rep(0,K)),n_param)
  
  ## mu_x,mu_y,tau2 : leave as zero ##
  
  ## rho0,sigma2,zeta,rho1,gamma,alpha ##
  
  ## Compute spectral density f(k_j) ##
  
  ## Compute k_j^T k_j
  w2 <- apply(wave^2,2,sum)
  
  ## Compute f(k_j)
  numerator <- 1
  denominator <- ((1/rho0)^2 + w2)^(nu + 1)
  spec <- numerator/denominator
  
  ## Adjustment for cosine only terms
  spec[1:n_cos] <- spec[1:n_cos]/2
  
  ## Normalise
  if(norm){
    spec_normalised <- spec*(n^2)/sum(spec)
  }
  else{
    spec_normalised <- spec/sum(spec)
  }
  
  ## Rescale by sigma^2
  spec_normalised <- sigma2 * spec_normalised
  
  ## Compute df(k_j)/d(sigma^2) & df(k_j)/d(rho_0)
  spec_normalised_grad_sigma2 <- 1/sigma2*spec_normalised
  spec_normalised_grad_rho0 <- 
    4/rho0^3*((w2+(1/rho0^2))^(-1)-sum((w2+(1/rho0^2))^(-1)*spec)/sum(spec))*spec_normalised
  
  ## Compute d^2f(k_j)/d(sigma^2)^2 & df^2(k_j)/d(rho_0)^2
  spec_normalised_grad2_sigma2 <- 0
  spec_normalised_grad2_rho0 <- 
    -12/rho0^4*((w2+(1/rho0^2))^(-1)-sum((w2+(1/rho0^2))^(-1)*spec)/sum(spec))*spec_normalised +
    24/rho0^6*((w2+(1/rho0^2))^(-2)-sum((w2+(1/rho0^2))^(-2)*spec)/sum(spec))*spec_normalised -
    32/rho0^6*((w2+(1/rho0^2))^(-1)*(sum((w2+(1/rho0^2))^(-1)*spec))/sum(spec)-(sum((w2+(1/rho0^2))^(-1)*spec)/sum(spec))^2)*spec_normalised
  
  
  ## Compute Sigma 
  if(rho1==0){
    Sig <- cbind(c(0,0),c(0,0))
  }
  else{
    Sig_Chol <- cbind(c(cos(alpha),-gamma*sin(alpha)),
                      c(sin(alpha),gamma*cos(alpha)))/rho1
    Sig <- solve(t(Sig_Chol)%*%Sig_Chol)
  }
  
  ## Compute Grad. Sigma
  Sigma_grad_rho1 <- 2 / rho1 * Sig
  
  Sigma_grad_gamma <- 
    2 * rho1^2 / gamma^3 * cbind(c(-sin(alpha)^2,sin(alpha)*cos(alpha)),
                                 c(sin(alpha)*cos(alpha),-cos(alpha)^2))
  
  Sigma_grad_alpha <- 
    rho1^2*(gamma^2-1) / gamma^2 * cbind(c(-2*sin(alpha)*cos(alpha),cos(alpha)^2-sin(alpha)^2),
                                         c(cos(alpha)^2-sin(alpha)^2,2*sin(alpha)*cos(alpha))) 
  
  ## Compute k_j^T * Grad Sigma * k_j
  rho1_multiplier <- apply(wave*Sigma_grad_rho1%*%wave,2,sum)
  gamma_multiplier <- apply(wave*Sigma_grad_gamma%*%wave,2,sum)
  alpha_multiplier <- apply(wave*Sigma_grad_alpha%*%wave,2,sum)
  
  
  ## Compute Grad^2 Sigma
  Sigma_grad2_rho1 <- 2 / rho1^2 * Sig
  
  Sigma_grad2_gamma <- 
    -6 * rho1^2 / gamma^4 * cbind(c(-sin(alpha)^2,sin(alpha)*cos(alpha)),
                                  c(sin(alpha)*cos(alpha),-cos(alpha)^2))
  
  Sigma_grad2_alpha <- 
    2*rho1^2*(gamma^2-1) / gamma^2 * cbind(c(sin(alpha)^2-cos(alpha)^2,-2*sin(alpha)*cos(alpha)),
                                           c(-2*sin(alpha)*cos(alpha),cos(alpha)^2-sin(alpha)^2))
  
  ## Compute k_j^T * Grad^2 Sigma * k_j
  rho1_multiplier_2 <- apply(wave*Sigma_grad2_rho1%*%wave,2,sum)
  gamma_multiplier_2 <- apply(wave*Sigma_grad2_gamma%*%wave,2,sum)
  alpha_multiplier_2 <- apply(wave*Sigma_grad2_alpha%*%wave,2,sum)
  
  
  ## Compute -Delta*(k_j^T * Sigma * k_j + zeta)
  ExpTerm <- -dt*apply(wave*Sig%*%wave,2,sum)-dt*rep(zeta,dim(wave)[2])
  
  
  ## Compute Remaining Gradients ##
  
  ## rho0 ##
  Q_grad2_vec[[1]] = spec_normalised_grad2_rho0*(1-exp(2*dt*ExpTerm))/(-2*ExpTerm)
  
  ## sigma2: leave as zero ##
  
  ## zeta ##
  numerator_zeta = rep(1,dim(wave)[2])-(2*dt^2*(-ExpTerm)^2+2*dt*(-ExpTerm)+1)*exp(2*dt*ExpTerm)
  denominator_zeta = (-ExpTerm)^3
  Q_grad2_vec[[3]] = spec_normalised*numerator_zeta/denominator_zeta
  
  ## rho1 ##
  numerator1_rho1 = ((2*dt*(-ExpTerm)+1)*exp(2*dt*ExpTerm)-rep(1,dim(wave)[2]))*rho1_multiplier_2
  denominator1_rho1 = 2*(-ExpTerm)^2
  numerator2_rho1 = (rep(1,dim(wave)[2])-(2*dt^2*(-ExpTerm)^2+2*dt*(-ExpTerm)+1)*exp(2*dt*ExpTerm))*rho1_multiplier^2
  denominator2_rho1 = (-ExpTerm)^3
  Q_grad2_vec[[4]] = spec_normalised*(numerator1_rho1/denominator1_rho1+numerator2_rho1/denominator2_rho1)
  
  ## gamma ##
  numerator1_gamma = ((2*dt*(-ExpTerm)+1)*exp(2*dt*ExpTerm)-rep(1,dim(wave)[2]))*gamma_multiplier_2
  denominator1_gamma = 2*(-ExpTerm)^2
  numerator2_gamma = (rep(1,dim(wave)[2])-(2*dt^2*(-ExpTerm)^2+2*dt*(-ExpTerm)+1)*exp(2*dt*ExpTerm))*gamma_multiplier^2
  denominator2_gamma = (-ExpTerm)^3
  Q_grad2_vec[[5]] = spec_normalised*(numerator1_gamma/denominator1_gamma+numerator2_gamma/denominator2_gamma)
  
  ## alpha ##
  numerator1_alpha = ((2*dt*(-ExpTerm)+1)*exp(2*dt*ExpTerm)-rep(1,dim(wave)[2]))*alpha_multiplier_2
  denominator1_alpha = 2*(-ExpTerm)^2
  numerator2_alpha = (rep(1,dim(wave)[2])-(2*dt^2*(-ExpTerm)^2+2*dt*(-ExpTerm)+1)*exp(2*dt*ExpTerm))*alpha_multiplier^2
  denominator2_alpha = (-ExpTerm)^3
  Q_grad2_vec[[6]] = spec_normalised*(numerator1_alpha/denominator1_alpha+numerator2_alpha/denominator2_alpha)
  
  ## Output ##
  return(Q_grad2_vec)
}

#######################

## ---

#########################
# TANGENT KALMAN FILTER #
#########################

TangentKalmanFilter <- function(y, lin_pred = matrix(0, dim(y)[1], dim(y)[2]), 
                                Phi, R, G, Q, Phi_grad, R_grad, G_grad, Q_grad,
                                n, t = dim(y)[1], 
                                K = dim(G)[1], log_lik = TRUE, mu0 = rep(0,K), 
                                Sigma0 = as.matrix(Q), mu0_grad = rep(0,K), 
                                Sigma0_grad = Q_grad, param,
                                lin_pred_grad = matrix(0,dim(y)[1],dim(y)[2])){
  
  ## Prelims
  Phi_T <- t(Phi)
  Phi_grad_T = t(Phi_grad)
  
  ## Run KalmanFilter ##
  KF = KalmanFilter(y = y, linear_pred = lin_pred, Phi = Phi, R = R, G = G,
                    Q = Q, N = dim(y)[2], t = t, K = K, log_lik = log_lik, m0 = mu0, 
                    R0 = Sigma0)
  
  ## KalmanFilter output ##
  mu_i_i1 = KF$m_i_i1
  mu_i_i = KF$m_i_i
  Sigma_i_i1 = KF$R_i_i1
  Sigma_i_i = KF$R_i_i
  innov_resid = KF$innov_resid
  innov_cov = KF$innov_cov
  ll = KF$ll
  
  ## Initialise gradient arrays ##
  
  # Predictions
  mu_i_i1_grad <- array(0, c(t, K)) # Mean Gradient
  Sigma_i_i1_grad <- array(0, c(t, K, K)) # Covariance Gradient
  
  # Updates
  mu_i_i_grad <- array(0, c(t+1, K)) # Mean Gradient
  Sigma_i_i_grad <- array(0, c(t+1, K, K)) # Covariance Gradient
  
  # Innovations
  innov_resid_grad <- array(0, c(t, dim(y)[2])) # Residual Gradient
  innov_cov_grad <- array(0, c(t, dim(y)[2], dim(y)[2])) # Covariance Gradient
  
  # Contional log-likelihood
  conditional_ll_grad <- array(0, t)
  
  # Log-likelihood
  ll_grad = 0
  
  
  ## Initial gradient values ##
  mu_i_i_grad[1,] <- mu0_grad # Mean
  Sigma_i_i_grad[1,,] <- Sigma0_grad # Covariance
  
  
  ## Tangent Kalman Filter ##
  
  for (i in 1:t) {
    
    ## Prediction
    
    ## -- Mean
    mu_i_i1_grad[i,] <- G %*% mu_i_i_grad[i,] + G_grad %*% mu_i_i[i,]
    
    ## -- Covariance
    Sigma_i_i1_grad[i,,] <- 
      G_grad %*% Sigma_i_i[i,,] %*% t(G) +
      G %*% Sigma_i_i_grad[i,,] %*% t(G) +
      G %*% Sigma_i_i[i,,] %*% t(G_grad) + 
      Q_grad
    Sigma_i_i1_grad[i,,] <- (Sigma_i_i1_grad[i,,] + t(Sigma_i_i1_grad[i,,]))/2 
    
    ## Update
    
    ## -- 'm'
    m_i_grad <- Phi_grad %*% mu_i_i1[i,] + Phi %*% mu_i_i1_grad[i,]
    
    ## -- Innovation Covariance
    S_i_grad <- 
      Phi_grad %*% Sigma_i_i1[i,,] %*% Phi_T +
      Phi %*% Sigma_i_i1_grad[i,,] %*% Phi_T + 
      Phi %*% Sigma_i_i1[i,,] %*% Phi_grad_T +
      R_grad
    S_i_grad_inv <- -innov_cov[i,,] %*% S_i_grad %*% innov_cov[i,,]
    innov_cov_grad[i,,] <- S_i_grad_inv
    
    ## -- Kalman Gain
    K_i <- Sigma_i_i1[i,,] %*% Phi_T %*% innov_cov[i,,]
    
    K_i_grad <- 
      Sigma_i_i1_grad[i,,] %*% Phi_T %*% innov_cov[i,,] +
      Sigma_i_i1[i,,] %*% Phi_grad_T %*% innov_cov[i,,] + 
      Sigma_i_i1[i,,] %*% Phi_T %*% innov_cov_grad[i,,]
    
    ## -- Innovation Residual
    innov_resid_grad[i,] <- - m_i_grad - lin_pred_grad[i,]
    
    ## -- Mean
    mu_i_i_grad[i+1,] <- 
      mu_i_i1_grad[i,] + 
      K_i_grad %*% innov_resid[i,] +
      K_i %*% innov_resid_grad[i,]
    
    ## -- Covariance
    Sigma_i_i_grad[i+1,,] <- 
      Sigma_i_i1_grad[i,,] - 
      K_i_grad %*% Phi %*% Sigma_i_i1[i,,] -
      K_i %*% Phi_grad %*% Sigma_i_i1[i,,] -
      K_i %*% Phi %*% Sigma_i_i1_grad[i,,]
    Sigma_i_i_grad[i+1,,] <- 
      (Sigma_i_i_grad[i+1,,] + t(Sigma_i_i_grad[i+1,,]))/2
    
    ## Log-likelihood gradient
    if (log_lik){
      
      conditional_ll_grad[i] = -
        0.5 * t(innov_resid_grad[i,]) %*% innov_cov[i,,] %*% innov_resid[i,] -
        0.5 * t(innov_resid[i,])  %*% innov_cov_grad[i,,] %*% innov_resid[i,] -
        0.5 * t(innov_resid[i,]) %*% innov_cov[i,,] %*% innov_resid_grad[i,] -
        0.5 * sum(diag(innov_cov[i,,] %*% S_i_grad))
      
      ll_grad <- ll_grad + conditional_ll_grad[i]
    }
  }
  
  ## Output
  output = list(mu_i_i = mu_i_i, Sigma_i_i = Sigma_i_i, mu_i_i1 = mu_i_i1, 
                Sigma_i_i1 = Sigma_i_i1, innov_resid = innov_resid, 
                innov_cov = innov_cov, K_i = K_i, mu_i_i_grad = mu_i_i_grad,
                Sigma_i_i_grad = Sigma_i_i_grad, mu_i_i1_grad = mu_i_i1_grad,
                Sigma_i_i1_grad = Sigma_i_i1_grad, 
                innov_resid_grad = innov_resid_grad, 
                innov_cov_grad = innov_cov_grad, K_i_grad = K_i_grad)
  
  if (log_lik){
    output = c(output, list(ll = ll, ll_grad = ll_grad,
                            conditional_ll_grad = conditional_ll_grad))
  }
  
  ## Print Output
  return(output)
}

MultiTangentKalmanFilter <- function(y, lin_pred = matrix(0, dim(y)[1], dim(y)[2]), 
                                     Phi, R, G, Q, Phi_grad, R_grad, G_grad, 
                                     Q_grad, n, t = dim(y)[1], 
                                     K = dim(G)[1], log_lik = TRUE, 
                                     mu0 = rep(0,K), Sigma0 = as.matrix(Q), 
                                     mu0_grad = rep(list(rep(0,K)),n_param),
                                     Sigma0_grad = Q_grad, n_param, 
                                     grad_params = 1:n_param, 
                                     lin_pred_grad = rep(list(lin_pred),n_param)){
  
  ## Preliminaries ##
  Phi_T <- t(Phi)
  Phi_grad_T = lapply(1:n_param, function(i) t(Phi_grad[[i]]))
  
  ## Run KalmanFilter ##
  KF = KalmanFilter(y = y, linear_pred = lin_pred, Phi = Phi, R = R, G = G,
                    Q = Q, N = dim(y)[2], t = t, K = K, log_lik = log_lik, m0 = mu0, 
                    R0 = Sigma0)
  
  ## KalmanFilter output ##
  mu_i_i1 = KF$m_i_i1
  mu_i_i = KF$m_i_i
  Sigma_i_i1 = KF$R_i_i1
  Sigma_i_i = KF$R_i_i
  innov_resid = KF$innov_resid
  innov_cov = KF$innov_cov
  ll = KF$ll
  
  ## Initialise gradient arrays ##
  
  # Predictions
  mu_i_i1_grad <- array(0, c(t, K)) # Mean Gradient
  Sigma_i_i1_grad <- array(0, c(t, K, K)) # Covariance Gradient
  
  mu_i_i1_grad = rep(list(mu_i_i1_grad),n_param)
  Sigma_i_i1_grad = rep(list(Sigma_i_i1_grad),n_param)
  
  # Updates
  mu_i_i_grad <- array(0, c(t+1, K)) # Mean Gradient
  Sigma_i_i_grad <- array(0, c(t+1, K, K)) # Covariance Gradient
  
  mu_i_i_grad = rep(list(mu_i_i_grad),n_param)
  Sigma_i_i_grad = rep(list(Sigma_i_i_grad),n_param)
  
  # Innovations
  innov_resid_grad <- array(0, c(t, dim(y)[2])) # Residual Gradient
  innov_cov_grad <- array(0, c(t, dim(y)[2], dim(y)[2])) # Covariance Gradient
  
  innov_resid_grad = rep(list(innov_resid_grad),n_param)
  innov_cov_grad = rep(list(innov_cov_grad),n_param)
  
  # Contional log-likelihood
  conditional_ll_grad <- array(0, t)
  conditional_ll_grad = rep(list(conditional_ll_grad),n_param)
  
  # Log-likelihood
  ll_grad = 0
  ll_grad = rep(list(ll_grad),n_param)
  
  
  ## Initial gradient values ##
  for (j in grad_params){
    mu_i_i_grad[[j]][1,] <- mu0_grad[[j]] # Mean
    Sigma_i_i_grad[[j]][1,,] <- Sigma0_grad[[j]] # Covariance
  }
  
  
  ## Tangent Kalman Filter ##
  for (j in grad_params){
    for (i in 1:t) {
      
      ## Prediction
      
      ## -- Mean
      mu_i_i1_grad[[j]][i,] <- 
        G %*% mu_i_i_grad[[j]][i,] + G_grad[[j]] %*% mu_i_i[i,]
      
      ## -- Covariance
      Sigma_i_i1_grad[[j]][i,,] <- 
        G_grad[[j]] %*% Sigma_i_i[i,,] %*% t(G) +
        G %*% Sigma_i_i_grad[[j]][i,,] %*% t(G) +
        G %*% Sigma_i_i[i,,] %*% t(G_grad[[j]]) + 
        Q_grad[[j]]
      Sigma_i_i1_grad[[j]][i,,] <- 
        (Sigma_i_i1_grad[[j]][i,,] + t(Sigma_i_i1_grad[[j]][i,,]))/2 
      
      ## Update
      
      ## -- 'm'
      m_i_grad <- 
        Phi_grad[[j]] %*% mu_i_i1[i,] + Phi %*% mu_i_i1_grad[[j]][i,]
      
      ## -- Innovation Covariance
      S_i_grad <- 
        Phi_grad[[j]] %*% Sigma_i_i1[i,,] %*% Phi_T +
        Phi %*% Sigma_i_i1_grad[[j]][i,,] %*% Phi_T + 
        Phi %*% Sigma_i_i1[i,,] %*% Phi_grad_T[[j]] +
        R_grad[[j]]
      S_i_grad_inv <- -innov_cov[i,,] %*% S_i_grad %*% innov_cov[i,,]
      innov_cov_grad[[j]][i,,] <- S_i_grad_inv
      
      ## -- Kalman Gain
      K_i <- Sigma_i_i1[i,,] %*% Phi_T %*% innov_cov[i,,]
      
      K_i_grad <- 
        Sigma_i_i1_grad[[j]][i,,] %*% Phi_T %*% innov_cov[i,,] +
        Sigma_i_i1[i,,] %*% Phi_grad_T[[j]] %*% innov_cov[i,,] + 
        Sigma_i_i1[i,,] %*% Phi_T %*% innov_cov_grad[[j]][i,,]
      
      ## -- Innovation Residual
      innov_resid_grad[[j]][i,] <- - m_i_grad - lin_pred_grad[[j]][i,]
      
      ## -- Mean
      mu_i_i_grad[[j]][i+1,] <- 
        mu_i_i1_grad[[j]][i,] + 
        K_i_grad %*% innov_resid[i,] +
        K_i %*% innov_resid_grad[[j]][i,]
      
      ## -- Covariance
      Sigma_i_i_grad[[j]][i+1,,] <- 
        Sigma_i_i1_grad[[j]][i,,] - 
        K_i_grad %*% Phi %*% Sigma_i_i1[i,,] -
        K_i %*% Phi_grad[[j]] %*% Sigma_i_i1[i,,] -
        K_i %*% Phi %*% Sigma_i_i1_grad[[j]][i,,]
      Sigma_i_i_grad[[j]][i+1,,] <- 
        (Sigma_i_i_grad[[j]][i+1,,] + t(Sigma_i_i_grad[[j]][i+1,,]))/2
      
      ## Log-likelihood gradient
      if (log_lik){
        conditional_ll_grad[[j]][i] = -
          0.5 * t(innov_resid_grad[[j]][i,]) %*% innov_cov[i,,] %*% innov_resid[i,] -
          0.5 * t(innov_resid[i,])  %*% S_i_grad_inv %*% innov_resid[i,] -
          0.5 * t(innov_resid[i,]) %*% innov_cov[i,,] %*% innov_resid_grad[[j]][i,] -
          0.5 * sum(diag(innov_cov[i,,] %*% S_i_grad))
        
        ll_grad[[j]] <- ll_grad[[j]] + conditional_ll_grad[[j]][i]
      }
    }
  }
  
  ## Output
  output = list(mu_i_i = mu_i_i, Sigma_i_i = Sigma_i_i, mu_i_i1 = mu_i_i1, 
                Sigma_i_i1 = Sigma_i_i1, innov_resid = innov_resid, 
                innov_cov = innov_cov, K_i = K_i, mu_i_i_grad = mu_i_i_grad,
                Sigma_i_i_grad = Sigma_i_i_grad, mu_i_i1_grad = mu_i_i1_grad,
                Sigma_i_i1_grad = Sigma_i_i1_grad, 
                innov_resid_grad = innov_resid_grad, 
                innov_cov_grad = innov_cov_grad, K_i_grad = K_i_grad)
  
  if (log_lik){
    output = c(output, list(ll = ll, ll_grad = ll_grad,
                            conditional_ll_grad = conditional_ll_grad))
  }
  
  ## Print Output
  return(output)
}

#########################

##################################
# SPECTRAL TANGENT KALMAN FILTER #
##################################

TKF_spectral <- function(w = NULL, w_FFT_real = NULL, param_obs = NULL, 
                         param_sig = NULL, n, t, K=n^2, log_lik = TRUE, 
                         dt = 1, m0 = rep(0,K), R0 = rep(0,K),
                         m0_grad = rep(0,K), R0_grad = rep(0,K), nu = 1,
                         grad_params, log_lik_grad = TRUE){
  
  ## If 'w_FFT_real' not given, compute as inverse-FFT of w
  if(is.null(w_FFT_real)){
    w_FFT_real <- FFT_real_ST(STmatrix_to_vec(w),n=n,t=t,inv=TRUE)
  }
  
  ## If 'w_FFT_real' given as matrix,convert into vector
  else{
    if(class(w_FFT_real)=="matrix"){
      w_FFT_real <- STmatrix_to_vec(w_FFT_real)
    }
  }
  
  ## Select required terms from w_FFT_real
  if(K<n^2){
    basis_indices <- spde_initialise(n,t,K)$basis_indices
    w_FFT_indices <- rep(seq(0,n^2*(t-1),n^2),each=K)+basis_indices
    w_FFT_real = w_FFT_real[w_FFT_indices]
  }
  
  ## Prelims
  spde_ft = spde_initialise(n,t,K)
  wave = spde_ft$wave
  cosine_indices = spde_ft$cosine_indices
  n_cos = spde_ft$n_cos
  n_cos_sin = length(cosine_indices)
  K = spde_ft$K
  
  
  ## SPDE matrices
  G_vec <- G_f_vec(param_sig,wave,cosine_indices,dt=dt,n_cos)
  Q_vec <- Q_f_vec(param_sig,wave,n,n_cos,nu,dt=dt,norm=TRUE)
  
  G_1 = G_vec$G_1; G_2 = G_vec$G_2; G_cos = G_vec$G_cos
  Q_vec_cos = Q_vec[1:n_cos]; Q_vec_cos_sin = Q_vec[cosine_indices]
  tau2 <- param_obs$tau2_vals[[1]]
  
  
  ## SPDE matrix gradients
  R_grad_vec = R_grad_f_vec(param_obs,param_sig,n,n_obs=n^2)[[grad_param]]
  G_grad_vec = G_grad_f_vec(param_obs,param_sig,n,K,wave,cosine_indices,dt=dt,nu,n_cos)[[grad_param]]
  Q_grad_vec = Q_grad_f_vec(param_obs,param_sig,n,K,wave,n_cos,nu,dt=dt,norm=TRUE)[[grad_param]]
  
  G_grad_cos = G_grad_vec$G_grad_cos
  G_grad_1 = G_grad_vec$G_grad_1
  G_grad_2 = G_grad_vec$G_grad_2
  
  
  ## Kalman Filter 
  KF <- kf_spectral(w_FFT_real = w_FFT_real,m_i_i1 = rep(0,K*t),
                    m_i_i=c(m0,rep(0,K*t)),R_i_i1=rep(0,K*t),
                    R_i_i=c(R0,rep(0,K*t)),G_cos = G_cos, 
                    G_1 = G_1,G_2 = G_2,spectrum = Q_vec, 
                    t = t, n_cos = n_cos,  n_cos_sin = n_cos_sin, 
                    tau2 = tau2)
  
  
  
  ## Tangent Kalman Filter
  TKF = tkf_spectral(w_FFT_real,KF$m_i_i1,KF$m_i_i,KF$R_i_i1,KF$R_i_i,
                     m_i_i1_grad=rep(0,K*t),
                     m_i_i_grad=c(m0_grad,rep(0,K*t)),
                     R_i_i1_grad=rep(0,K*t),
                     R_i_i_grad=c(R0_grad,rep(0,K*t)),
                     G_cos,G_1,G_2,G_grad_cos,G_grad_1,G_grad_2,Q_vec,
                     Q_grad_vec,t, n_cos, n_cos_sin, tau2, R_grad_vec)
  
  ## Prepare output
  mu_i_i = vec_to_STmatrix(KF$m_i_i,t=t+1)
  mu_i_i1 = vec_to_STmatrix(KF$m_i_i1,t=t)
  Sigma_i_i = vec_to_STmatrix(KF$R_i_i,t=t+1)
  Sigma_i_i1 = vec_to_STmatrix(KF$R_i_i1,t=t)
  
  mu_i_i_grad = vec_to_STmatrix(TKF$mu_i_i_grad,t=t+1)
  mu_i_i1_grad = vec_to_STmatrix(TKF$mu_i_i1_grad,t=t)
  Sigma_i_i_grad = vec_to_STmatrix(TKF$Sigma_i_i_grad,t=t+1)
  Sigma_i_i1_grad = vec_to_STmatrix(TKF$Sigma_i_i1_grad,t=t)
  
  
  ## Output
  output <- list(mu_i_i = mu_i_i, Sigma_i_i = Sigma_i_i, 
                 mu_i_i1 = mu_i_i1, Sigma_i_i1 = Sigma_i_i1,
                 mu_i_i_grad = mu_i_i_grad, Sigma_i_i_grad = Sigma_i_i_grad, 
                 mu_i_i1_grad = mu_i_i1_grad, Sigma_i_i1_grad = Sigma_i_i1_grad)
  
  ## If 'log_lik'=TRUE: compute likelihood  & amend output
  if(log_lik){
    
    if(K==n^2){
      ll = ll_spectral(log_lik=0,w_FFT_real, KF$m_i_i1, KF$R_i_i1, t, K, 
                       tau2)$ll
    }
    
    if(K<n^2){
      ll = ll_non_spectral(log_lik=0,w,KF$m_i_i1,KF$R_i_i1,t,n,K,tau2)$ll
    }
    
    output = c(output, list(ll=ll))
  }
  
  ## If 'log_lik_grad'=TRUE: compute grad log-likelihood & amend output
  if(log_lik_grad){
    if(K==n^2){
      loglik_grad = ll_grad_spectral(log_lik_grad=0,
                                     conditional_log_lik_grad = rep(0,t),
                                     w_FFT_real,KF$m_i_i1,KF$R_i_i1,
                                     TKF$mu_i_i1_grad,TKF$Sigma_i_i1_grad,t,K,
                                     tau2,R_grad_vec)
      
      ll_grad = loglik_grad$ll_grad
      conditional_ll_grad = loglik_grad$conditional_ll_grad
    }
    
    if(K<n^2){
      loglik_grad = ll_grad_non_spectral(log_lik_grad=0,
                                         conditional_log_lik_grad=rep(0,t),
                                         w=w,m_i_i1=KF$m_i_i1,R_i_i1=KF$R_i_i1,
                                         m_i_i1_grad=TKF$mu_i_i1_grad,
                                         R_i_i1_grad=TKF$Sigma_i_i1_grad,
                                         t=t,n=n,K=K,tau2=tau2,
                                         R_grad = diag(R_grad_vec))
      ll_grad = loglik_grad$ll_grad
      conditional_ll_grad = loglik_grad$conditional_ll_grad
    }
    
    
    output = c(output,list(ll_grad=ll_grad,
                           conditional_ll_grad = conditional_ll_grad))
  }
  
  ## Return output
  return(output)
}

Multi_TKF_spectral <- function(w = NULL, w_FFT_real = NULL, param_obs = NULL, 
                               param_sig = NULL, n, t, K=n^2, log_lik = TRUE,
                               dt = 1, m0 = rep(0,K), R0 = rep(0,K),
                               m0_grad = rep(list(rep(0,K)),9), 
                               R0_grad = rep(list(rep(0,K)),9), 
                               nu = 1, grad_params, log_lik_grad = TRUE){
  
  ## If 'w_FFT_real' not given, compute as inverse-FFT of w
  if(is.null(w_FFT_real)){
    w_FFT_real <- FFT_real_ST(STmatrix_to_vec(w),n=n,t=t,inv=TRUE)
  }
  
  ## If 'w_FFT_real' given as matrix,convert into vector
  else{
    if(class(w_FFT_real)=="matrix"){
      w_FFT_real <- STmatrix_to_vec(w_FFT_real)
    }
  }
  
  ## Select required terms from w_FFT_real
  if(K<n^2){
    basis_indices <- spde_initialise(n,t,K)$basis_indices
    w_FFT_indices <- rep(seq(0,n^2*(t-1),n^2),each=K)+basis_indices
    w_FFT_real = w_FFT_real[w_FFT_indices]
  }
  
  ## Prelims
  spde_ft = spde_initialise(n,t,K)
  wave = spde_ft$wave
  cosine_indices = spde_ft$cosine_indices
  n_cos = spde_ft$n_cos
  n_cos_sin = length(cosine_indices)
  K = spde_ft$K
  
  
  ## SPDE matrices
  G_vec <- G_f_vec(param_sig,wave,cosine_indices,dt=dt,n_cos)
  Q_vec <- Q_f_vec(param_sig,wave,n,n_cos,nu,dt=dt,norm=TRUE)
  
  G_1 = G_vec$G_1; G_2 = G_vec$G_2; G_cos = G_vec$G_cos
  Q_vec_cos = Q_vec[1:n_cos]; Q_vec_cos_sin = Q_vec[cosine_indices]
  tau2 <- param_obs$tau2_vals[[1]]
  
  ## Kalman Filter 
  KF <- kf_spectral(w_FFT_real = w_FFT_real,m_i_i1 = rep(0,K*t),
                    m_i_i=c(m0,rep(0,K*t)),R_i_i1=rep(0,K*t),
                    R_i_i=c(R0,rep(0,K*t)),G_cos = G_cos, 
                    G_1 = G_1,G_2 = G_2,spectrum = Q_vec, 
                    t = t, n_cos = n_cos, n_cos_sin = n_cos_sin, 
                    tau2 = tau2)
  
  ## Kalman Filter Outputs
  mu_i_i = vec_to_STmatrix(KF$m_i_i,t=t+1)
  mu_i_i1 = vec_to_STmatrix(KF$m_i_i1,t=t)
  Sigma_i_i = vec_to_STmatrix(KF$R_i_i,t=t+1)
  Sigma_i_i1 = vec_to_STmatrix(KF$R_i_i1,t=t)
  
  output = list(mu_i_i = mu_i_i, Sigma_i_i = Sigma_i_i,mu_i_i1 = mu_i_i1, 
                Sigma_i_i1 = Sigma_i_i1)
  
  
  ## Initialise gradient arrays
  
  # Predictions
  mu_i_i1_grad = rep(list(array(0, c(t, K))),9) # Mean gradients
  Sigma_i_i1_grad = rep(list(array(0, c(t, K))),9) # Cov. gradients
  
  # Updates
  mu_i_i_grad = rep(list(array(0, c(t+1, K))),9) # Mean gradients
  Sigma_i_i_grad = rep(list(array(0, c(t+1, K))),9) #C ov. gradients
  
  # Log-likelihood
  ll_grad = rep(list(0),9)
  conditional_ll_grad = rep(list(rep(0,t)),9)
  
  
  ## SPDE matrix gradients
  R_grad_vec = R_grad_f_vec(param_obs,param_sig,n,n_obs=n^2)
  G_grad_vec = G_grad_f_vec(param_obs,param_sig,n,K,wave,cosine_indices,dt=dt,nu,n_cos)
  Q_grad_vec = Q_grad_f_vec(param_obs,param_sig,n,K,wave,n_cos,nu,dt=dt,norm=TRUE)
  
  G_grad_cos = lapply(1:9, function(i) G_grad_vec[[i]]$G_grad_cos)
  G_grad_1 = lapply(1:9, function(i) G_grad_vec[[i]]$G_grad_1)
  G_grad_2 = lapply(1:9, function(i) G_grad_vec[[i]]$G_grad_2)
  
  
  for (j in grad_params){
    
    ## Tangent Kalman Filter
    TKF = tkf_spectral(w_FFT_real,KF$m_i_i1,KF$m_i_i,KF$R_i_i1,KF$R_i_i,
                       m_i_i1_grad=rep(0,K*t),
                       m_i_i_grad=c(m0_grad[[j]],rep(0,K*t)),
                       R_i_i1_grad=rep(0,K*t),
                       R_i_i_grad=c(R0_grad[[j]],rep(0,K*t)),
                       G_cos,G_1,G_2,G_grad_cos[[j]],G_grad_1[[j]],
                       G_grad_2[[j]],Q_vec,Q_grad_vec[[j]],t, n_cos, n_cos_sin, 
                       tau2, R_grad_vec[[j]])
    
    mu_i_i_grad[[j]] = vec_to_STmatrix(TKF$mu_i_i_grad,t=t+1)
    mu_i_i1_grad[[j]] = vec_to_STmatrix(TKF$mu_i_i1_grad,t=t)
    Sigma_i_i_grad[[j]] = vec_to_STmatrix(TKF$Sigma_i_i_grad,t=t+1)
    Sigma_i_i1_grad[[j]] = vec_to_STmatrix(TKF$Sigma_i_i1_grad,t=t)
    
    if (log_lik_grad){
      if(K==n^2){
        loglik_grad = ll_grad_spectral(log_lik_grad=0,
                                       conditional_log_lik_grad = rep(0,t),
                                       w_FFT_real,KF$m_i_i1,KF$R_i_i1,
                                       TKF$mu_i_i1_grad,TKF$Sigma_i_i1_grad,t,K,
                                       tau2,R_grad_vec[[j]])
        
        ll_grad[[j]] = loglik_grad$ll_grad
        conditional_ll_grad[[j]] = loglik_grad$conditional_ll_grad
      }
      
      if(K<n^2){
        loglik_grad = ll_grad_non_spectral(log_lik_grad=0,
                                           conditional_log_lik_grad=rep(0,t),
                                           w=w,m_i_i1=KF$m_i_i1,R_i_i1=KF$R_i_i1,
                                           m_i_i1_grad=TKF$mu_i_i1_grad,
                                           R_i_i1_grad=TKF$Sigma_i_i1_grad,
                                           t=t,n=n,K=K,tau2=tau2,
                                           R_grad = diag(R_grad_vec[[j]]))
        ll_grad[[j]] = loglik_grad$ll_grad
        conditional_ll_grad[[j]] = loglik_grad$conditional_ll_grad
      }
    }
  }
  
  
  ## Output
  output <- c(output, (list(mu_i_i_grad = mu_i_i_grad, 
                            Sigma_i_i_grad = Sigma_i_i_grad, 
                            mu_i_i1_grad = mu_i_i1_grad, 
                            Sigma_i_i1_grad = Sigma_i_i1_grad)))
  
  ## If 'log_lik'=TRUE: compute likelihood  & amend output
  if(log_lik){
    
    if(K==n^2){
      ll = ll_spectral(log_lik=0,w_FFT_real, KF$m_i_i1, KF$R_i_i1, t, K, 
                       tau2)$ll
    }
    
    if(K<n^2){
      ll = ll_non_spectral(log_lik=0,w,KF$m_i_i1,KF$R_i_i1,t,n,K,tau2)$ll
    }
    
    output = c(output, list(ll=ll))
  }
  
  ## If 'log_lik_grad'=TRUE: amend output
  if(log_lik_grad){
    output = c(output,list(ll_grad=ll_grad,
                           conditional_ll_grad = conditional_ll_grad))
  }
  
  ## Return output
  return(output)
}

tkf_spectral = function(w_FFT_real,m_i_i1,m_i_i,R_i_i1,R_i_i,
                        m_i_i1_grad,m_i_i_grad,R_i_i1_grad,R_i_i_grad,
                        G_cos,G_1,G_2,G_grad_cos,G_grad_1,G_grad_2,Q_vec,
                        Q_grad_vec, t, n_cos, n_cos_sin, tau2, R_grad){
  
  tkf = .C("tkf_spectral", y_FFT = as.double(w_FFT_real),
           mu_i_i1=as.double(m_i_i1),mu_i_i=as.double(m_i_i),
           Sigma_i_i1=as.double(R_i_i1),Sigma_i_i = as.double(R_i_i), 
           mu_i_i1_grad = as.double(m_i_i1_grad),
           mu_i_i_grad = as.double(m_i_i_grad),
           Sigma_i_i1_grad = as.double(R_i_i1_grad), 
           Sigma_i_i_grad = as.double(R_i_i_grad), 
           G_cos = as.double(G_cos), G_1 = as.double(G_1), 
           G_2 = as.double(G_2),G_grad_cos = as.double(G_grad_cos), 
           G_grad_1 = as.double(G_grad_1),G_grad_2 = as.double(G_grad_2), 
           Q = as.double(Q_vec), Q_grad = as.double(Q_grad_vec),
           T = as.integer(t), n_cos = as.integer(n_cos), 
           n_cos_sin=as.integer(n_cos_sin),tau2 = as.double(tau2), 
           R_grad = as.double(R_grad))
  
  
  output = list(mu_i_i1_grad = tkf$mu_i_i1_grad, 
                mu_i_i_grad = tkf$mu_i_i_grad, 
                Sigma_i_i1_grad = tkf$Sigma_i_i1_grad,
                Sigma_i_i_grad = tkf$Sigma_i_i_grad)
  
  return(output)
}

##################################

## ---

#############################
# 2ND TANGENT KALMAN FILTER #
#############################

TangentKalmanFilter2 <- function(y, lin_pred = matrix(0, dim(y)[1], dim(y)[2]),
                                 Phi, R, G, Q, Phi_grad, R_grad, G_grad, Q_grad,
                                 Phi_grad2, R_grad2, G_grad2, Q_grad2,
                                 n, t = dim(y)[1], 
                                 K = dim(G)[1], log_lik = TRUE, mu0 = rep(0,K),
                                 Sigma0 = as.matrix(Q), mu0_grad = rep(0,K),
                                 Sigma0_grad = Q_grad, mu0_grad2 = rep(0,K),
                                 Sigma0_grad2 = Q_grad2, param,
                                 lin_pred_grad = matrix(0,dim(y)[1],dim(y)[2]),
                                 lin_pred_grad2 = matrix(0,dim(y)[1],dim(y)[2])){
  
  ## Prelims
  Phi_T <- t(Phi)
  Phi_grad_T = t(Phi_grad)
  Phi_grad2_T = t(Phi_grad2)
  
  ## Run KalmanFilter ##
  KF = KalmanFilter(y = y, linear_pred = lin_pred, Phi = Phi, R = R, G = G,
                    Q = Q, N = dim(y)[2], t = t, K = K, log_lik = log_lik, m0 = mu0, 
                    R0 = Sigma0)
  
  ## KalmanFilter output ##
  mu_i_i1 = KF$m_i_i1
  mu_i_i = KF$m_i_i
  Sigma_i_i1 = KF$R_i_i1
  Sigma_i_i = KF$R_i_i
  innov_resid = KF$innov_resid
  innov_cov = KF$innov_cov
  ll = KF$ll
  
  ## Run Tangent Kalman Filter
  TKF = TangentKalmanFilter(y = y, lin_pred = lin_pred, Phi = Phi, R = R, G = G,
                            Q = Q, Phi_grad = Phi_grad, R_grad = R_grad,
                            G_grad = G_grad, Q_grad = Q_grad, n = n, t = t, 
                            K = K, log_lik = log_lik, mu0 = mu0, Sigma0 = Sigma0,
                            mu0_grad = mu0_grad, Sigma0_grad = Sigma0_grad, 
                            param = param, lin_pred_grad = lin_pred_grad)
  
  mu_i_i1_grad = TKF$mu_i_i1_grad
  mu_i_i_grad = TKF$mu_i_i_grad
  Sigma_i_i1_grad = TKF$Sigma_i_i1_grad
  Sigma_i_i_grad = TKF$Sigma_i_i_grad
  innov_resid_grad = TKF$innov_resid_grad
  innov_cov_grad = TKF$innov_cov_grad
  ll_grad = TKF$ll_grad   
  conditional_ll_grad = TKF$conditional_ll_grad
  
  ## Initialise 2nd Derivative Arrays ##
  
  # Predictions
  mu_i_i1_grad2 <- array(0, c(t, K)) # Mean 2nd Derivative 
  Sigma_i_i1_grad2 <- array(0, c(t, K, K)) # Covariance 2nd Derivative
  
  # Updates
  mu_i_i_grad2 <- array(0, c(t+1, K)) # Mean 2nd Derivative
  Sigma_i_i_grad2 <- array(0, c(t+1, K, K)) # Covariance 2nd Derivative
  
  # Innovations
  innov_resid_grad2 <- array(0, c(t, n^2)) # Residual 2nd Derivative
  innov_cov_grad2 <- array(0, c(t, n^2, n^2)) # Covariance 2nd Derivative
  
  # Contional log-likelihood
  conditional_ll_grad2 <- array(0, t)
  
  # Log-likelihood
  ll_grad2 = 0
  
  
  ## Initial 2nd derivative values ##
  mu_i_i_grad2[1,] <- mu0_grad2 # Mean
  Sigma_i_i_grad2[1,,] <- Sigma0_grad2 # Covariance
  
  
  ## Second Order Tangent Kalman Filter ##
  
  for (i in 1:t) {
    
    ## Prediction
    
    ## -- Mean
    mu_i_i1_grad2[i,] <- 
      G_grad2 %*% mu_i_i[i,] + 
      2 * G_grad %*% mu_i_i_grad[i,] +
      G %*% mu_i_i_grad2[i,]
    
    ## -- Covariance
    Sigma_i_i1_grad2[i,,] <- 
      G_grad2 %*% Sigma_i_i[i,,] %*% t(G) + 
      G %*% Sigma_i_i_grad2[i,,] %*% t(G) +
      G %*% Sigma_i_i[i,,] %*% t(G_grad2) + 
      2 * G_grad %*% Sigma_i_i_grad[i,,] %*% t(G) + 
      2 * G_grad %*% Sigma_i_i[i,,] %*% t(G_grad) + 
      2 * G %*% Sigma_i_i_grad[i,,] %*% t(G_grad) + 
      Q_grad2
    
    Sigma_i_i1_grad2[i,,] <- 
      (Sigma_i_i1_grad2[i,,] + t(Sigma_i_i1_grad2[i,,]))/2 
    
    ## Update
    
    ## -- 'm'
    m_i_grad2 <- 
      Phi_grad2 %*% mu_i_i1[i,] + 
      2 * Phi_grad %*% mu_i_i1_grad[i,] +
      Phi %*% mu_i_i1_grad2[i,]
    
    ## -- Innovation Covariance
    S_i <- solve(innov_cov[i,,])
    
    S_i_grad2 <- 
      Phi_grad2 %*% Sigma_i_i1[i,,] %*% Phi_T +
      Phi %*% Sigma_i_i1_grad2[i,,] %*% Phi_T + 
      Phi %*% Sigma_i_i1[i,,] %*% Phi_grad2_T +
      2 * Phi_grad %*% Sigma_i_i1_grad[i,,] %*% Phi_T + 
      2 * Phi_grad %*% Sigma_i_i1[i,,] %*% Phi_grad_T + 
      2 * Phi %*% Sigma_i_i1_grad[i,,] %*% Phi_grad_T + 
      R_grad2
    
    S_i_grad2_inv <- 
      2 * innov_cov_grad[i,,] %*% S_i %*% innov_cov_grad[i,,] -
      innov_cov[i,,] %*% S_i_grad2 %*% innov_cov[i,,]
    
    innov_cov_grad2[i,,] <- S_i_grad2_inv
    
    ## -- Kalman Gain
    K_i <- Sigma_i_i1[i,,] %*% Phi_T %*% innov_cov[i,,]
    
    K_i_grad <- 
      Sigma_i_i1_grad[i,,] %*% Phi_T %*% innov_cov[i,,] +
      Sigma_i_i1[i,,] %*% Phi_grad_T %*% innov_cov[i,,] + 
      Sigma_i_i1[i,,] %*% Phi_T %*% innov_cov_grad[i,,]
    
    K_i_grad2 <- 
      Sigma_i_i1_grad2[i,,] %*% Phi_T %*% innov_cov[i,,] +
      Sigma_i_i1[i,,] %*% Phi_grad2_T %*% innov_cov[i,,] + 
      Sigma_i_i1[i,,] %*% Phi_T %*% innov_cov_grad2[i,,] +
      2 * Sigma_i_i1_grad[i,,] %*% Phi_grad_T %*% innov_cov[i,,] + 
      2 * Sigma_i_i1_grad[i,,] %*% Phi_T %*% innov_cov_grad[i,,] + 
      2 * Sigma_i_i1[i,,] %*% Phi_grad_T %*% innov_cov_grad[i,,]
      
    
    ## -- Innovation Residual
    innov_resid_grad2[i,] <- - m_i_grad2 - lin_pred_grad2[i,]
    
    ## -- Mean
    mu_i_i_grad2[i+1,] <- 
      mu_i_i1_grad2[i,] + 
      K_i_grad2 %*% innov_resid[i,] + 
      2 * K_i_grad %*% innov_resid_grad[i,] +
      K_i %*% innov_resid_grad2[i,]
    
    ## -- Covariance
    Sigma_i_i_grad2[i+1,,] <- 
      Sigma_i_i1_grad2[i,,] -
      K_i_grad2 %*% Phi %*% Sigma_i_i1[i,,] -
      K_i %*% Phi_grad2 %*% Sigma_i_i1[i,,] -
      K_i %*% Phi %*% Sigma_i_i1_grad2[i,,] -
      2 * K_i_grad %*% Phi_grad %*% Sigma_i_i1[i,,] -
      2 * K_i_grad %*% Phi %*% Sigma_i_i1_grad[i,,] - 
      2 * K_i %*% Phi_grad %*% Sigma_i_i1_grad[i,,]

    Sigma_i_i_grad2[i+1,,] <- 
      (Sigma_i_i_grad2[i+1,,] + t(Sigma_i_i_grad2[i+1,,]))/2
    
    ## Log-likelihood gradient
    if (log_lik){
      
      conditional_ll_grad2[i] = -
        0.5 * t(innov_resid_grad2[i,]) %*% innov_cov[i,,] %*% innov_resid[i,] -
        0.5 * t(innov_resid[i,]) %*% innov_cov_grad2[i,,] %*% innov_resid[i,] -
        0.5 * t(innov_resid[i,]) %*% innov_cov[i,,] %*% innov_resid_grad2[i,] -
        t(innov_resid_grad[i,]) %*% innov_cov_grad[i,,] %*% innov_resid[i,] -
        t(innov_resid_grad[i,]) %*% innov_cov[i,,] %*% innov_resid_grad[i,] -
        t(innov_resid[i,]) %*% innov_cov_grad[i,,] %*% innov_resid_grad[i,] -
        0.5 * sum(diag(innov_cov[i,,] %*% S_i_grad2)) +
        0.5 * sum(diag(innov_cov_grad[i,,] %*% S_i %*% innov_cov_grad[i,,] %*% S_i))
        
    ll_grad2 <- ll_grad2 + conditional_ll_grad2[i]
    }
  }
  
  ## Output
  output = list(mu_i_i = mu_i_i, Sigma_i_i = Sigma_i_i, mu_i_i1 = mu_i_i1, 
                Sigma_i_i1 = Sigma_i_i1, innov_resid = innov_resid, 
                innov_cov = innov_cov, K_i = K_i, mu_i_i_grad = mu_i_i_grad,
                Sigma_i_i_grad = Sigma_i_i_grad, mu_i_i1_grad = mu_i_i1_grad,
                Sigma_i_i1_grad = Sigma_i_i1_grad, 
                innov_resid_grad = innov_resid_grad, 
                innov_cov_grad = innov_cov_grad, K_i_grad = K_i_grad,
                mu_i_i_grad2 = mu_i_i_grad2,
                Sigma_i_i_grad2 = Sigma_i_i_grad2, 
                mu_i_i1_grad2 = mu_i_i1_grad2,
                Sigma_i_i1_grad2 = Sigma_i_i1_grad2, 
                innov_resid_grad2 = innov_resid_grad2, 
                innov_cov_grad2 = innov_cov_grad2, K_i_grad2 = K_i_grad2)
  
  if (log_lik){
    output = c(output, list(ll = ll, ll_grad = ll_grad, 
                            conditional_ll_grad = conditional_ll_grad,
                            ll_grad2 = ll_grad2,
                            conditional_ll_grad2 = conditional_ll_grad2))
  }
  
  ## Print Output
  return(output)
}

MultiTangentKalmanFilter2 <- function(y, lin_pred = matrix(0, dim(y)[1], dim(y)[2]), 
                                      Phi, R, G, Q, Phi_grad, R_grad, G_grad, 
                                      Q_grad, Phi_grad2, R_grad2, G_grad2, 
                                      Q_grad2, n, t = dim(y)[1], 
                                      K = dim(G)[1], log_lik = TRUE, 
                                      mu0 = rep(0,K), Sigma0 = as.matrix(Q), 
                                      mu0_grad = rep(list(rep(0,K)),n_param),
                                      Sigma0_grad = Q_grad, 
                                      mu0_grad2 = rep(list(rep(0,K)),n_param),
                                      Sigma0_grad2 = Q_grad2, n_param,
                                      grad_params = 1:n_param,
                                      lin_pred_grad = rep(list(lin_pred),n_param),
                                      lin_pred_grad2 = rep(list(lin_pred),n_param)){
  
  ## Preliminaries ##
  Phi_T <- t(Phi)
  Phi_grad_T = lapply(1:n_param, function(i) t(Phi_grad[[i]]))
  Phi_grad2_T = lapply(1:n_param, function(i) t(Phi_grad2[[i]]))
  
  ## Run KalmanFilter ##
  KF = KalmanFilter(y = y, linear_pred = lin_pred, Phi = Phi, R = R, G = G,
                    Q = Q, N = dim(y)[2], t = t, K = K, log_lik = log_lik, m0 = mu0, 
                    R0 = Sigma0)
  
  ## KalmanFilter output ##
  mu_i_i1 = KF$m_i_i1
  mu_i_i = KF$m_i_i
  Sigma_i_i1 = KF$R_i_i1
  Sigma_i_i = KF$R_i_i
  innov_resid = KF$innov_resid
  innov_cov = KF$innov_cov
  ll = KF$ll
  
  ## Run Tangent KalmanFilter ##
  TKF = MultiTangentKalmanFilter(y = y, lin_pred = lin_pred, Phi = Phi, R = R, 
                                 G = G,Q = Q, Phi_grad = Phi_grad, 
                                 R_grad = R_grad,G_grad = G_grad, 
                                 Q_grad = Q_grad, n = n, t = t, K = K, 
                                 log_lik = log_lik, mu0 = mu0, 
                                 Sigma0 = Sigma0, mu0_grad = mu0_grad, 
                                 Sigma0_grad = Sigma0_grad, 
                                 n_param = n_param,
                                 lin_pred_grad = lin_pred_grad)
  
  ## TangentKalmanFilter output ##
  mu_i_i1_grad = TKF$mu_i_i1_grad
  mu_i_i_grad = TKF$mu_i_i_grad
  Sigma_i_i1_grad = TKF$Sigma_i_i1_grad
  Sigma_i_i_grad = TKF$Sigma_i_i_grad
  innov_resid_grad = TKF$innov_resid_grad
  innov_cov_grad = TKF$innov_cov_grad
  ll_grad = TKF$ll_grad 
  conditional_ll_grad = TKF$conditional_ll_grad
  
  
  ## Initialise grad2nd derivative arrays ##
  
  # Predictions
  mu_i_i1_grad2 <- array(0, c(t, K)) # Mean 2nd Derivative
  Sigma_i_i1_grad2 <- array(0, c(t, K, K)) # Covariance 2nd Derivative
  
  mu_i_i1_grad2 = rep(list(mu_i_i1_grad2),n_param)
  Sigma_i_i1_grad2 = rep(list(Sigma_i_i1_grad2),n_param)
  
  # Updates
  mu_i_i_grad2 <- array(0, c(t+1, K)) # Mean 2nd Derivative
  Sigma_i_i_grad2 <- array(0, c(t+1, K, K)) # Covariance 2nd Derivative
  
  mu_i_i_grad2 = rep(list(mu_i_i_grad2),n_param)
  Sigma_i_i_grad2 = rep(list(Sigma_i_i_grad2),n_param)
  
  # Innovations
  innov_resid_grad2 <- array(0, c(t, dim(y)[2])) # Residual 2nd Derivative
  innov_cov_grad2 <- array(0, c(t, dim(y)[2], dim(y)[2])) # Covariance 2nd Derivative
  
  innov_resid_grad2 = rep(list(innov_resid_grad2),n_param)
  innov_cov_grad2 = rep(list(innov_cov_grad2),n_param)
  
  # Contional log-likelihood
  conditional_ll_grad2 <- array(0, t)
  conditional_ll_grad2 = rep(list(conditional_ll_grad2),n_param)
  
  # Log-likelihood
  ll_grad2 = 0
  ll_grad2 = rep(list(ll_grad2),n_param)
  
  
  ## Initial gradient values ##
  for (j in grad_params){
    mu_i_i_grad2[[j]][1,] <- mu0_grad2[[j]] # Mean
    Sigma_i_i_grad2[[j]][1,,] <- Sigma0_grad2[[j]] # Covariance
  }
  
  ## Tangent Kalman Filter ##
  for (j in grad_params){
    for (i in 1:t) {
      
      ## Prediction
      
      ## -- Mean
      mu_i_i1_grad2[[j]][i,] <- 
        G_grad2[[j]] %*% mu_i_i[i,] + 
        2 * G_grad[[j]] %*% mu_i_i_grad[[j]][i,] +
        G %*% mu_i_i_grad2[[j]][i,]
      
      ## -- Covariance
      Sigma_i_i1_grad2[[j]][i,,] <- 
        G_grad2[[j]] %*% Sigma_i_i[i,,] %*% t(G) + 
        G %*% Sigma_i_i_grad2[[j]][i,,] %*% t(G) +
        G %*% Sigma_i_i[i,,] %*% t(G_grad2[[j]]) + 
        2 * G_grad[[j]] %*% Sigma_i_i_grad[[j]][i,,] %*% t(G) + 
        2 * G_grad[[j]] %*% Sigma_i_i[i,,] %*% t(G_grad[[j]]) + 
        2 * G %*% Sigma_i_i_grad[[j]][i,,] %*% t(G_grad[[j]]) + 
        Q_grad2[[j]]
      
      Sigma_i_i1_grad2[[j]][i,,] <- 
        (Sigma_i_i1_grad2[[j]][i,,] + t(Sigma_i_i1_grad2[[j]][i,,]))/2 
      
      ## Update
      
      ## -- 'm'
      m_i_grad2 <- 
        Phi_grad2[[j]] %*% mu_i_i1[i,] + 
        2 * Phi_grad[[j]] %*% mu_i_i1_grad[[j]][i,] +
        Phi %*% mu_i_i1_grad2[[j]][i,]
      
      ## -- Innovation Covariance
      S_i <- solve(innov_cov[i,,])
      
      S_i_grad2 <- 
        Phi_grad2[[j]] %*% Sigma_i_i1[i,,] %*% Phi_T +
        Phi %*% Sigma_i_i1_grad2[[j]][i,,] %*% Phi_T + 
        Phi %*% Sigma_i_i1[i,,] %*% Phi_grad2_T[[j]] +
        2 * Phi_grad[[j]] %*% Sigma_i_i1_grad[[j]][i,,] %*% Phi_T + 
        2 * Phi_grad[[j]] %*% Sigma_i_i1[i,,] %*% Phi_grad_T[[j]] + 
        2 * Phi %*% Sigma_i_i1_grad[[j]][i,,] %*% Phi_grad_T[[j]] + 
        R_grad2[[j]]
      
      S_i_grad2_inv <- 
        2 * innov_cov_grad[[j]][i,,] %*% S_i %*% innov_cov_grad[[j]][i,,] -
        innov_cov[i,,] %*% S_i_grad2 %*% innov_cov[i,,]
      
      innov_cov_grad2[[j]][i,,] <- S_i_grad2_inv
      
      ## -- Kalman Gain
      K_i <- Sigma_i_i1[i,,] %*% Phi_T %*% innov_cov[i,,]
      
      K_i_grad <- 
        Sigma_i_i1_grad[[j]][i,,] %*% Phi_T %*% innov_cov[i,,] +
        Sigma_i_i1[i,,] %*% Phi_grad_T[[j]] %*% innov_cov[i,,] + 
        Sigma_i_i1[i,,] %*% Phi_T %*% innov_cov_grad[[j]][i,,]
      
      K_i_grad2 <- 
        Sigma_i_i1_grad2[[j]][i,,] %*% Phi_T %*% innov_cov[i,,] +
        Sigma_i_i1[i,,] %*% Phi_grad2_T[[j]] %*% innov_cov[i,,] + 
        Sigma_i_i1[i,,] %*% Phi_T %*% innov_cov_grad2[[j]][i,,] +
        2 * Sigma_i_i1_grad[[j]][i,,] %*% Phi_grad_T[[j]] %*% innov_cov[i,,] + 
        2 * Sigma_i_i1_grad[[j]][i,,] %*% Phi_T %*% innov_cov_grad[[j]][i,,] + 
        2 * Sigma_i_i1[i,,] %*% Phi_grad_T[[j]] %*% innov_cov_grad[[j]][i,,]
      
      
      ## -- Innovation Residual
      innov_resid_grad2[[j]][i,] <- - m_i_grad2 - lin_pred_grad2[[j]][i,]
      
      ## -- Mean
      mu_i_i_grad2[[j]][i+1,] <- 
        mu_i_i1_grad2[[j]][i,] + 
        K_i_grad2 %*% innov_resid[i,] + 
        2 * K_i_grad %*% innov_resid_grad[[j]][i,] +
        K_i %*% innov_resid_grad2[[j]][i,]
      
      ## -- Covariance
      Sigma_i_i_grad2[[j]][i+1,,] <- 
        Sigma_i_i1_grad2[[j]][i,,] -
        K_i_grad2 %*% Phi %*% Sigma_i_i1[i,,] -
        K_i %*% Phi_grad2[[j]] %*% Sigma_i_i1[i,,] -
        K_i %*% Phi %*% Sigma_i_i1_grad2[[j]][i,,] -
        2 * K_i_grad %*% Phi_grad[[j]] %*% Sigma_i_i1[i,,] -
        2 * K_i_grad %*% Phi %*% Sigma_i_i1_grad[[j]][i,,] - 
        2 * K_i %*% Phi_grad[[j]] %*% Sigma_i_i1_grad[[j]][i,,]
      
      Sigma_i_i_grad2[[j]][i+1,,] <- 
        (Sigma_i_i_grad2[[j]][i+1,,] + t(Sigma_i_i_grad2[[j]][i+1,,]))/2
      
      ## Log-likelihood gradient
      if (log_lik){
        
        conditional_ll_grad2[[j]][i] = -
          0.5 * t(innov_resid_grad2[[j]][i,]) %*% innov_cov[i,,] %*% innov_resid[i,] -
          0.5 * t(innov_resid[i,]) %*% innov_cov_grad2[[j]][i,,] %*% innov_resid[i,] -
          0.5 * t(innov_resid[i,]) %*% innov_cov[i,,] %*% innov_resid_grad2[[j]][i,] -
          t(innov_resid_grad[[j]][i,]) %*% innov_cov_grad[[j]][i,,] %*% innov_resid[i,] -
          t(innov_resid_grad[[j]][i,]) %*% innov_cov[i,,] %*% innov_resid_grad[[j]][i,] -
          t(innov_resid[i,]) %*% innov_cov_grad[[j]][i,,] %*% innov_resid_grad[[j]][i,] -
          0.5 * sum(diag(innov_cov[i,,] %*% S_i_grad2)) +
          0.5 * sum(diag(innov_cov_grad[[j]][i,,] %*% S_i %*% innov_cov_grad[[j]][i,,] %*% S_i))
        
        ll_grad2[[j]] = ll_grad2[[j]] + conditional_ll_grad2[[j]][i]
      }
    }
  }
  
  ## Output
  output = list(mu_i_i = mu_i_i, Sigma_i_i = Sigma_i_i, mu_i_i1 = mu_i_i1, 
                Sigma_i_i1 = Sigma_i_i1, innov_resid = innov_resid, 
                innov_cov = innov_cov, K_i = K_i, mu_i_i_grad = mu_i_i_grad,
                Sigma_i_i_grad = Sigma_i_i_grad, mu_i_i1_grad = mu_i_i1_grad,
                Sigma_i_i1_grad = Sigma_i_i1_grad, 
                innov_resid_grad = innov_resid_grad, 
                innov_cov_grad = innov_cov_grad, K_i_grad = K_i_grad,
                mu_i_i_grad2 = mu_i_i_grad2,Sigma_i_i_grad2 = Sigma_i_i_grad2, 
                mu_i_i1_grad2 = mu_i_i1_grad2,
                Sigma_i_i1_grad2 = Sigma_i_i1_grad2, 
                innov_resid_grad2 = innov_resid_grad2, 
                innov_cov_grad2 = innov_cov_grad2, K_i_grad = K_i_grad2)
  
  if (log_lik){
    output = c(output, list(ll = ll, ll_grad = ll_grad,
                            conditional_ll_grad = conditional_ll_grad,
                            ll_grad2 = ll_grad2,
                            conditional_ll_grad2 = conditional_ll_grad2))
  }
  
  ## Print Output
  return(output)
}

#############################

######################################
# SPECTRAL 2ND TANGENT KALMAN FILTER #
######################################

TKF2_spectral <- function(w = NULL, w_FFT_real = NULL, param_obs = NULL, 
                          param_sig = NULL, n, t, K=n^2, log_lik = TRUE, 
                          dt = 1, m0 = rep(0,K), R0 = rep(0,K),
                          m0_grad = rep(0,K), R0_grad = rep(0,K), 
                          m0_grad2 = rep(0,K), R0_grad2 = rep(0,K), nu = 1,
                          grad_param, log_lik_grad = TRUE, 
                          log_lik_grad2 = TRUE){
  
  ## If 'w_FFT_real' not given, compute as inverse-FFT of w
  if(is.null(w_FFT_real)){
    w_FFT_real <- FFT_real_ST(STmatrix_to_vec(w),n=n,t=t,inv=TRUE)
  }
  
  ## If 'w_FFT_real' given as matrix,convert into vector
  else{
    if(class(w_FFT_real)=="matrix"){
      w_FFT_real <- STmatrix_to_vec(w_FFT_real)
    }
  }
  
  ## Select required terms from w_FFT_real
  if(K<n^2){
    basis_indices <- spde_initialise(n,t,K)$basis_indices
    w_FFT_indices <- rep(seq(0,n^2*(t-1),n^2),each=K)+basis_indices
    w_FFT_real = w_FFT_real[w_FFT_indices]
  }
  
  ## Prelims
  spde_ft = spde_initialise(n,t,K)
  wave = spde_ft$wave
  cosine_indices = spde_ft$cosine_indices
  n_cos = spde_ft$n_cos
  n_cos_sin = length(cosine_indices)
  K = spde_ft$K
  
  
  ## SPDE matrices
  G_vec <- G_f_vec(param_sig,wave,cosine_indices,dt=dt,n_cos)
  Q_vec <- Q_f_vec(param_sig,wave,n,n_cos,nu,dt=dt,norm=TRUE)
  
  G_1 = G_vec$G_1; G_2 = G_vec$G_2; G_cos = G_vec$G_cos
  Q_vec_cos = Q_vec[1:n_cos]; Q_vec_cos_sin = Q_vec[cosine_indices]
  tau2 <- param_obs$tau2_vals[[1]]
  
  
  ## SPDE matrix gradients
  R_grad_vec = R_grad_f_vec(param_obs,param_sig,n,n_obs=n^2)[[grad_param]]
  G_grad_vec = G_grad_f_vec(param_obs,param_sig,n,K,wave,cosine_indices,dt=dt,nu,n_cos)[[grad_param]]
  Q_grad_vec = Q_grad_f_vec(param_obs,param_sig,n,K,wave,n_cos,nu,dt=dt,norm=TRUE)[[grad_param]]
  
  G_grad_cos = G_grad_vec$G_grad_cos
  G_grad_1 = G_grad_vec$G_grad_1
  G_grad_2 = G_grad_vec$G_grad_2
  
  ## SPDE matrix 2nd derivatives
  G_grad2_vec = G_grad2_f_vec(param_obs,param_sig,n,K,wave,cosine_indices,dt=dt,nu,n_cos)[[grad_param]]
  Q_grad2_vec = Q_grad2_f_vec(param_obs,param_sig,n,K,wave,n_cos,nu,dt=dt,norm=TRUE)[[grad_param]]
  
  G_grad2_cos = G_grad2_vec$G_grad2_cos
  G_grad2_1 = G_grad2_vec$G_grad2_1
  G_grad2_2 = G_grad2_vec$G_grad2_2
  
  
  ## Kalman Filter 
  KF <- kf_spectral(w_FFT_real = w_FFT_real,m_i_i1 = rep(0,K*t),
                    m_i_i=c(m0,rep(0,K*t)),R_i_i1=rep(0,K*t),
                    R_i_i=c(R0,rep(0,K*t)),G_cos = G_cos, 
                    G_1 = G_1,G_2 = G_2,spectrum = Q_vec, 
                    t = t, n_cos = n_cos,  n_cos_sin = n_cos_sin, 
                    tau2 = tau2)
  
  
  
  ## Tangent Kalman Filter
  TKF = tkf_spectral(w_FFT_real,KF$m_i_i1,KF$m_i_i,KF$R_i_i1,KF$R_i_i,
                     m_i_i1_grad=rep(0,K*t),
                     m_i_i_grad=c(m0_grad,rep(0,K*t)),
                     R_i_i1_grad=rep(0,K*t),
                     R_i_i_grad=c(R0_grad,rep(0,K*t)),
                     G_cos,G_1,G_2,G_grad_cos,G_grad_1,G_grad_2,Q_vec,
                     Q_grad_vec,t, n_cos, n_cos_sin, tau2, R_grad_vec)
  
  
  ## Second Order Tangent Kalman Filter
  TKF2 = tkf2_spectral(w_FFT_real, KF$m_i_i1, KF$m_i_i, KF$R_i_i1, KF$R_i_i,
                      TKF$mu_i_i1_grad, TKF$mu_i_i_grad, TKF$Sigma_i_i1_grad,
                      TKF$Sigma_i_i_grad, m_i_i1_grad2 = rep(0,K*t),
                      m_i_i_grad2 = c(m0_grad2,rep(0,K*t)),
                      R_i_i1_grad2 = rep(0,K*t),
                      R_i_i_grad2 = c(R0_grad2,rep(0,K*t)), G_cos, G_1, G_2, 
                      G_grad_cos, G_grad_1, G_grad_2, G_grad2_cos, G_grad2_1,
                      G_grad2_2, Q_vec, Q_grad_vec, Q_grad2_vec, t, n_cos, 
                      n_cos_sin, tau2, R_grad_vec)
  
  ## Prepare output
  mu_i_i = vec_to_STmatrix(KF$m_i_i,t=t+1)
  mu_i_i1 = vec_to_STmatrix(KF$m_i_i1,t=t)
  Sigma_i_i = vec_to_STmatrix(KF$R_i_i,t=t+1)
  Sigma_i_i1 = vec_to_STmatrix(KF$R_i_i1,t=t)
  
  mu_i_i_grad = vec_to_STmatrix(TKF$mu_i_i_grad,t=t+1)
  mu_i_i1_grad = vec_to_STmatrix(TKF$mu_i_i1_grad,t=t)
  Sigma_i_i_grad = vec_to_STmatrix(TKF$Sigma_i_i_grad,t=t+1)
  Sigma_i_i1_grad = vec_to_STmatrix(TKF$Sigma_i_i1_grad,t=t)
  
  mu_i_i_grad2 = vec_to_STmatrix(TKF2$mu_i_i_grad2,t=t+1)
  mu_i_i1_grad2 = vec_to_STmatrix(TKF2$mu_i_i1_grad2,t=t)
  Sigma_i_i_grad2 = vec_to_STmatrix(TKF2$Sigma_i_i_grad2,t=t+1)
  Sigma_i_i1_grad2 = vec_to_STmatrix(TKF2$Sigma_i_i1_grad2,t=t)
  
  ## Output
  output <- list(mu_i_i = mu_i_i, Sigma_i_i = Sigma_i_i, 
                 mu_i_i1 = mu_i_i1, Sigma_i_i1 = Sigma_i_i1,
                 mu_i_i_grad = mu_i_i_grad, Sigma_i_i_grad = Sigma_i_i_grad, 
                 mu_i_i1_grad = mu_i_i1_grad, Sigma_i_i1_grad = Sigma_i_i1_grad,
                 mu_i_i_grad2 = mu_i_i_grad2, Sigma_i_i_grad2 = Sigma_i_i_grad2, 
                 mu_i_i1_grad2 = mu_i_i1_grad2, Sigma_i_i1_grad2 = Sigma_i_i1_grad2)
  
  ## If 'log_lik'=TRUE: compute likelihood  & amend output
  if(log_lik){
    
    if(K==n^2){
      ll = ll_spectral(log_lik=0,w_FFT_real, KF$m_i_i1, KF$R_i_i1, t, K, 
                       tau2)$ll
    }
    
    if(K<n^2){
      ll = ll_non_spectral(log_lik=0,w,KF$m_i_i1,KF$R_i_i1,t,n,K,tau2)$ll
    }
    
    output = c(output, list(ll=ll))
  }
  
  ## If 'log_lik_grad'=TRUE: compute grad log-likelihood & amend output
  if(log_lik_grad){
    if(K==n^2){
      loglik_grad = ll_grad_spectral(log_lik_grad=0,
                                     conditional_log_lik_grad = rep(0,t),
                                     w_FFT_real,KF$m_i_i1,KF$R_i_i1,
                                     TKF$mu_i_i1_grad,TKF$Sigma_i_i1_grad,t,K,
                                     tau2,R_grad_vec)
      
      ll_grad = loglik_grad$ll_grad
      conditional_ll_grad = loglik_grad$conditional_ll_grad
    }
    
    if(K<n^2){
      loglik_grad = ll_grad_non_spectral(log_lik_grad=0,
                                         conditional_log_lik_grad=rep(0,t),
                                         w=w,m_i_i1=KF$m_i_i1,R_i_i1=KF$R_i_i1,
                                         m_i_i1_grad=TKF$mu_i_i1_grad,
                                         R_i_i1_grad=TKF$Sigma_i_i1_grad,
                                         t=t,n=n,K=K,tau2=tau2,
                                         R_grad = diag(R_grad_vec))
      ll_grad = loglik_grad$ll_grad
      conditional_ll_grad = loglik_grad$conditional_ll_grad
    }
  }
  
  ## If 'log_lik_grad2'=TRUE: compute 2nd derivative of log-lik & amend output
  if(log_lik_grad2){
    if(K==n^2){
      loglik_grad2 = ll_grad2_spectral(log_lik_grad2=0,
                                       conditional_log_lik_grad2 = rep(0,t),
                                       w_FFT_real,KF$m_i_i1,KF$R_i_i1,
                                       TKF$mu_i_i1_grad,TKF$Sigma_i_i1_grad,
                                       TKF2$mu_i_i1_grad2,TKF2$Sigma_i_i1_grad2,
                                       t,K,tau2,R_grad_vec)
      ll_grad2 = loglik_grad2$ll_grad2
      conditional_ll_grad2 = loglik_grad2$conditional_ll_grad2
    }
    
    if(K<n^2){
      loglik_grad2 = ll_grad2_non_spectral(log_lik_grad2=0,
                                           conditional_log_lik_grad2 = rep(0,t),
                                           w,KF$m_i_i1,KF$R_i_i1,
                                           TKF$mu_i_i1_grad,TKF$Sigma_i_i1_grad,
                                           TKF2$mu_i_i1_grad2,TKF2$Sigma_i_i1_grad2,
                                           t,n,K,tau2,diag(R_grad_vec))
      ll_grad2 = loglik_grad2$ll_grad2
      conditional_ll_grad2 = loglik_grad2$conditional_ll_grad2
    }
    
    output = c(output,list(ll_grad2=ll_grad2,
                           conditional_ll_grad2 = conditional_ll_grad2))
  }
  
  ## Return output
  return(output)
}

Multi_TKF2_spectral <- function(w = NULL, w_FFT_real = NULL, param_obs = NULL,
                                param_sig = NULL, n, t, K=n^2, log_lik = TRUE, 
                                dt = 1, m0 = rep(0,K), R0 = rep(0,K),
                                m0_grad = rep(list(rep(0,K)),9), 
                                R0_grad = rep(list(rep(0,K)),9), 
                                m0_grad2 = rep(list(rep(0,K)),9),
                                R0_grad2 = rep(list(rep(0,K)),9),
                                nu = 1, grad_params, log_lik_grad = TRUE,
                                log_lik_grad2 = TRUE){
  
  ## If 'w_FFT_real' not given, compute as inverse-FFT of w
  if(is.null(w_FFT_real)){
    w_FFT_real <- FFT_real_ST(STmatrix_to_vec(w),n=n,t=t,inv=TRUE)
  }
  
  ## If 'w_FFT_real' given as matrix,convert into vector
  else{
    if(class(w_FFT_real)=="matrix"){
      w_FFT_real <- STmatrix_to_vec(w_FFT_real)
    }
  }
  
  ## Select required terms from w_FFT_real
  if(K<n^2){
    basis_indices <- spde_initialise(n,t,K)$basis_indices
    w_FFT_indices <- rep(seq(0,n^2*(t-1),n^2),each=K)+basis_indices
    w_FFT_real = w_FFT_real[w_FFT_indices]
  }
  
  ## Prelims
  spde_ft = spde_initialise(n,t,K)
  wave = spde_ft$wave
  cosine_indices = spde_ft$cosine_indices
  n_cos = spde_ft$n_cos
  n_cos_sin = length(cosine_indices)
  K = spde_ft$K
  
  
  ## SPDE matrices
  G_vec <- G_f_vec(param_sig,wave,cosine_indices,dt=dt,n_cos)
  Q_vec <- Q_f_vec(param_sig,wave,n,n_cos,nu,dt=dt,norm=TRUE)
  
  G_1 = G_vec$G_1; G_2 = G_vec$G_2; G_cos = G_vec$G_cos
  Q_vec_cos = Q_vec[1:n_cos]; Q_vec_cos_sin = Q_vec[cosine_indices]
  tau2 <- param_obs$tau2_vals[[1]]
  
  ## Kalman Filter 
  KF <- kf_spectral(w_FFT_real = w_FFT_real,m_i_i1 = rep(0,K*t),
                    m_i_i=c(m0,rep(0,K*t)),R_i_i1=rep(0,K*t),
                    R_i_i=c(R0,rep(0,K*t)),G_cos = G_cos, 
                    G_1 = G_1,G_2 = G_2,spectrum = Q_vec, 
                    t = t, n_cos = n_cos,  n_cos_sin = n_cos_sin, 
                    tau2 = tau2)
  
  ## Kalman Filter Outputs
  mu_i_i = vec_to_STmatrix(KF$m_i_i,t=t+1)
  mu_i_i1 = vec_to_STmatrix(KF$m_i_i1,t=t)
  Sigma_i_i = vec_to_STmatrix(KF$R_i_i,t=t+1)
  Sigma_i_i1 = vec_to_STmatrix(KF$R_i_i1,t=t)
  
  output = list(mu_i_i = mu_i_i, Sigma_i_i = Sigma_i_i,mu_i_i1 = mu_i_i1, 
                Sigma_i_i1 = Sigma_i_i1)
  
  
  ## Initialise gradient arrays
  
  # Predictions
  mu_i_i1_grad = rep(list(array(0, c(t, K))),9) # Mean gradients
  Sigma_i_i1_grad = rep(list(array(0, c(t, K))),9) # Cov. gradients
  
  # Updates
  mu_i_i_grad = rep(list(array(0, c(t+1, K))),9) # Mean gradients
  Sigma_i_i_grad = rep(list(array(0, c(t+1, K))),9) #Cov. gradients
  
  # Log-likelihood
  ll_grad = rep(list(0),9)
  
  # Conditional log-likelihood
  conditional_ll_grad = rep(list(rep(0,t)),9)
  
  
  ## SPDE matrix gradients
  R_grad_vec = R_grad_f_vec(param_obs,param_sig,n,n_obs=n^2)
  G_grad_vec = G_grad_f_vec(param_obs,param_sig,n,K,wave,cosine_indices,dt=dt,nu,n_cos)
  Q_grad_vec = Q_grad_f_vec(param_obs,param_sig,n,K,wave,n_cos,nu,dt=dt,norm=TRUE)
  
  G_grad_cos = lapply(1:9, function(i) G_grad_vec[[i]]$G_grad_cos)
  G_grad_1 = lapply(1:9, function(i) G_grad_vec[[i]]$G_grad_1)
  G_grad_2 = lapply(1:9, function(i) G_grad_vec[[i]]$G_grad_2)
  
  
  ## Initialise 2nd derivative arrays
  
  # Predictions
  mu_i_i1_grad2 = rep(list(array(0, c(t, K))),9) # Mean 2nd deriv
  Sigma_i_i1_grad2 = rep(list(array(0, c(t, K))),9) # Cov. 2nd deriv
  
  # Updates
  mu_i_i_grad2 = rep(list(array(0, c(t+1, K))),9) # Mean 2nd deriv
  Sigma_i_i_grad2 = rep(list(array(0, c(t+1, K))),9) #C ov. 2nd deriv
  
  # Log-likelihood
  ll_grad2 = rep(list(0),9)
  conditional_ll_grad2 = rep(list(rep(0,t)),9)
  
  
  ## SPDE matrix 2nd derivatives
  G_grad2_vec = G_grad2_f_vec(param_obs,param_sig,n,K,wave,cosine_indices,dt=dt,nu,n_cos)
  Q_grad2_vec = Q_grad2_f_vec(param_obs,param_sig,n,K,wave,n_cos,nu,dt=dt,norm=TRUE)
  
  G_grad2_cos = lapply(1:9, function(i) G_grad2_vec[[i]]$G_grad2_cos)
  G_grad2_1 = lapply(1:9, function(i) G_grad2_vec[[i]]$G_grad2_1)
  G_grad2_2 = lapply(1:9, function(i) G_grad2_vec[[i]]$G_grad2_2)
  
  
  for (j in grad_params){
    
    ## Tangent Kalman Filter
    TKF = tkf_spectral(w_FFT_real,KF$m_i_i1,KF$m_i_i,KF$R_i_i1,KF$R_i_i,
                       m_i_i1_grad=rep(0,K*t),
                       m_i_i_grad=c(m0_grad[[j]],rep(0,K*t)),
                       R_i_i1_grad=rep(0,K*t),
                       R_i_i_grad=c(R0_grad[[j]],rep(0,K*t)),
                       G_cos,G_1,G_2,G_grad_cos[[j]],G_grad_1[[j]],
                       G_grad_2[[j]],Q_vec,Q_grad_vec[[j]],t, n_cos, n_cos_sin, 
                       tau2, R_grad_vec[[j]])
    
    ## Tangent Kalman Filter Ouputs
    mu_i_i_grad[[j]] = vec_to_STmatrix(TKF$mu_i_i_grad,t=t+1)
    mu_i_i1_grad[[j]] = vec_to_STmatrix(TKF$mu_i_i1_grad,t=t)
    Sigma_i_i_grad[[j]] = vec_to_STmatrix(TKF$Sigma_i_i_grad,t=t+1)
    Sigma_i_i1_grad[[j]] = vec_to_STmatrix(TKF$Sigma_i_i1_grad,t=t)
    
    
    ## Log-likelihood gradient
    if (log_lik_grad){
      if(K==n^2){
        loglik_grad = ll_grad_spectral(log_lik_grad=0,
                                       conditional_log_lik_grad = rep(0,t),
                                       w_FFT_real,KF$m_i_i1,KF$R_i_i1,
                                       TKF$mu_i_i1_grad,TKF$Sigma_i_i1_grad,t,K,
                                       tau2,R_grad_vec[[j]])
        
        ll_grad[[j]] = loglik_grad$ll_grad
        conditional_ll_grad[[j]] = loglik_grad$conditional_ll_grad
      }
      
      if(K<n^2){
        loglik_grad = ll_grad_non_spectral(log_lik_grad=0,
                                           conditional_log_lik_grad=rep(0,t),
                                           w=w,m_i_i1=KF$m_i_i1,R_i_i1=KF$R_i_i1,
                                           m_i_i1_grad=TKF$mu_i_i1_grad,
                                           R_i_i1_grad=TKF$Sigma_i_i1_grad,
                                           t=t,n=n,K=K,tau2=tau2,
                                           R_grad = diag(R_grad_vec[[j]]))
        ll_grad[[j]] = loglik_grad$ll_grad
        conditional_ll_grad[[j]] = loglik_grad$conditional_ll_grad
      }
    }
    
    
    ## Second Order Tangent Kalman Filter
    TKF2 = tkf2_spectral(w_FFT_real,KF$m_i_i1,KF$m_i_i,KF$R_i_i1,KF$R_i_i,
                         TKF$mu_i_i1_grad,TKF$mu_i_i_grad,TKF$Sigma_i_i1_grad,
                         TKF$Sigma_i_i_grad,
                         m_i_i1_grad2=rep(0,K*t),
                         m_i_i_grad2=c(m0_grad2[[j]],rep(0,K*t)),
                         R_i_i1_grad2=rep(0,K*t),
                         R_i_i_grad2=c(R0_grad2[[j]],rep(0,K*t)),
                         G_cos,G_1,G_2,G_grad_cos[[j]],G_grad_1[[j]],
                         G_grad_2[[j]],G_grad2_cos[[j]],G_grad2_1[[j]],
                         G_grad2_2[[j]],Q_vec,Q_grad_vec[[j]],Q_grad2_vec[[j]],
                         t, n_cos, n_cos_sin, tau2, R_grad_vec[[j]])
    
    ## Second Order Tangent Kalman Filter Outputs
    mu_i_i_grad2[[j]] = vec_to_STmatrix(TKF2$mu_i_i_grad2,t=t+1)
    mu_i_i1_grad2[[j]] = vec_to_STmatrix(TKF2$mu_i_i1_grad2,t=t)
    Sigma_i_i_grad2[[j]] = vec_to_STmatrix(TKF2$Sigma_i_i_grad2,t=t+1)
    Sigma_i_i1_grad2[[j]] = vec_to_STmatrix(TKF2$Sigma_i_i1_grad2,t=t)
    
    ## Log-likelihood 2nd Derivative
    if (log_lik_grad2){
      if(K==n^2){
        loglik_grad2 = ll_grad2_spectral(log_lik_grad2=0,
                                        conditional_log_lik_grad2 = rep(0,t),
                                        w_FFT_real,KF$m_i_i1,
                                        KF$R_i_i1,TKF$mu_i_i1_grad,
                                        TKF$Sigma_i_i1_grad,TKF2$mu_i_i1_grad2,
                                        TKF2$Sigma_i_i1_grad2,t,K,tau2,
                                        R_grad_vec[[j]])
        
        ll_grad2[[j]] = loglik_grad2$ll_grad2
        conditional_ll_grad2[[j]] = loglik_grad2$conditional_ll_grad2
      }
      if(K<n^2){
        loglik_grad2 = ll_grad2_non_spectral(log_lik_grad2=0,
                                             conditional_log_lik_grad2 = rep(0,t),
                                             w,KF$m_i_i1, KF$R_i_i1,
                                             TKF$mu_i_i1_grad,
                                             TKF$Sigma_i_i1_grad,
                                             TKF2$mu_i_i1_grad2,
                                             TKF2$Sigma_i_i1_grad2,t,n,K,tau2,
                                             diag(R_grad_vec[[j]]))
        
        ll_grad2[[j]] = loglik_grad2$ll_grad2
        conditional_ll_grad2[[j]] = loglik_grad2$conditional_ll_grad2
      }
    }
  }
  
  ## Output
  output <- c(output, (list(mu_i_i_grad = mu_i_i_grad, 
                            Sigma_i_i_grad = Sigma_i_i_grad, 
                            mu_i_i1_grad = mu_i_i1_grad, 
                            Sigma_i_i1_grad = Sigma_i_i1_grad)))
  
  output <- c(output,(list(mu_i_i_grad2 = mu_i_i_grad2, 
                           Sigma_i_i_grad2 = Sigma_i_i_grad2, 
                           mu_i_i1_grad2 = mu_i_i1_grad2, 
                           Sigma_i_i1_grad2 = Sigma_i_i1_grad2)))
  
  ## If 'log_lik'=TRUE: compute likelihood  & amend output
  if(log_lik){
    ll = ll_spectral(log_lik=0,w_FFT_real, KF$m_i_i1, KF$R_i_i1, t, K, 
                     tau2)$ll
    
    output = c(output, list(ll=ll))
  }
  
  ## If 'log_lik_grad'=TRUE: amend output
  if(log_lik_grad){
    output = c(output,list(ll_grad=ll_grad,
                           conditional_ll_grad = conditional_ll_grad))
  }
  
  ## If 'log_lik_grad2'=TRUE: amend output
  if(log_lik_grad2){
    output = c(output,list(ll_grad2=ll_grad2,
                           conditional_ll_grad2 = conditional_ll_grad2))
  }
  
  ## Return output
  return(output)
}

tkf2_spectral = function(w_FFT_real,m_i_i1,m_i_i,R_i_i1,R_i_i,
                         m_i_i1_grad,m_i_i_grad,R_i_i1_grad,R_i_i_grad,
                         m_i_i1_grad2,m_i_i_grad2,R_i_i1_grad2,R_i_i_grad2,
                         G_cos,G_1,G_2,G_grad_cos,G_grad_1,G_grad_2,
                         G_grad2_cos,G_grad2_1,G_grad2_2, Q_vec, Q_grad_vec, 
                         Q_grad2_vec, t, n_cos, n_cos_sin, tau2, R_grad){
  
  tkf2 = .C("tkf2_spectral", y_FFT = as.double(w_FFT_real),
            mu_i_i1=as.double(m_i_i1),mu_i_i=as.double(m_i_i),
            Sigma_i_i1=as.double(R_i_i1),Sigma_i_i = as.double(R_i_i), 
            mu_i_i1_grad = as.double(m_i_i1_grad),
            mu_i_i_grad = as.double(m_i_i_grad),
            Sigma_i_i1_grad = as.double(R_i_i1_grad), 
            Sigma_i_i_grad = as.double(R_i_i_grad), 
            mu_i_i1_grad2 = as.double(m_i_i1_grad2),
            mu_i_i_grad2 = as.double(m_i_i_grad2),
            Sigma_i_i1_grad2 = as.double(R_i_i1_grad2), 
            Sigma_i_i_grad2 = as.double(R_i_i_grad2),
            G_cos = as.double(G_cos), G_1 = as.double(G_1), 
            G_2 = as.double(G_2), G_grad_cos = as.double(G_grad_cos), 
            G_grad_1 = as.double(G_grad_1), G_grad_2 = as.double(G_grad_2), 
            G_grad2_cos = as.double(G_grad2_cos), 
            G_grad2_1 = as.double(G_grad2_1), G_grad_2 = as.double(G_grad2_2),
            Q = as.double(Q_vec), Q_grad = as.double(Q_grad_vec),
            Q_grad2 = as.double(Q_grad2_vec), T = as.integer(t), 
            n_cos = as.integer(n_cos), n_cos_sin=as.integer(n_cos_sin), 
            tau2 = as.double(tau2), R_grad = as.double(R_grad))
  
  
  output = list(mu_i_i1_grad2 = tkf2$mu_i_i1_grad2, 
                mu_i_i_grad2 = tkf2$mu_i_i_grad2, 
                Sigma_i_i1_grad2 = tkf2$Sigma_i_i1_grad2,
                Sigma_i_i_grad2 = tkf2$Sigma_i_i_grad2)
  
  return(output)
}

######################################

## ---

#######################
# GRAD LOG-LIKELIHOOD #
#######################

log_likelihood_grad <- function(param_obs = NULL, param_sig = NULL,
                                w = NULL, w_FFT_real = NULL, n, t,
                                K=n^2, dt = 1, m0 = rep(0,K), R0 = rep(0,K), 
                                m0_grad = rep(list(rep(0,K)),9), 
                                R0_grad = rep(list(rep(0,K)),9),
                                nu = 1, grad_params=1:9, x=NULL, log_scale=FALSE, 
                                log_indices = c(1,2,3,4,5,9), negative=FALSE){
  
  ## Parameters ##
  param = c(param_sig,as.numeric(unlist(param_obs$tau2_vals)))
  
  ## If parameters provided on log-scale, rescale
  if(log_scale){
    param[log_indices] <- exp(param[log_indices])
  }
  
  ## Signal and observation parameters ##
  param_sig = param[1:8]
  param_obs = param_obs_func(m=1,m_indices=list(1:n^2),
                             tau2_vals=as.numeric(list(param[9])))
  
  ## Compute linear predictor X*beta
  if(!is.null(x)){
    linear_pred <- apply(x,c(2,3),linear_predictor,beta = param[-c(1:9)])
  }
  else{
    linear_pred <- 0
  }
  
  ## If 'w_FFT_real' not given, compute as inverse-FFT of (w-x*beta)
  if(is.null(w_FFT_real)){
    w_FFT_real <- FFT_real_ST(STmatrix_to_vec(w-linear_pred),n=n,t=t,inv=TRUE)
  }
  
  ## Compute log-likelihood using 'Multi_TKF_spectral'
  ll_grad <- unlist(Multi_TKF_spectral(w=w,w_FFT_real=w_FFT_real,
                                       param_obs=param_obs,param_sig=param_sig,
                                       n=n,t=t,K=K,log_lik=TRUE,dt=dt,m0=m0,
                                       R0=R0,m0_grad=m0_grad,R0_grad=R0_grad,
                                       nu=nu,grad_params=grad_params,
                                       log_lik_grad=TRUE)$ll_grad)
  
  ## Output
  if(negative){
    return(-ll_grad)
  }
  else{
    return(ll_grad)
  }
}

ll_grad_non_spectral = function(log_lik_grad = 0,
                                conditional_log_lik_grad = rep(0,t),w,m_i_i1,
                                R_i_i1,m_i_i1_grad,R_i_i1_grad,t,n,K,tau2,
                                R_grad){
  
  ## SPDE_FT object
  SPDE_FT = spde_initialise(n,t,K)
  wave = SPDE_FT$wave
  cosine_indices = SPDE_FT$cosine_indices
  n_cos = SPDE_FT$n_cos
  K = SPDE_FT$K
  
  ## FFT matrix
  Phi = FFT_real_matrix(wave,cosine_indices,n_cos,n)
  Phi_T = t(Phi)
  
  ## R matrix
  R = diag(rep(tau2,n^2))
  
  ## FFT matrix grad
  Phi_grad = matrix(0,nrow=n^2,ncol=K)
  Phi_grad_T = t(Phi_grad)
  
  ## Initialise
  ll_grad = log_lik_grad
  conditional_ll_grad = conditional_log_lik_grad
  
  for (i in 1:t){
    S_i = R + Phi%*%diag(R_i_i1[(i-1)*K + 1:K])%*%Phi_T
    S_i_inv = solve(S_i)
    resid = w[i,] - Phi%*%m_i_i1[(i-1)*K + 1:K]
    
    S_i_grad = 
      Phi_grad %*% diag(R_i_i1[(i-1)*K + 1:K]) %*% Phi_T +
      Phi %*% diag(R_i_i1_grad[(i-1)*K + 1:K]) %*% Phi_T + 
      Phi %*% diag(R_i_i1[(i-1)*K + 1:K]) %*% Phi_grad_T +
      R_grad
    S_i_grad_inv <- -S_i_inv %*% S_i_grad %*% S_i_inv
    
    resid_grad = - Phi_grad %*% m_i_i1[(i-1)*K + 1:K] - 
      Phi %*% m_i_i1_grad[(i-1)*K + 1:K]
    
    conditional_ll_grad[i] = -
      0.5 * t(resid_grad) %*% S_i_inv %*% resid -
      0.5 * t(resid)  %*% S_i_grad_inv %*% resid -
      0.5 * t(resid) %*% S_i_inv %*% resid_grad -
      0.5 * sum(diag(S_i_inv %*% S_i_grad))
    
    ll_grad = ll_grad + conditional_ll_grad[i]
  }
  
  return(list(ll_grad=ll_grad,
              conditional_ll_grad = conditional_ll_grad))
}

ll_grad_spectral = function(log_lik_grad=0,conditional_log_lik_grad = rep(0,t),
                            w_FFT_real,m_i_i1,R_i_i1,m_i_i1_grad,R_i_i1_grad,
                            t,K,tau2,R_grad){
  
  log_lik_grad = .C("ll_grad_spectral",ll_grad=as.double(log_lik_grad), 
                    conditional_ll_grad = as.double(conditional_log_lik_grad),
                    y_FFT=as.double(w_FFT_real), 
                    mu_i_i1=as.double(m_i_i1),
                    Sigma_i_i1=as.double(R_i_i1), 
                    mu_i_i1_grad=as.double(m_i_i1_grad),
                    Sigma_i_i1_grad=as.double(R_i_i1_grad), 
                    T=as.integer(t),K = as.integer(K), 
                    tau2 = as.double(tau2), 
                    R_grad = as.double(R_grad))
  
  output = list(ll_grad = log_lik_grad$ll_grad,
                conditional_ll_grad = log_lik_grad$conditional_ll_grad)
  return(output)
}

#######################

## ---

###########################
# 2nd GRAD LOG-LIKELIHOOD #
###########################

log_likelihood_grad2 <- function(param_obs = NULL, param_sig = NULL,
                                 w = NULL, w_FFT_real = NULL, n, t,
                                 K=n^2, dt = 1, m0 = rep(0,K), R0 = rep(0,K), 
                                 m0_grad = rep(list(rep(0,K)),9), 
                                 R0_grad = rep(list(rep(0,K)),9),
                                 m0_grad2 = rep(list(rep(0,K)),9), 
                                 R0_grad2 = rep(list(rep(0,K)),9),
                                 nu = 1, grad_params = 1:9, x = NULL,
                                 log_scale = FALSE, 
                                 log_indices=c(1,2,3,4,5,9),negative=FALSE){
  
  ## Parameters ##
  param = c(param_sig,as.numeric(unlist(param_obs$tau2_vals)))
  
  ## If parameters provided on log-scale, rescale
  if(log_scale){
    param[log_indices] <- exp(param[log_indices])
  }
  
  ## Signal and observation parameters ##
  param_sig = param[1:8]
  param_obs = param_obs_func(1,list(1:n^2),list(param[9]))
  
  
  ## Compute linear predictor X*beta
  if(!is.null(x)){
    linear_pred <- apply(x,c(2,3),linear_predictor,beta = param[-c(1:9)])
  }
  else{
    linear_pred <- 0
  }
  
  ## If 'w_FFT_real' not given, compute as inverse-FFT of (w-x*beta)
  if(is.null(w_FFT_real)){
    w_FFT_real <- FFT_real_ST(STmatrix_to_vec(w-linear_pred),n=n,t=t,inv=TRUE)
  }
  
  ## Compute log-likelihood using 'Multi_TKF_spectral'
  ll_grad2 <- unlist(Multi_TKF2_spectral(w=w,w_FFT_real=w_FFT_real,
                                         param_obs = param_obs,
                                         param_sig = param_sig,
                                         n=n,t=t,K=K,log_lik=TRUE,dt=dt,m0=m0,
                                         R0=R0,m0_grad=m0_grad,R0_grad=R0_grad,
                                         m0_grad2=m0_grad2,R0_grad2=R0_grad2,
                                         nu=nu,grad_params=grad_params,
                                         log_lik_grad=TRUE,
                                         log_lik_grad2=TRUE)$ll_grad2)
  
  ## Output
  if(negative){
    return(-ll_grad2)
  }
  else{
    return(ll_grad2)
  }
}

ll_grad2_non_spectral = function(log_lik_grad2=0,
                                 conditional_log_lik_grad2 = rep(0,t),w,m_i_i1,
                                 R_i_i1,m_i_i1_grad,R_i_i1_grad,m_i_i1_grad2,
                                 R_i_i1_grad2,t,n,K,tau2,R_grad){
  
  ## SPDE_FT object
  SPDE_FT = spde_initialise(n,t,K)
  wave = SPDE_FT$wave
  cosine_indices = SPDE_FT$cosine_indices
  n_cos = SPDE_FT$n_cos
  K = SPDE_FT$K
  
  ## FFT matrix
  Phi = FFT_real_matrix(wave,cosine_indices,n_cos,n)
  Phi_T = t(Phi)
  
  ## R matrix
  R = diag(rep(tau2,n^2))
  
  ## FFT matrix grad
  Phi_grad = matrix(0,nrow=n^2,ncol=K)
  Phi_grad_T = t(Phi_grad)
  
  ## FFT matrix grad2
  Phi_grad2 = matrix(0,nrow=n^2,ncol=K)
  Phi_grad2_T = t(Phi_grad2)
  
  ## R matrix grad2
  R_grad2 = diag(rep(0,n^2))
  
  ## Initialise
  ll_grad2 = log_lik_grad2
  conditional_ll_grad2 = conditional_log_lik_grad2
  
  for (i in 1:t){
    S_i = R + Phi%*%diag(R_i_i1[(i-1)*K + 1:K])%*%Phi_T
    S_i_inv = solve(S_i)
    resid = w[i,] - Phi%*%m_i_i1[(i-1)*K + 1:K]
    
    S_i_grad = 
      Phi_grad %*% diag(R_i_i1[(i-1)*K + 1:K]) %*% Phi_T +
      Phi %*% diag(R_i_i1_grad[(i-1)*K + 1:K]) %*% Phi_T + 
      Phi %*% diag(R_i_i1[(i-1)*K + 1:K]) %*% Phi_grad_T +
      R_grad
    S_i_grad_inv <- -S_i_inv %*% S_i_grad %*% S_i_inv
    
    resid_grad = - Phi_grad %*% m_i_i1[(i-1)*K + 1:K] - 
      Phi %*% m_i_i1_grad[(i-1)*K + 1:K]
    
    S_i_grad2 = 
      Phi_grad2 %*% diag(R_i_i1[(i-1)*K + 1:K]) %*% Phi_T +
      Phi %*% diag(R_i_i1_grad2[(i-1)*K + 1:K]) %*% Phi_T + 
      Phi %*% diag(R_i_i1[(i-1)*K + 1:K]) %*% Phi_grad2_T +
      2 * Phi_grad %*% diag(R_i_i1_grad[(i-1)*K + 1:K]) %*% Phi_T + 
      2 * Phi_grad %*% diag(R_i_i1[(i-1)*K + 1:K]) %*% Phi_grad_T + 
      2 * Phi %*% diag(R_i_i1_grad[(i-1)*K + 1:K]) %*% Phi_grad_T + 
      R_grad2
    
    S_i_grad2_inv = 
      2 * S_i_grad_inv %*% S_i %*% S_i_grad_inv -
      S_i_inv %*% S_i_grad2 %*% S_i_inv
    
    resid_grad2 = - Phi_grad2 %*% m_i_i1[(i-1)*K + 1:K] - 
      2 * Phi_grad %*% m_i_i1_grad[(i-1)*K + 1:K] -
      Phi %*% m_i_i1_grad2[(i-1)*K + 1:K]
    
    conditional_ll_grad2[i] = -
      0.5 * t(resid_grad2) %*% S_i_inv %*% resid -
      0.5 * t(resid) %*% S_i_grad2_inv %*% resid -
      0.5 * t(resid) %*% S_i_inv %*% resid_grad2 -
      t(resid_grad) %*% S_i_grad_inv %*% resid -
      t(resid_grad) %*% S_i_inv %*% resid_grad -
      t(resid) %*% S_i_grad_inv %*% resid_grad -
      0.5 * sum(diag(S_i_inv %*% S_i_grad2)) +
      0.5 * sum(diag(S_i_grad_inv %*% S_i %*% S_i_grad_inv %*% S_i))
    
    ll_grad2 = ll_grad2 + conditional_ll_grad2[i]
  }
  
  return(list(ll_grad2=ll_grad2,
              conditional_ll_grad2 = conditional_ll_grad2))
}

ll_grad2_spectral = function(log_lik_grad2=0,conditional_log_lik_grad2=rep(0,t),
                             w_FFT_real,m_i_i1,R_i_i1,m_i_i1_grad,R_i_i1_grad,
                             m_i_i1_grad2,R_i_i1_grad2,t,K,tau2,R_grad){
  
  log_lik_grad2 = .C("ll_grad2_spectral",ll_grad2=as.double(log_lik_grad2),
                     conditional_ll_grad2=as.double(conditional_log_lik_grad2),
                     y_FFT=as.double(w_FFT_real), 
                     mu_i_i1=as.double(m_i_i1),
                     Sigma_i_i1=as.double(R_i_i1), 
                     mu_i_i1_grad=as.double(m_i_i1_grad),
                     Sigma_i_i1_grad=as.double(R_i_i1_grad), 
                     mu_i_i1_grad2=as.double(m_i_i1_grad2),
                     Sigma_i_i1_grad2=as.double(R_i_i1_grad2),
                     T=as.integer(t),K = as.integer(K), 
                     tau2 = as.double(tau2), 
                     R_grad = as.double(R_grad))
  
  output = list(ll_grad2 = log_lik_grad2$ll_grad2,
                conditional_ll_grad2 = log_lik_grad2$conditional_ll_grad2)
  return(output)
}

###########################

## ---
###################################
# CONDITIONAL GRAD LOG-LIKELIHOOD #
###################################

conditional_log_likelihood_grad <- function(param_obs = NULL, param_sig = NULL,
                                            w = NULL, w_FFT_real = NULL, n, t,
                                            K=n^2, dt = 1, m0 = rep(0,K), 
                                            R0 = rep(0,K), 
                                            m0_grad = rep(list(rep(0,K)),9), 
                                            R0_grad = rep(list(rep(0,K)),9),
                                            nu = 1, grad_params=1:9, x=NULL,
                                            log_scale = FALSE, 
                                            log_indices = c(1,2,3,4,5,9),
                                            negative = FALSE){
  
  ## Parameters ##
  param = c(param_sig,as.numeric(unlist(param_obs$tau2_vals)))
  
  ## If parameters provided on log-scale, rescale
  if(log_scale){
    param[log_indices] <- exp(param[log_indices])
  }
  
  ## Signal and observation parameters ##
  param_sig = param[1:8]
  param_obs = param_obs_func(1,list(1:n^2),list(param[9]))
  
  ## Compute linear predictor X*beta
  if(!is.null(x)){
    linear_pred <- apply(x,c(2,3),linear_predictor,beta = param[-c(1:9)])
  }
  else{
    linear_pred <- 0
  }
  
  ## If 'w_FFT_real' not given, compute as inverse-FFT of (w-x*beta)
  if(is.null(w_FFT_real)){
    w_FFT_real <- FFT_real_ST(STmatrix_to_vec(w-linear_pred),n=n,t=t,inv=TRUE)
  }
  
  ## Compute log-likelihood using 'Multi_TKF_spectral'
  conditional_ll_grad <- Multi_TKF_spectral(w=w,w_FFT_real=w_FFT_real,
                                            param_obs = param_obs,
                                            param_sig = param_sig,n=n,t=t,K=n^2,
                                            log_lik=TRUE,dt=dt,m0=m0,R0=R0,
                                            m0_grad=m0_grad,R0_grad=R0_grad,
                                            nu=nu,grad_params=grad_params,
                                            log_lik_grad=TRUE)$conditional_ll_grad
  
  ## Output
  if(negative){
    return(lapply(conditional_ll_grad, function(x) -x))
  }
  else{
    return(conditional_ll_grad)
  }
}

###################################

#######################################
# CONDITIONAL 2nd GRAD LOG-LIKELIHOOD #
#######################################

conditional_log_likelihood_grad2 <- function(param_obs = NULL, param_sig = NULL,
                                             w = NULL, w_FFT_real = NULL, n, t,
                                             K = n^2, dt = 1, m0 = rep(0,K), 
                                             R0 = rep(0,K), 
                                             m0_grad = rep(list(rep(0,K)),9), 
                                             R0_grad = rep(list(rep(0,K)),9),
                                             m0_grad2 = rep(list(rep(0,K)),9), 
                                             R0_grad2 = rep(list(rep(0,K)),9),
                                             nu = 1, grad_params=1:9, x=NULL,
                                             log_scale=FALSE, 
                                             log_indices=c(1,2,3,4,5,9),
                                             negative=FALSE){
  
  ## Parameters ##
  param = c(param_sig,as.numeric(unlist(param_obs$tau2_vals)))
  
  ## If parameters provided on log-scale, rescale
  if(log_scale){
    param[log_indices] <- exp(param[log_indices])
  }
  
  ## Signal and observation parameters ##
  param_sig = param[1:8]
  param_obs = param_obs_func(1,list(1:n^2),list(param[9]))
  
  ## Compute linear predictor X*beta
  if(!is.null(x)){
    linear_pred <- apply(x,c(2,3),linear_predictor,beta = param[-c(1:9)])
  }
  else{
    linear_pred <- 0
  }
  
  ## If 'w_FFT_real' not given, compute as inverse-FFT of (w-x*beta)
  if(is.null(w_FFT_real)){
    w_FFT_real <- FFT_real_ST(STmatrix_to_vec(w-linear_pred),n=n,t=t,inv=TRUE)
  }
  
  ## Compute log-likelihood using 'Multi_TKF_spectral'
  conditional_ll_grad2 <- Multi_TKF2_spectral(w=w,w_FFT_real=w_FFT_real,
                                              param_obs=param_obs,
                                              param_sig=param_sig,
                                              n=n,t=t,K=n^2,
                                              log_lik=TRUE,dt=dt,m0=m0,R0=R0,
                                              m0_grad=m0_grad,R0_grad=R0_grad,
                                              m0_grad2=m0_grad2,
                                              R0_grad2=R0_grad2,nu=nu,
                                              grad_params=grad_params,
                                              log_lik_grad=TRUE,
                                              log_lik_grad2=TRUE)$conditional_ll_grad2
  
  ## Output
  if(negative){
    return(lapply(conditional_ll_grad2, function(x) -x))
  }
  else{
    return(conditional_ll_grad2)
  }
}

#######################################

# ----

############################
# BATCH MAXIMUM LIKELIHOOD #
############################

BML <- function(y, n = sqrt(dim(y)[2]), t = dim(y)[1], K = dim(G)[1], 
                param_obs0, param_sig0, param_bias0, n_param, grad_params, 
                step_size, log_scale = FALSE, log_indices, 
                hessian_indices = c(3,4,6:8), n_iterations, n_obs, coords, 
                spec_weights_mat = NULL, plot = FALSE, plot_limits = NULL, 
                param_obs_true = NULL, param_sig_true = NULL, 
                param_bias_true = NULL,param_obs_names = NULL, 
                param_sig_names = NULL, param_bias_names = NULL, 
                param_est_starts, param_est_ends, plot_point_freq = 1, 
                plot_freq = 10, save_plot = FALSE, 
                plot_filename = "BML_plot.pdf", save = FALSE,
                filename = "BML.pdf", mean_param = FALSE, 
                print_iter_freq = 10, iter_numbers = 1:n_iterations,
                mu0 = NULL, Sigma0 = NULL, mu0_grad_param = NULL, 
                Sigma0_grad_param = NULL, mu0_grad_space = NULL,
                Sigma0_grad_space = NULL,mu0_grad2_param = NULL, 
                Sigma0_grad2_param = NULL, back_track = FALSE,
                back_track_tol = 1e-3, back_track_independent = TRUE,
                back_track_time_independent = FALSE, back_track_scale = 0.5,
                back_track_param = 1:n_param, plot_vertical = TRUE,
                W_fourier = NULL, dt = 1){
  
  ## Number of parameters ##
  n_param_sig = length(param_sig)
  n_param_obs = param_obs0$m
  n_param_bias = param_bias0$m
  
  ## Parameter indices
  param_sig_indices = 1:n_param_sig
  param_obs_indices = (n_param_sig+1):(n_param_sig+n_param_obs)
  if(n_param_bias>0){param_bias_indices = (n_param_sig+n_param_obs+1):n_param}
  
  ## Matrix of parameter estimates ##
  param_est = matrix(0,ncol=n_param,nrow=n_iterations+1)
  
  ## Mean parameter estimates ##
  if(mean_param){
    
    ## Indices to compute mean parameter estimate(s) ##
    mean_param_est_indices = rep(list(c()),n_param)
    
    ## Mean parameter estimate(s) ##
    mean_param_est = rep(list(rep(0,100)),n_param)
    
    ## Start indices for mean parameter estimate(s) ##
    which_start = rep(list(c()),n_param)
    
    ## Start(s) for mean parameter estimate plot(s) ##
    mean_param_est_plot_start = rep(list(rep(0,100)),n_param)
    
  }
  
  ## Initial parameter estimate ##
  param0 = c(param_sig0,as.numeric(unlist(param_obs0$tau2_vals)))
  if(n_param_bias>0){param0 = c(param0,as.numeric(unlist(param_bias0$bias_vals)))}
  
  ## Rescale if param0 on log scale ##
  if(log_scale){param0[log_indices] <- exp(param0[log_indices])}
  
  ## Add initial parameter estimate to matrix of parameter estimates ##
  param_est[1,] = param0
  
  ## Matrices of log-lik,log-lik gradients, 2nd derivatives ##
  ll = rep(0,n_iterations)
  ll_grads_param = matrix(0,ncol=n_param,nrow=n_iterations)
  ll_grads_space = rep(list(matrix(0,ncol=2,nrow=n_iterations)),n_obs)
  ll_grads2_param = matrix(0,ncol=n_param,nrow=n_iterations)
  
  ## Matrices of obj. func., conditional obj. func. grads ##
  obj_func = matrix(0,ncol=1,nrow=n_iterations)
  obj_func_grads_param = matrix(0,ncol=n_param,nrow=n_iterations)
  obj_func_grads_space = rep(list(matrix(0,ncol=2,nrow=n_iterations)),n_obs)
  
  
  ## Hessian & Non-Hessian parameters to estimate ##
  grad_params_hessian = intersect(hessian_indices,grad_params)
  grad_params_n_hessian = c(setdiff(hessian_indices,grad_params),
                            setdiff(grad_params,hessian_indices))
  
  ## SPDE_FT object ##
  SPDE_FT = spde_initialise(n,t,K)
  wave = SPDE_FT$wave
  cosine_indices = SPDE_FT$cosine_indices
  n_cos = SPDE_FT$n_cos
  n_cos_sin = length(cosine_indices)
  K = SPDE_FT$K
  
  ## Initial values ##
  if(is.null(mu0)){
    mu0 = rep(0,K)
  }
  if(is.null(Sigma0)){
    Sigma0 = Q_f(param_sig0,wave,n,n_cos=n_cos,nu=1,dt=dt,norm=TRUE,
                 spec_weights_mat = spec_weights_mat)
  }
  if(is.null(mu0_grad_param)){
    mu0_grad_param = rep(list(rep(0,K)),n_param)
  }
  if(is.null(Sigma0_grad_param)){
    Sigma0_grad_param = Q_grad_f(param_obs0,param_sig0,n,K,wave,n_cos,nu,dt=dt,norm=TRUE,
                           param_bias0,spec_weights_mat = spec_weights_mat)
  }
  if(is.null(mu0_grad_space)){
    mu0_grad_space = rep(list(rep(0,K)),2*n_obs)
  }
  if(is.null(Sigma0_grad_space)){
    Sigma0_grad_space = Q_grad_space_f(n_obs,K)
  }
  if(!is.null(hessian_indices)){
    if(is.null(mu0_grad2_param)){
      mu0_grad2_param = rep(list(rep(0,K)),n_param)
    }
    if(is.null(Sigma0_grad2_param)){
      Sigma0_grad2_param = Q_grad2_f(param_obs0,param_sig0,n,K,wave,n_cos,nu,dt=dt,
                               norm=TRUE,param_bias0,
                               spec_weights_mat = spec_weights_mat)
    }
  }
  
  ## Compute Phi, Phi_grad, Phi_grad2 ##
  Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs,coords)
  Phi_grad_param = Phi_grad_f(param_obs0,param_sig0,n,K,param_bias0,n_obs)
  Phi_grad_space = Phi_grad_space_f(wave,cosine_indices,n_cos,n,n_obs,coords=coords,grad_coord_index = 1:n_obs)
  Phi_grad2_param = Phi_grad2_f(param_obs0,param_sig0,n,K,param_bias0,n_obs)
  
  
  ### BML ###
  for (i in 1:n_iterations){
    
    ## Current obs & sig parameters ##
    param_obs = param_obs_func(param_obs0$m,param_obs0$m_indices,
                               param_est[i,param_obs_indices])
    param_sig = param_est[i,param_sig_indices]
    
    ## Current bias parameters ##
    if(n_param_bias>0){
      param_bias = param_bias_func(param_bias0$m,param_bias0$m_indices,
                                   param_est[i,param_bias_indices])
    }
    else{
      param_bias = param_bias_func()
    }
    
    ## LGSSM matrices ##
    XBeta = XBeta_f(param_bias,n,t,n_obs)
    R = R_f(param_obs,n,n_obs)
    G = G_f(param_sig,wave,cosine_indices,dt=dt,n_cos=n_cos)
    Q = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=dt,norm=TRUE,
            spec_weights_mat = spec_weights_mat)
    
    ## LGSSM (parameter) gradient arrays ##
    XBeta_grad_param = XBeta_grad_f(param_obs,param_sig,n,t,param_bias,n_obs)
    R_grad_param = R_grad_f(param_obs,param_sig,n,param_bias,n_obs)
    G_grad_param = G_grad_f(param_obs,param_sig,n,K,wave,cosine_indices,dt=dt,nu,n_cos,
                      param_bias)
    Q_grad_param = Q_grad_f(param_obs,param_sig,n,K,wave,n_cos,nu,dt=dt,norm=TRUE,
                      param_bias,spec_weights_mat = spec_weights_mat)
    
    ## LGSSM (space) gradient arrays ##
    XBeta_grad_space = XBeta_grad_space_f(n_obs,t)
    R_grad_space = R_grad_space_f(n_obs)
    G_grad_space = G_grad_space_f(n_obs,K)
    Q_grad_space = Q_grad_space_f(n_obs,K)
    
    ## LGSSM 2nd (parameter) derivative arrays ##
    if(!is.null(hessian_indices)){
      XBeta_grad2_param = XBeta_grad2_f(param_obs,param_sig,n,t,param_bias,n_obs)
      R_grad2_param = R_grad2_f(param_obs,param_sig,n,param_bias,n_obs)
      G_grad2_param = G_grad2_f(param_obs,param_sig,n,K,wave,cosine_indices,dt=dt,nu,
                          n_cos,param_bias)
      Q_grad2_param = Q_grad2_f(param_obs,param_sig,n,K,wave,n_cos,nu,dt=dt,norm=TRUE,
                          param_bias,spec_weights_mat = spec_weights_mat)
    }
    
    ## Run TKF (parameters) ##
    
    if(!is.null(hessian_indices)){
      
      ## Run MTKF2 (parameters) ##
      MTKF_param = MultiTangentKalmanFilter2(y = y, Phi = Phi, 
                                       R = R, G = G, Q = Q, 
                                       Phi_grad = Phi_grad_param, 
                                       R_grad = R_grad_param, 
                                       G_grad = G_grad_param, 
                                       Q_grad = Q_grad_param, 
                                       Phi_grad2 = Phi_grad2_param, 
                                       R_grad2 = R_grad2_param, 
                                       G_grad2 = G_grad2_param, 
                                       Q_grad2 = Q_grad2_param, n = n_obs, 
                                       t = t, K = K, log_lik = TRUE, mu0 = mu0, 
                                       Sigma0 = Sigma0, 
                                       mu0_grad = mu0_grad_param,
                                       Sigma0_grad = Sigma0_grad_param, 
                                       mu0_grad2 = mu0_grad2_param,
                                       Sigma0_grad2 = Sigma0_grad2_param,
                                       n_param = n_param,
                                       grad_params = grad_params_hessian,
                                       lin_pred = XBeta,
                                       lin_pred_grad = XBeta_grad_param,
                                       lin_pred_grad2 = XBeta_grad2_param)
    }
    
    else{
      
      ## Run MTKF (parameters) ##
      MTKF_param = MultiTangentKalmanFilter(y = y, Phi = Phi, 
                                      R = R, G = G, Q = Q, 
                                      Phi_grad = Phi_grad_param, 
                                      R_grad = R_grad_param, 
                                      G_grad = G_grad_param, 
                                      Q_grad = Q_grad_param, 
                                      n = n_obs, t = t, K = K, log_lik = TRUE, 
                                      mu0 = mu0, Sigma0 = Sigma0, 
                                      mu0_grad = mu0_grad_param,
                                      Sigma0_grad = Sigma0_grad_param, 
                                      n_param = n_param,
                                      grad_params = grad_params,
                                      lin_pred = XBeta,
                                      lin_pred_grad = XBeta_grad_param)
      
    }
      
    ## Log-likelihood ##
    ll[i] = unlist(MTKF_param$ll)
    
    ## Log-likelihood (parameter) gradient ##
    ll_grads_param[i,] = unlist(MTKF_param$ll_grad)
    
    ## Log likelihood second (parameter) gradient ##
    if(!is.null(hessian_indices)){
      ll_grads2_param[i,grad_params_hessian] = 
        unlist(MTKF_param$ll_grad2[grad_params_hessian])
    }
    
    ## Objective function ##
    obj_func[i] = objective_function(t,MTKF_param$Sigma_i_i,W_fourier)
    
    # Objective function (parameter) gradient ##
    obj_func_grads_param[i,] = 
      objective_function_grad_param(n_param,t,MTKF_param$Sigma_i_i_grad,
                                    W_fourier)
    
    
    ## Run TKF (space) ##
    
    ## Run MTKF (space) ##
    MTKF_space = MultiTangentKalmanFilter(y = y, Phi = Phi, 
                                          R = R, G = G, Q = Q, 
                                          Phi_grad = Phi_grad_space, 
                                          R_grad = R_grad_space, 
                                          G_grad = G_grad_space, 
                                          Q_grad = Q_grad_space, 
                                          n = n_obs, t = t, K = K, log_lik = TRUE, 
                                          mu0 = mu0, Sigma0 = Sigma0, 
                                          mu0_grad = mu0_grad_space,
                                          Sigma0_grad = Sigma0_grad_space, 
                                          n_param = 2*n_obs,
                                          grad_params = 1:(2*n_obs),
                                          lin_pred = XBeta,
                                          lin_pred_grad = XBeta_grad_space)
    
    ## Log likelihood (spatial) gradient ##
    for (k in 1:n_obs){
      ll_grads_space[[k]][i,] = unlist(MTKF_space$ll_grad[(2*k-1):(2*k)])
    }
    
    ## Objective function (spatial) gradient ##
    obj_func_grad = objective_function_grad_space(n_obs,t,MTKF_space$Sigma_i_i_grad,
                                                  W_fourier)
    for (k in 1:n_obs){
      obj_func_grads_space[[k]][i,] = obj_func_grad[[k]]
    }
    
    #########################################################
    ## ADD GRADIENT OF LIKELIHOOD AND OBJ. FUNC. WRT SPACE ##
    #########################################################
    
    ## Parameter update ##
    
    ## Rescale to log-scale ##
    param_est[i+1,] = param_est[i,]
    param_est[i+1,log_indices] = log(param_est[i+1,log_indices])
    
    ## Gradient ascent ##
    param_est[i+1,grad_params_n_hessian] = 
      param_est[i+1,grad_params_n_hessian] + 
      step_size[i,grad_params_n_hessian]*ll_grads_param[i,grad_params_n_hessian]
    
    ## Gradient ascent (w/ Hessian) ##
    param_est[i+1,grad_params_hessian] =   
      param_est[i+1,grad_params_hessian] - 
      step_size[i,grad_params_hessian]*
      ll_grads2_param[i,grad_params_hessian]^(-1)*
      ll_grads_param[i,grad_params_hessian]
    
    ## Rescale to standard scale ##
    param_est[i+1,log_indices] = exp(param_est[i+1,log_indices])
    
    ## Ensure alpha is in [0,pi/2] ##
    if(param_est[i+1,6]<0){param_est[i+1,6] = param_est[i,6]/2}
    if(param_est[i+1,6]>=pi/2){param_est[i+1,6] = (param_est[i,6]/2+pi/2)/2}
    
    ## Ensure mu_x is in [-0.5,0.5] ##
    if(param_est[i+1,7]<(-0.5)){param_est[i+1,7] = param_est[i,7]/2}
    if(param_est[i+1,7]>0.5){param_est[i+1,7] = param_est[i,7]/2}
    
    ## Ensure mu_y is in [-0.5,0.5] ##
    if(param_est[i+1,8]<(-0.5)){param_est[i+1,8] = param_est[i,8]/2}
    if(param_est[i+1,8]>0.5){param_est[i+1,8] = param_est[i,8]/2}
    
    
    
    ## Back Tracking Line Search ##
    if(back_track){
      
      ## Back Tracking Line Search For Each Parameter ##
      if(back_track_independent){
        
        checked = rep(FALSE,n_param)
        
        ## Loop over parameters ##
        for(l in back_track_param){
          
          while(!checked[l]){
            
            ## Parameter vector with only l^th parameter updated ##
            param_check = param_est[i,]
            param_check[l] = param_est[i+1,l]
            
            ## Observation, signal, bias parameters ##
            param_obs_new = param_obs_func(param_obs0$m,param_obs0$m_indices,
                                           param_check[param_obs_indices])
            param_sig_new = param_check[param_sig_indices]
            
            if(n_param_bias>0){
              param_bias_new = param_bias_func(param_bias0$m,
                                               param_bias0$m_indices,
                                               param_check[param_bias_indices])
            }
            else{
              param_bias_new = param_bias_func()
            }
            
            ## LGGSM matrices for this parameter vector ##
            Phi_new = Phi_f(wave,cosine_indices,n_cos,n,n_obs,coords)
            XBeta_new = XBeta_f(param_bias_new,n,t,n_obs)
            R_new = R_f(param_obs_new,n,n_obs)
            G_new = G_f(param_sig_new,wave,cosine_indices,dt=dt,n_cos=n_cos)
            Q_new = Q_f(param_sig_new,wave,n,n_cos=n_cos,nu=1,dt=dt,norm=TRUE,
                        spec_weights_mat = spec_weights_mat)
            
            ## Compute likelihood for this parameter vector ##
            ll_update = KalmanFilter(y = y, Phi = Phi_new, R = R_new, 
                                     G = G_new, Q = Q_new, m0 = mu0, 
                                     R0 = Sigma0, linear_pred = XBeta_new, 
                                     N = n_obs, t = t, K = K, 
                                     log_lik = TRUE)$ll
            
            ## If new likelihood > old likelihood, accept new parameter ##
            if(ll_update > ll[i]){
              checked[l] = TRUE
            }
            
            checked[l] = TRUE
            
            ## If old likelihood >~ new likelihood ##
            if((ll[i]-ll_update)>back_track_tol){
              
              ## Half future step-sizes ##
              if(!back_track_time_independent){
                step_size[(i+1):(n_iterations),l] = 
                  back_track_scale*step_size[(i+1):(n_iterations),l]
              }
              
              ## Re-update new parameter ##
              param_est[i+1,l] = 
                back_track_scale*((1/back_track_scale-1)*param_est[i,l] + param_est[i+1,l])
              
              ## Check again ##
              checked[l] = FALSE
            }
          }
        }
      }
      
      ## Back Tracking Line Search For All Parameters ##
      if(!back_track_independent){
        
        checked = FALSE
        
        while(!checked){
          
          ## Parameter vector with all elements updated ##
          param_check = param_est[i+1,]
          
          ## Observation, signal, bias parameters ##
          param_obs_new = param_obs_func(param_obs0$m,param_obs0$m_indices,
                                         param_check[param_obs_indices])
          param_sig_new = param_check[param_sig_indices]
          
          if(n_param_bias>0){
            param_bias_new = param_bias_func(param_bias0$m,
                                             param_bias0$m_indices,
                                             param_check[param_bias_indices])
          }
          else{
            param_bias_new = param_bias_func()
          }
          
          ## LGGSM matrices for this parameter vector ##
          Phi_new = Phi_f(wave,cosine_indices,n_cos,n,n_obs,coords)
          XBeta_new = XBeta_f(param_bias_new,n,t,n_obs)
          R_new = R_f(param_obs_new,n,n_obs)
          G_new = G_f(param_sig_new,wave,cosine_indices,dt=dt,n_cos=n_cos)
          Q_new = Q_f(param_sig_new,wave,n,n_cos=n_cos,nu=1,dt=dt,norm=TRUE,
                      spec_weights_mat = spec_weights_mat)
          
          ## Compute likelihood for this parameter vector ##
          ll_update = KalmanFilter(y = y, Phi = Phi_new, R = R_new, 
                                   G = G_new, Q = Q_new, m0 = mu0, 
                                   R0 = Sigma0, linear_pred = XBeta_new, 
                                   N = n_obs, t = t, K = K, log_lik = TRUE)$ll
          
          ## If new likelihood > old likelihood, accept new parameter ##
          if(ll_update > ll[i]){
            checked = TRUE
          }
          
          checked = TRUE
          
          ## If old likelihood >~ new likelihood ## 
          if((ll[i]-ll_update)>back_track_tol){
            
            ## Half step-size ##
            if(!back_track_time_independent){
              step_size[(i+1):(n_iterations),] = 
                back_track_scale*step_size[(i+1):(n_iterations),]
            }
            
            ## Re-update new parameter ##
            param_est[i+1,] = 
              back_track_scale*((1/back_track_scale-1)*param_est[i,] + param_est[i+1,])
            
            ## Check again ##
            checked = FALSE
          }
        }
      }
    }
    
    
    ## Update 'initial' values for TKF ##
    mu0 = MTKF_param$mu_i_i[t+1,]
    Sigma0 = MTKF_param$Sigma_i_i[t+1,,]
    mu0_grad_param = lapply(MTKF_param$mu_i_i_grad, function(A) A[t+1,])
    Sigma0_grad_param = lapply(MTKF_param$Sigma_i_i_grad, function(A) A[t+1,,])
    mu0_grad_space = lapply(MTKF_space$mu_i_i_grad, function(A) A[t+1,])
    Sigma0_grad_space = lapply(MTKF_space$Sigma_i_i_grad, function(A) A[t+1,,])
    if(!is.null(hessian_indices)){
      mu0_grad2_param = lapply(MTKF_param$mu_i_i_grad2, function(A) A[t+1,])
      Sigma0_grad2_param = lapply(MTKF_param$Sigma_i_i_grad2, function(A) A[t+1,,])
    }
    
    ## Print iteration ##
    if(!is.null(print_iter_freq)){
      if(iter_numbers[i]%%print_iter_freq==0){
        cat("Iteration ", iter_numbers[i], "\n")
      }
    }
    
    ## Mean parameter estimates ##
    if(mean_param){
      
      if(!is.null(param_est_starts)){
        
        for (k in 1:n_param){
          
          if(i>(param_est_starts[[k]][1]+1)){
            
            ## Which start value to use to compute mean parameter estimate ##
            which_start[[k]] = max(which(i>(param_est_starts[[k]]+1)))
            
            ## Which start value to use to plot mean parameter estimate ##
            if (k %in% param_sig_indices){
              mean_param_est_plot_start[[k]][which_start[[k]]] = param_sig_t_dep$t_vals[[max(which(unlist(param_sig_true$t_vals)<param_est_starts[[k]][which_start[[k]]]))]]
            }
            if (k %in% param_obs_indices){
              mean_param_est_plot_start[[k]][which_start[[k]]] = param_obs_t_dep$t_vals[[max(which(unlist(param_obs_true$t_vals)<param_est_starts[[k]][which_start[[k]]]))]]
            }
            if(n_param_bias>0){
              if (k %in% param_bias_indices){
                mean_param_est_plot_start[[k]][which_start[[k]]] = param_bias_t_dep$t_vals[[max(which(unlist(param_bias_true$t_vals)<param_est_starts[[k]][which_start[[k]]]))]]
              }
            }
            
            ## Compute mean parameter estimate indices ##
            mean_param_est_indices[[k]] = 
              (param_est_starts[[k]][which_start[[k]]]+1):(min(i,param_est_ends[[k]][which_start[[k]]]))
            
            ## Compute mean parameter estimate ##
            mean_param_est[[k]][which_start[[k]]] = 
              mean(param_est[mean_param_est_indices[[k]],k])
          }
        }
      }
    }
    
    ## Plots ##
    if(plot){
      
      if(i%%plot_freq==0){
        
        ## Open plotting device ##
        if(save_plot){
          if(i==n_iterations){
            setwd(fig_wd)
            pdf(plot_filename,width=10,height=6.5)
            par(mar=c(2,4.5,1.5,1.5),mgp=c(1.5,.5,0))
          }
        }
        
        ## Number of Plots ##
        if(n_param_bias==0){par(mfrow=c(3,3))}
        if(n_param_bias>0){par(mfrow=c(3,4))}
        
        ## Plot Points ##
        plot_points = seq(1,i,plot_point_freq)
        
        ## Signal Parameters ##
        for (k in param_sig_indices){
          
          y_lim = NULL
          y_lab = expression(theta)
          main = NULL
          
          if(!is.null(plot_limits)){
            y_lim = plot_limits[[k]]
          }
          
          if(!is.null(param_sig_names)){
            y_lab = parse(text=paste("hat(",toString(param_sig_names[[k]]),")",sep=""))
            main = param_sig_names[[k]]
          }
          
          plot(plot_points, param_est[plot_points,k], cex=0.4,
               xlim = c(1,n_iterations+1), xlab = "t",ylab = y_lab,
               ylim = y_lim, main = main)
          
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
                       n_iterations+1,
                       param_sig_true$param_vals[[length(param_sig_true$t_vals)]][k],
                       col="red",lty=2)
            }
          }
          
          if(mean_param){
            if(i>(param_est_starts[[k]][1]+1)){
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
              y_lab = parse(text=paste("hat(",toString(param_obs_names[[k-8]]),")",sep=""))
              main = param_obs_names[[k-8]]
            }
            
            plot(plot_points, param_est[plot_points,k], cex=0.4,
                 xlim = c(1,n_iterations+1) ,xlab="t", ylab=expression(theta),
                 ylim = y_lim, main = main)
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
                       n_iterations+1,
                       param_obs_true$param_vals[[length(param_obs_true$t_vals)]]$tau2_vals[[k-8]],
                       col="red",lty=2)
            }
          }
          
          if(mean_param){
            if(i>(param_est_starts[[k]][1]+1)){
              for (m in 1:(which_start[[k]])){
                if(plot_vertical){
                  abline(v=param_est_starts[[k]][m],col="black",lty=2)
                }
                if(length(param_obs_true$t_vals)==1){
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
        
        ## Bias Parameters ##
        if(n_param_bias>0){
          for (k in param_bias_indices){
            if(k==param_bias_indices[[1]]){
              
              y_lim = NULL
              y_lab = expression(theta)
              main = NULL
              
              if(!is.null(plot_limits)){
                y_lim = plot_limits[[10]]
              }
              
              if(!is.null(param_bias_names)){
                y_lab = parse(text=paste("hat(",toString(param_obs_names[[k-(n_param_sig+n_param_obs)]]),")",sep=""))
                main = param_bias_names[[k-(n_param_sig+n_param_obs)]]
              }
              
              plot(plot_points, param_est[plot_points,k], cex=0.4,
                   xlim = c(1,n_iterations+1) ,xlab="t", ylab=expression(theta),
                   ylim = y_lim, main = main)
            }
            
            if(k>param_bias_indices[[1]]){
              points(plot_points,param_est[plot_points,k],pch=19,cex=0.4)
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
                         n_iterations+1,
                         param_bias_true$param_vals[[length(param_bias_true$t_vals)]]$bias_vals[[k-(n_param_sig+n_param_obs)]],
                         col="red",lty=2)
              }
            }
            
            if(mean_param){
              if(i>(param_est_starts[[k]][1]+1)){
                for (m in 1:(which_start[[k]])){
                  if(plot_vertical){
                    abline(v=param_est_starts[[k]][m],col="black",lty=2)
                  }
                  if(length(param_bias_true$t_vals)==1){
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
        }
        
        ## Close plotting device ##
        if(save_plot){
          if(i==n_iterations){
            dev.off()
          }
        }
      }
    }
  }
  
  ## Ouput ##
  output = list(param_est = param_est,ll_grads_param = ll_grads_param,
                ll_grads_space = ll_grads_space, 
                obj_func_grads_param = obj_func_grads_param,
                obj_func_grads_space = obj_func_grads_space)
  
  if(!is.null(hessian_indices)){
    output = c(output, list(ll_grads2_param = ll_grads2_param))
  }
  
  output = c(output,list(mu0 = mu0, Sigma0 = Sigma0, 
                         mu0_grad_param = mu0_grad_param,
                         Sigma0_grad_param = Sigma0_grad_param,
                         mu0_grad_space = mu0_grad_space,
                         Sigma0_grad_space = Sigma0_grad_space))
  
  if(!is.null(grad_params_hessian)){
    output = c(output,list(mu0_grad2_param = mu0_grad2_param, 
                           Sigma0_grad2_param = Sigma0_grad2_param))
  }
  
  if(back_track){
    output = c(output,list(step_size = step_size))
  }
  
  ## Save ##
  if(save){
    saveRDS(output, file = filename)
  }

  return(output)
}

############################

#####################################
# SPECTRAL BATCH MAXIMUM LIKELIHOOD #
#####################################

BML_spectral <- function(w, n, t, K = n^2, x = NULL, param_obs0, param_sig0,
                         n_param, grad_params = 1:9, step_size, 
                         log_scale = FALSE, log_indices, 
                         hessian_indices = c(3,4,6:8), n_iterations,
                         plot = FALSE, plot_limits = NULL, 
                         param_obs_true = NULL, param_sig_true = NULL, 
                         param_obs_names = NULL, param_sig_names = NULL, 
                         param_est_starts, param_est_ends, 
                         plot_point_freq = 1, plot_freq = 10, save_plot = FALSE,
                         plot_filename = "BML_plot.pdf", save = FALSE,
                         filename = "BML.pdf", mean_param = FALSE, 
                         print_iter_freq = 100,
                         iter_numbers = 1:n_iterations, m0 = NULL, 
                         R0 = NULL, m0_grad_param = NULL, R0_grad_param = NULL,
                         m0_grad2_param = NULL, R0_grad2_param = NULL, 
                         back_track = FALSE,
                         back_track_tol = 1e-3, back_track_independent = TRUE,
                         back_track_time_independent = FALSE, back_track_scale = 0.5,
                         back_track_param = 1:n_param, plot_vertical = TRUE,
                         W_fourier = NULL, dt = 1){
  
  ## Number of parameters ##
  n_param_sig = length(param_sig0)
  n_param_obs = param_obs0$m
  
  ## Parameter indices ##
  param_sig_indices = 1:n_param_sig
  param_obs_indices = (n_param_sig+1):(n_param_sig+n_param_obs)
  
  ## Matrix of parameter estimates ##
  param_est = matrix(0,ncol=n_param,nrow=n_iterations+1)
  
  ## Mean parameter estimates ##
  if(mean_param){
    
    ## Indices to compute mean parameter estimate(s) ##
    mean_param_est_indices = rep(list(c()),n_param)
    
    ## Mean parameter estimate(s) ##
    mean_param_est = rep(list(rep(0,100)),n_param)
    
    ## Start indices for mean parameter estimate(s) ##
    which_start = rep(list(c()),n_param)
    
    ## Start(s) for mean parameter estimate plot(s) ##
    mean_param_est_plot_start = rep(list(rep(0,100)),n_param)
    
  }
  
  ## Initial parameter estimate ##
  param0 = c(param_sig0,as.numeric(unlist(param_obs0$tau2_vals)))
  
  ## Rescale if param0 on log scale ##
  if(log_scale){param0[log_indices] <- exp(param0[log_indices])}
  
  ## Add initial parameter estimate to matrix of parameter estimates ##
  param_est[1,] = param0
  
  ## Matrices of log-lik, log-lik gradients, 2nd derivatives ##
  ll = rep(0,n_iterations)
  ll_grads_param = matrix(0,ncol=n_param,nrow=n_iterations)
  ll_grads_space = rep(list(matrix(0,ncol=2,nrow=n_iterations)),n_obs)
  ll_grads2_param = matrix(0,ncol=n_param,nrow=n_iterations)
  
  ## Matrices of obj. func., conditional obj. func. grads ##
  obj_func = matrix(0,nrow=n_iterations,ncol=1)
  obj_func_grads_param = matrix(0,ncol=n_param,nrow=n_iterations)
  obj_func_grads_space = rep(list(matrix(0,ncol=2,nrow=n_iterations)),n_obs)
  
  
  ## Hessian & Non-Hessian parameters to estimate ##
  grad_params_hessian = intersect(hessian_indices,grad_params)
  grad_params_n_hessian = c(setdiff(hessian_indices,grad_params),
                           setdiff(grad_params,hessian_indices))
  
  ## SPDE_FT object ##
  SPDE_FT = spde_initialise(n,t,K)
  wave = SPDE_FT$wave
  cosine_indices = SPDE_FT$cosine_indices
  n_cos = SPDE_FT$n_cos
  n_cos_sin = length(cosine_indices)
  K = SPDE_FT$K
  
  ## Initial values ##
  if(is.null(m0)){
    m0 = rep(0,K)
  }
  if(is.null(R0)){
    R0 = Q_f_vec(param_sig0,wave=wave,n=n,n_cos=n_cos,nu=1,dt=dt,norm=TRUE)
  }
  if(is.null(m0_grad_param)){
    m0_grad_param = rep(list(rep(0,K)),n_param)
  }
  if(is.null(R0_grad_param)){
    R0_grad_param = Q_grad_f_vec(param_obs0,param_sig0,n=n,K=K,wave=wave,n_cos=n_cos,
                           nu=1,dt=dt,norm=TRUE)
  }
  if(is.null(m0_grad2_param)){
    m0_grad2_param = rep(list(rep(0,K)),n_param)
  }
  if(is.null(R0_grad2_param)){
    R0_grad2_param = Q_grad2_f_vec(param_obs0,param_sig0,n=n,K=K,wave=wave,
                             n_cos=n_cos,nu=1,dt=dt,norm=TRUE)
  }
  
  
  ## Run BML ##
  
  for (i in 1:n_iterations){
    
    ## Current obs & sig parameters ##
    param_obs = param_obs_func(param_obs0$m,param_obs0$m_indices,
                               param_est[i,param_obs_indices])
    param_sig = param_est[i,param_sig_indices]
    
    if(length(grad_params_hessian)>0){
      
      ## Run MTKF2 ##
      MTKF2_param = Multi_TKF2_spectral(w = w, param_obs = param_obs,
                                  param_sig = param_sig, n = n, t = t, K = K,
                                  log_lik = TRUE, dt = dt, m0 = m0, R0 = R0,
                                  m0_grad = m0_grad_param, 
                                  R0_grad = R0_grad_param,
                                  m0_grad2 = m0_grad2_param, 
                                  R0_grad2 = R0_grad2_param,
                                  grad_params = grad_params_hessian,
                                  log_lik_grad = TRUE, log_lik_grad2 = TRUE)
      
      ## Log Likelihood ##
      ll[i] = unlist(MTKF2_param$ll)
      
      ## Log Likelihood Gradient ##
      ll_grads_param[i,grad_params_hessian] = 
        unlist(MTKF2_param$ll_grad)[grad_params_hessian]
      
      ## Log Likelihood 2nd Gradient ##
      ll_grads2_param[i,grad_params_hessian] = 
        unlist(MTKF2_param$ll_grad2)[grad_params_hessian]
      
      ## Objective function ##
      Sigma_i_i_array = array(0,dim=c(t+1,K,K))
      for (j in 1:(t+1)){Sigma_i_i_array[j,,]=diag(MTKF2_param$Sigma_i_i[j,])}
      obj_func[i] = objective_function(t,Sigma_i_i_array,W_fourier)
      
      ## Objective function (parameter) gradient ##
      Sigma_i_i_grad_param_array = rep(list(array(0,dim=c(t+1,K,K))),n_param)
      for (k in 1:n_param){
        for (j in 1:(t+1)){
          Sigma_i_i_grad_param_array[[k]][j,,]=diag(MTKF2_param$Sigma_i_i_grad[[k]][j,])
        }
      }
      obj_func_grads_param[i,] = 
        objective_function_grad_param(n_param,t,Sigma_i_i_grad_param_array,
                                      W_fourier)
      
      MTKF_param = MTKF2_param
    }
    
    if(length(grad_params_n_hessian)>0){
      
      ## Run MTKF ##
      MTKF_param = Multi_TKF_spectral(w=w,param_obs=param_obs,param_sig=param_sig,
                                n=n,t=t,K=K,log_lik=TRUE,dt=dt,m0=m0,R0=R0,
                                m0_grad=m0_grad_param,R0_grad=R0_grad_param,
                                grad_params=grad_params_n_hessian,
                                log_lik_grad=TRUE)
      
      ## Log Likelihood ##
      ll[i] = unlist(MTKF_param$ll)
      
      ## Log Likelihood Gradient ##
      ll_grads_param[i,grad_params_n_hessian] = 
        unlist(MTKF_param$ll_grad)[grad_params_n_hessian]
      
      ## Objective function ##
      Sigma_i_i_array = array(0,dim=c(t+1,K,K))
      for (j in 1:(t+1)){Sigma_i_i_array[j,,]=diag(MTKF_param$Sigma_i_i[j,])}
      obj_func[i] = objective_function(t,Sigma_i_i_array,W_fourier)
      
      ## Objective function (parameter) gradient ##
      Sigma_i_i_grad_param_array = rep(list(array(0,dim=c(t+1,K,K))),n_param)
      for (k in 1:n_param){
        for (j in 1:(t+1)){
          Sigma_i_i_grad_param_array[[k]][j,,]=diag(MTKF_param$Sigma_i_i_grad[[k]][j,])
        }
      }
      obj_func_grads_param[i,] = 
        objective_function_grad_param(n_param,t,Sigma_i_i_grad_param_array,
                                      W_fourier)
    }
    
    ## Parameter update ##
    
    ## Rescale to log-scale ##
    param_est[i+1,] = param_est[i,]
    param_est[i+1,log_indices] = log(param_est[i+1,log_indices])
    
    ## Gradient ascent ##
    param_est[i+1,grad_params_n_hessian] = 
      param_est[i+1,grad_params_n_hessian] + 
      step_size[i,grad_params_n_hessian]*ll_grads_param[i,grad_params_n_hessian]
    
    ## Gradient ascent (w/ Hessian) ##
    param_est[i+1,grad_params_hessian] =   
      param_est[i+1,grad_params_hessian] - 
      step_size[i,grad_params_hessian]*
      ll_grads2_param[i,grad_params_hessian]^(-1)*
      ll_grads_param[i,grad_params_hessian]
    
    ## Rescale to standard scale ##
    param_est[i+1,log_indices] = exp(param_est[i+1,log_indices])
    
    ## Ensure alpha is in [0,pi/2] ##
    if(param_est[i+1,6]<0){param_est[i+1,6] = param_est[i,6]/2}
    if(param_est[i+1,6]>=pi/2){param_est[i+1,6] = (param_est[i,6]/2+pi/2)/2}
    
    ## Ensure mu_x is in [-0.5,0.5] ##
    if(param_est[i+1,7]<(-0.5)){param_est[i+1,7] = param_est[i,7]/2}
    if(param_est[i+1,7]>0.5){param_est[i+1,7] = param_est[i,7]/2}
    
    ## Ensure mu_y is in [-0.5,0.5] ##
    if(param_est[i+1,8]<(-0.5)){param_est[i+1,8] = param_est[i,8]/2}
    if(param_est[i+1,8]>0.5){param_est[i+1,8] = param_est[i,8]/2}
    
    
    ## Back Tracking Line Search ##
    if(back_track){
      
      ## Back Tracking Line Search For Each Parameter ##
      if(back_track_independent){
        
        checked = rep(FALSE,n_param)
        
        ## Loop over parameters ##
        for(l in back_track_param){
          
          while(!checked[l]){
            
            ## Parameter vector with only l^th parameter updated ##
            param_check = param_est[i,]
            param_check[l] = param_est[i+1,l]
            
            ## Observation, signal, bias parameters ##
            param_obs_new = param_obs_func(param_obs0$m,param_obs0$m_indices,
                                           param_check[param_obs_indices])
            param_sig_new = param_check[param_sig_indices]
            
            ## Compute likelihood for this parameter vector ##
            ll_update = KF_spectral(w=w,param_obs=param_obs_new,
                                    param_sig=param_sig_new, n=n,t=t,K=K,
                                    log_lik=TRUE,m0=m0,R0=R0)$ll
            
            ## If new likelihood > old likelihood, accept new parameter ##
            if(ll_update > ll[i]){
              checked[l] = TRUE
            }
            
            checked[l] = TRUE
            
            ## If old likelihood >~ new likelihood ##
            if((ll[i]-ll_update)>back_track_tol){
              
              ## Half future step-sizes ##
              if(!back_track_time_independent){
                step_size[(i+1):(n_iterations),l] = 
                  back_track_scale*step_size[(i+1):(n_iterations),l]
              }
              
              ## Re-update new parameter ##
              param_est[i+1,l] = 
                back_track_scale*((1/back_track_scale-1)*param_est[i,l] + param_est[i+1,l])
              
              ## Check again ##
              checked[l] = FALSE
            }
          }
        }
      }
      
      ## Back Tracking Line Search For All Parameters ##
      if(!back_track_independent){
        
        checked = FALSE
        
        while(!checked){
          
          ## Parameter vector with all elements updated ##
          param_check = param_est[i+1,]
          
          ## Observation, signal, bias parameters ##
          param_obs_new = param_obs_func(param_obs0$m,param_obs0$m_indices,
                                         param_check[param_obs_indices])
          param_sig_new = param_check[param_sig_indices]
          
          ## Compute likelihood for this parameter vector ##
          ll_update = KF_spectral(w=w,param_obs=param_obs_new,
                                  param_sig=param_sig_new, n=n,t=t,K=K,
                                  log_lik=TRUE,m0=m0,R0=R0)$ll
          
          ## If new likelihood > old likelihood, accept new parameter ##
          if(ll_update > ll[i]){
            checked = TRUE
          }
          
          checked = TRUE
          
          ## If old likelihood >~ new likelihood ## 
          if((ll[i]-ll_update)>back_track_tol){
            
            ## Half step-size ##
            if(!back_track_time_independent){
              step_size[(i+1):(n_iterations),] = 
                back_track_scale*step_size[(i+1):(n_iterations),]
            }
            
            ## Re-update new parameter ##
            param_est[i+1,] = 
              back_track_scale*((1/back_track_scale-1)*param_est[i,] + param_est[i+1,])
            
            ## Check again ##
            checked = FALSE
          }
        }
      }
    }
    
    
    ## Update 'initial' values for KF ##
    m0 = MTKF_param$mu_i_i[t+1,]
    R0 = MTKF_param$Sigma_i_i[t+1,]
    m0_grad_param = lapply(MTKF_param$mu_i_i_grad, function(A) A[t+1,])
    R0_grad_param = lapply(MTKF_param$Sigma_i_i_grad, function(A) A[t+1,])
    
    if(length(grad_params_hessian)>0){
      m0_grad2_param = lapply(MTKF2_param$mu_i_i_grad2, function(A) A[t+1,])
      R0_grad2_param = lapply(MTKF2_param$Sigma_i_i_grad2, function(A) A[t+1,])
    }
    
    ## Print iteration ##
    if(!is.null(print_iter_freq)){
      if(iter_numbers[i]%%print_iter_freq==0){
        cat("Iteration ", iter_numbers[i], "\n")
      }
    }
    
    ## Mean parameter estimates ##
    if(mean_param){
      
      if(!is.null(param_est_starts)){
        
        for (k in 1:n_param){
          
          if(i>(param_est_starts[[k]][1]+1)){
            
            ## Which start value to use to compute mean parameter estimate ##
            which_start[[k]] = max(which(i>(param_est_starts[[k]]+1)))
            
            ## Which start value to use to plot mean parameter estimate ##
            if (k %in% param_sig_indices){
              mean_param_est_plot_start[[k]][which_start[[k]]] = param_sig_t_dep$t_vals[[max(which(unlist(param_sig_true$t_vals)<param_est_starts[[k]][which_start[[k]]]))]]
            }
            if (k %in% param_obs_indices){
              mean_param_est_plot_start[[k]][which_start[[k]]] = param_obs_t_dep$t_vals[[max(which(unlist(param_obs_true$t_vals)<param_est_starts[[k]][which_start[[k]]]))]]
            }
            
            ## Compute mean parameter estimate indices ##
            mean_param_est_indices[[k]] = 
              (param_est_starts[[k]][which_start[[k]]]+1):(min(i,param_est_ends[[k]][which_start[[k]]]))
            
            ## Compute mean parameter estimate ##
            mean_param_est[[k]][which_start[[k]]] = 
              mean(param_est[mean_param_est_indices[[k]],k])
          }
        }
      }
    }
    
    ## Plots ##
    if(plot){
      
      if(i%%plot_freq==0){
        
        ## Open plotting device ##
        if(save_plot){
          if(i==n_iterations){
            setwd(fig_wd)
            pdf(plot_filename,width=10,height=6.5)
            par(mar=c(2,4.5,1.5,1.5),mgp=c(1.5,.5,0))
          }
        }
        
        ## Number of Plots ##
        par(mfrow=c(3,3))
        
        ## Plot Points ##
        plot_points = seq(1,i,plot_point_freq)
        
        ## Signal Parameters ##
        for (k in param_sig_indices){
          
          y_lim = NULL
          y_lab = expression(theta)
          main = NULL
          
          if(!is.null(plot_limits)){
            y_lim = plot_limits[[k]]
          }
          
          if(!is.null(param_sig_names)){
            y_lab = parse(text=paste("hat(",toString(param_sig_names[[k]]),")",sep=""))
            main = param_sig_names[[k]]
          }
          
          plot(plot_points, param_est[plot_points,k], cex=0.4,
               xlim = c(1,n_iterations+1), xlab = "t",ylab = y_lab,
               ylim = y_lim, main = main)
          
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
                       n_iterations+1,
                       param_sig_true$param_vals[[length(param_sig_true$t_vals)]][k],
                       col="red",lty=2)
            }
          }
          
          if(mean_param){
            if(i>(param_est_starts[[k]][1]+1)){
              for (m in 1:(which_start[[k]])){
                if(plot_vertical){
                  abline(v=param_est_starts[[k]][m],col="black",lty=2)
                }
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
              y_lab = parse(text=paste("hat(",toString(param_obs_names[[k-8]]),")",sep=""))
              main = param_obs_names[[k-8]]
            }
            
            plot(plot_points, param_est[plot_points,k], cex=0.4,
                 xlim = c(1,n_iterations+1) ,xlab="t", ylab=expression(theta),
                 ylim = y_lim, main = main)
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
                       n_iterations+1,
                       param_obs_true$param_vals[[length(param_obs_true$t_vals)]]$tau2_vals[[k-8]],
                       col="red",lty=2)
            }
          }
          
          if(mean_param){
            if(i>(param_est_starts[[k]][1]+1)){
              for (m in 1:(which_start[[k]])){
                if(plot_vertical){
                  abline(v=param_est_starts[[k]][m],col="black",lty=2)
                }
                if(length(param_obs_true$t_vals)==1){
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
        
        ## Close plotting device ##
        if(save_plot){
          if(i==n_iterations){
            dev.off()
          }
        }
      }
    }
  }
  
  ## Ouput ##
  output = list(param_est=param_est,ll_grads_param=ll_grads_param)
  
  if(!is.null(hessian_indices)){
    output = c(output,list(ll_grads2_param = ll_grads2_param))
  }
  
  output = c(output,list(mu0 = m0, Sigma0 = R0, mu0_grad_param = m0_grad_param,
                         Sigma0_grad_param = R0_grad_param))
  
  if(!is.null(grad_params_hessian)){
    output = c(output,list(mu0_grad2_param = m0_grad2_param, 
                           Sigma0_grad2_param = R0_grad2_param))
  }
  
  if(back_track){
    output = c(output,list(step_size = step_size))
  }
  
  ## Save ##
  if(save){
    saveRDS(output, file = filename)
  }
  
  ## Output ##
  return(output)
  
}


#####################################

# ----

################################
# RECURSIVE MAXIMUM LIKELIHOOD #
################################

RML <- function(y, n = sqrt(dim(y)[2]), t = dim(y)[1], K = dim(G)[1], 
                param_obs0, param_sig0, param_bias0, n_param, grad_params, 
                step_size, log_scale = FALSE, log_indices, 
                hessian_indices = NULL, n_iterations = 1, n_obs, coords, 
                spec_weights_mat = NULL, plot = FALSE, plot_limits = NULL, 
                param_obs_true = NULL, param_sig_true = NULL, 
                param_bias_true = NULL, param_obs_names = NULL, 
                param_sig_names = NULL, param_bias_names = NULL, 
                param_est_starts, param_est_ends, plot_point_freq = 1, 
                plot_freq = 200, save_plot = FALSE,
                plot_filename = "RML_plot.pdf", save = FALSE,
                filename = "RML.pdf", mean_param = FALSE, 
                print_iter_freq = 100, 
                iter_numbers = 1:(n_iterations*t),
                mu0 = NULL, Sigma0 = NULL, mu0_grad_param = NULL, 
                Sigma0_grad_param = NULL, mu0_grad_space = NULL,
                Sigma0_grad_space = NULL, mu0_grad2_param = NULL, 
                Sigma0_grad2_param = NULL, back_track = FALSE,
                back_track_tol = 1e-3, back_track_independent = TRUE,
                back_track_time_independent = FALSE, back_track_scale = 0.5,
                back_track_param = 1:n_param,
                plot_vertical = TRUE, W_fourier = NULL, dt = 1){
  
  ## Number of parameters ##
  n_param_sig = length(param_sig0)
  n_param_obs = param_obs0$m
  n_param_bias = param_bias0$m
  
  ## Parameter indices ##
  param_sig_indices = 1:n_param_sig
  param_obs_indices = (n_param_sig+1):(n_param_sig+n_param_obs)
  if(n_param_bias>0){param_bias_indices = (n_param_sig+n_param_obs+1):n_param}
  
  ## Matrix of parameter estimates ##
  param_est = matrix(0,ncol=n_param,nrow=t*n_iterations+1)
  
  ## Matrix of KF means ##
  kf_means = matrix(0,ncol=K,nrow=t*n_iterations)
  
  ## Mean parameter estimates ##
  if(mean_param){
  
    ## Indices to compute mean parameter estimate(s) ##
    mean_param_est_indices = rep(list(c()),n_param)
  
    ## Mean parameter estimate(s) ##
    mean_param_est = rep(list(rep(0,100)),n_param)
  
    ## Start indices for mean parameter estimate(s) ##
    which_start = rep(list(c()),n_param)
  
    ## Start(s) for mean parameter estimate plot(s) ##
    mean_param_est_plot_start = rep(list(rep(0,100)),n_param)
    
  }
  
  ## Initial parameter estimate ##
  param0 = c(param_sig0,as.numeric(unlist(param_obs0$tau2_vals)))
  if(n_param_bias>0){param0 = c(param0,as.numeric(unlist(param_bias0$bias_vals)))}
  
  ## Rescale if param0 on log scale
  if(log_scale){param0[log_indices] <- exp(param0[log_indices])}
  
  ## Add initial param. estimate to matrix of param. estimates
  param_est[1,] = param0
  
  ## Matrices of log-lik, conditional log-lik gradients, 2nd derivatives ##
  ll = rep(0,t*n_iterations)
  ll_grads_param = matrix(0,ncol=n_param,nrow=t*n_iterations)
  ll_grads_space = rep(list(matrix(0,ncol=2,nrow=t*n_iterations)),n_obs)
  ll_grads2_param = matrix(0,ncol=n_param,nrow=t*n_iterations)
  
  ## Matrices of obj. func., conditional obj. func. grads ##
  obj_func = matrix(0,ncol=1,nrow=t*n_iterations)
  obj_func_grads_param = matrix(0,ncol=n_param,nrow=t*n_iterations)
  obj_func_grads_space = rep(list(matrix(0,ncol=2,nrow=t*n_iterations)),n_obs)
  
  ## Hessian & Non-Hessian parameters to estimate ##
  grad_params_hessian = intersect(hessian_indices,grad_params)
  grad_params_n_hessian = c(setdiff(hessian_indices,grad_params),
                            setdiff(grad_params,hessian_indices))
  
  ## SPDE_FT object ##
  SPDE_FT = spde_initialise(n,t,K)
  wave = SPDE_FT$wave
  cosine_indices = SPDE_FT$cosine_indices
  n_cos = SPDE_FT$n_cos
  n_cos_sin = length(cosine_indices)
  K = SPDE_FT$K
  
  ## Initial values ##
  if(is.null(mu0)){
    mu0 = rep(0,K)
  }
  if(is.null(Sigma0)){
    Sigma0 = Q_f(param_sig0,wave,n,n_cos=n_cos,nu=1,dt=dt,norm=TRUE,
                 spec_weights_mat = spec_weights_mat)
  }
  if(is.null(mu0_grad_param)){
    mu0_grad_param = rep(list(rep(0,K)),n_param)
  }
  if(is.null(Sigma0_grad_param)){
    Sigma0_grad_param = Q_grad_f(param_obs0,param_sig0,n,K,wave,n_cos,nu,dt=dt,
                                 norm=TRUE,param_bias0,
                                 spec_weights_mat = spec_weights_mat)
  }
  if(is.null(mu0_grad_space)){
    mu0_grad_space = rep(list(rep(0,K)),2*n_obs)
  }
  if(is.null(Sigma0_grad_space)){
    Sigma0_grad_space = Q_grad_space_f(n_obs,K)
  }
  if(length(grad_params_hessian)>0){
    if(is.null(mu0_grad2_param)){
      mu0_grad2_param = rep(list(rep(0,K)),n_param)
    }
    if(is.null(Sigma0_grad2)){
      Sigma0_grad2_param = Q_grad2_f(param_obs0,param_sig0,n,K,wave,n_cos,nu,
                                     dt=dt,norm=TRUE,param_bias0,
                                     spec_weights_mat = spec_weights_mat)
    }
  }
  
  ### RML ###
  for(j in 1:n_iterations){
    
    for (i in 1:t){
      
      ## Current obs & sig parameters ##
      param_obs = param_obs_func(param_obs0$m,param_obs0$m_indices,
                                 param_est[(j-1)*t+i,param_obs_indices])
      param_sig = param_est[(j-1)*t+i,param_sig_indices]
      
      if(n_param_bias>0){
        param_bias = param_bias_func(param_bias0$m,param_bias0$m_indices,
                                    param_est[(j-1)*t+i,param_bias_indices])
      }
      else{
        param_bias = param_bias_func()
      }
      
      ## LGSSM matrices ##
      Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs,coords)
      XBeta = XBeta_f(param_bias,n,t,n_obs)
      R = R_f(param_obs,n,n_obs)
      G = G_f(param_sig,wave,cosine_indices,dt=dt,n_cos=n_cos)
      Q = Q_f(param_sig,wave,n,n_cos=n_cos,nu=1,dt=dt,norm=TRUE,
              spec_weights_mat = spec_weights_mat)
      
      ## LGSSM (parameter) gradient arrays ##
      Phi_grad_param = Phi_grad_f(param_obs,param_sig,n,K,param_bias,n_obs)
      XBeta_grad_param = XBeta_grad_f(param_obs,param_sig,n,t,param_bias,n_obs)
      R_grad_param = R_grad_f(param_obs,param_sig,n,param_bias,n_obs)
      G_grad_param = G_grad_f(param_obs,param_sig,n,K,wave,cosine_indices,dt=dt,nu,n_cos,
                        param_bias)
      Q_grad_param = Q_grad_f(param_obs,param_sig,n,K,wave,n_cos,nu,dt=dt,norm=TRUE,
                        param_bias,spec_weights_mat = spec_weights_mat)
      
      
      ## LGSSM (space) gradient arrays ##
      Phi_grad_space = Phi_grad_space_f(wave,cosine_indices,n_cos,n,n_obs,coords = coords,grad_coord_index = 1:n_obs)
      XBeta_grad_space = XBeta_grad_space_f(n_obs,t)
      R_grad_space = R_grad_space_f(n_obs)
      G_grad_space = G_grad_space_f(n_obs,K)
      Q_grad_space = Q_grad_space_f(n_obs,K)
      
      
      ## LGSSM 2nd (parameter) derivative arrays ##
      if(!is.null(hessian_indices)){
        Phi_grad2_param = Phi_grad2_f(param_obs,param_sig,n,K,param_bias,n_obs)
        XBeta_grad2_param = XBeta_grad2_f(param_obs,param_sig,n,t,param_bias,n_obs)
        R_grad2_param = R_grad2_f(param_obs,param_sig,n,param_bias,n_obs)
        G_grad2_param = G_grad2_f(param_obs,param_sig,n,K,wave,cosine_indices,dt=dt,nu,
                            n_cos,param_bias)
        Q_grad2_param = Q_grad2_f(param_obs,param_sig,n,K,wave,n_cos,nu,dt=dt,norm=TRUE,
                            param_bias,spec_weights_mat = spec_weights_mat)
      }
      
      ## Current observation ##
      y_obs = matrix(y[i,],nrow=1)
      
      ## Linear predictor at current time ##
      XBeta_obs = matrix(XBeta[i,],nrow=1)
      XBeta_grad_param_obs = lapply(XBeta_grad_param,function(A) matrix(A[i,],nrow=1))
      XBeta_grad_space_obs = lapply(XBeta_grad_space,function(A) matrix(A[i,],nrow=1))
      if(!is.null(hessian_indices)){
        XBeta_grad2_param_obs = lapply(XBeta_grad2_param,function(A) matrix(A[i,],nrow=1))
      }
      
      
      ## Run TKF (parameters) for 1 time step ##
      
      if(length(grad_params_hessian)>0){
        
        ## Run MTKF2 (parameters) ##
        MTKF2_param = MultiTangentKalmanFilter2(y = y_obs, Phi = Phi, R = R, G = G, 
                                          Q = Q, Phi_grad = Phi_grad_param, 
                                          R_grad = R_grad_param, 
                                          G_grad = G_grad_param, 
                                          Q_grad = Q_grad_param, 
                                          Phi_grad2 = Phi_grad2_param, 
                                          R_grad2 = R_grad2_param, 
                                          G_grad2 = G_grad2_param, 
                                          Q_grad2 = Q_grad2_param, 
                                          mu0 = mu0, Sigma0 = Sigma0, 
                                          mu0_grad = mu0_grad_param, 
                                          Sigma0_grad = Sigma0_grad_param, 
                                          mu0_grad2 = mu0_grad2_param, 
                                          Sigma0_grad2 = Sigma0_grad2_param, 
                                          n = n_obs, t = 1, K = K, 
                                          log_lik = TRUE, n_param = n_param, 
                                          grad_params = grad_params_hessian,
                                          lin_pred = XBeta_obs,
                                          lin_pred_grad = XBeta_grad_param_obs,
                                          lin_pred_grad2 = XBeta_grad2_param_obs)
        
        ## Log likelihood ##
        ll[(j-1)*t+i] = unlist(MTKF2_param$ll)
        
        ## Log likelihood (parameter) gradient ##
        ll_grads_param[(j-1)*t+i,grad_params_hessian] = 
          unlist(MTKF2_param$ll_grad)[grad_params_hessian]
        
        ## Log likelihood second (parameter) gradient ##
        ll_grads2_param[(j-1)*t+i,grad_params_hessian] = 
          unlist(MTKF2_param$ll_grad2)[grad_params_hessian]
        
        ## Objective function ##
        obj_func[(j-1)*t+i] = objective_function(t=1,MTKF2_param$Sigma_i_i,
                                                 W_fourier)
        
        ## Objective function (parameter) gradient ##
        obj_func_grads_param[(j-1)*t+i,] = 
          objective_function_grad_param(n_param,t=1,MTKF2_param$Sigma_i_i_grad,
                                        W_fourier)
        
        ## KF means ##
        kf_means[(j-1)*t+i,] = MTKF2_param$mu_i_i[2,]
        
        MTKF_param = MTKF2_param
      }
      
      if(length(grad_params_n_hessian)>0){
        
        ## Run MTKF (parameters) ##
        MTKF_param = MultiTangentKalmanFilter(y = y_obs, Phi = Phi,R = R, G = G, 
                                        Q = Q, Phi_grad = Phi_grad_param, 
                                        R_grad = R_grad_param, 
                                        G_grad = G_grad_param, 
                                        Q_grad = Q_grad_param, mu0 = mu0, 
                                        Sigma0 = Sigma0, 
                                        mu0_grad = mu0_grad_param, 
                                        Sigma0_grad = Sigma0_grad_param, 
                                        n = n_obs, 
                                        t = 1, K = K, log_lik = TRUE, 
                                        n_param = n_param,
                                        grad_params = grad_params_n_hessian,
                                        lin_pred = XBeta_obs,
                                        lin_pred_grad = XBeta_grad_param_obs)
        
        ## Log likelihood ##
        ll[(j-1)*t+i] = unlist(MTKF_param$ll)
        
        ## Log likelihood (parameter) gradient ##
        ll_grads_param[(j-1)*t+i,grad_params_n_hessian] = 
          unlist(MTKF_param$ll_grad)[grad_params_n_hessian]
        
        ## Objective function ##
        obj_func[(j-1)*t+i] = objective_function(t=1,MTKF_param$Sigma_i_i,
                                                 W_fourier)
        
        ## Objective function (parameter) gradient ##
        obj_func_grads_param[(j-1)*t+i,] = 
          objective_function_grad_param(n_param,t=1,MTKF_param$Sigma_i_i_grad,
                                        W_fourier)
        
        ## KF means ##
        kf_means[(j-1)*t+i,] = MTKF_param$mu_i_i[2,]
      }
      
      
      ## Run TKF (space) for 1 time step ##
      
      ## Run MTKF (space) ##
      MTKF_space = MultiTangentKalmanFilter(y = y_obs, Phi = Phi,R = R, G = G, 
                                            Q = Q, Phi_grad = Phi_grad_space, 
                                            R_grad = R_grad_space, 
                                            G_grad = G_grad_space, 
                                            Q_grad = Q_grad_space, mu0 = mu0, 
                                            Sigma0 = Sigma0, 
                                            mu0_grad = mu0_grad_space, 
                                            Sigma0_grad = Sigma0_grad_space, 
                                            n = n_obs, t = 1, K = K, 
                                            log_lik = TRUE, n_param = 2*n_obs,
                                            grad_params = 1:(2*n_obs),
                                            lin_pred = XBeta_obs,
                                            lin_pred_grad = XBeta_grad_space_obs)
      
      ## Log likelihood (spatial) gradient ##
      for (k in 1:n_obs){
        ll_grads_space[[k]][(j-1)*t+i,] = unlist(MTKF_space$ll_grad[(2*k-1):(2*k)])
      }
      
      ## Objective function (spatial) gradient ##
      obj_func_grad = objective_function_grad_space(n_obs,t=1,MTKF_space$Sigma_i_i_grad,W_fourier)
      for (k in 1:n_obs){
        obj_func_grads_space[[k]][(j-1)*t+i,] = obj_func_grad[[k]]
      }
      
      ## Parameter update ##
      
      ## Rescale to log-scale ##
      param_est[(j-1)*t+i+1,] = param_est[(j-1)*t+i,]
      param_est[(j-1)*t+i+1,log_indices] = 
        log(param_est[(j-1)*t+i+1,log_indices])
      
      ## Gradient ascent ##
      param_est[(j-1)*t+i+1,grad_params_n_hessian] = 
        param_est[(j-1)*t+i+1,grad_params_n_hessian] + 
        step_size[(j-1)*t+i,grad_params_n_hessian]*
        ll_grads_param[(j-1)*t+i,grad_params_n_hessian]
      
      ## Gradient ascent (w/ Hessian) ##
      param_est[(j-1)*t+i+1,grad_params_hessian] =   
        param_est[(j-1)*t+i+1,grad_params_hessian] - 
        step_size[(j-1)*t+i,grad_params_hessian]*
        ll_grads2_param[(j-1)*t+i,grad_params_hessian]^(-1)*
        ll_grads_param[(j-1)*t+i,grad_params_hessian]
      
      ## Rescale to standard scale ##
      param_est[(j-1)*t+i+1,log_indices] = 
        exp(param_est[(j-1)*t+i+1,log_indices])
      
      ## Ensure alpha is in [0,pi/2] ##
      if(param_est[(j-1)*t+i+1,6]<0){
        param_est[(j-1)*t+i+1,6] = param_est[(j-1)*t+i,6]/2
      }
      if(param_est[(j-1)*t+i+1,6]>=pi/2){
        param_est[(j-1)*t+i+1,6] = (param_est[(j-1)*t+i,6]/2+pi/2)/2
      }
      
      ## Ensure mu_x is in [-0.5,0.5] ##
      if(param_est[(j-1)*t+i+1,7]<(-0.5)){
        param_est[(j-1)*t+i+1,7] = param_est[(j-1)*t+i,7]/2
      }
      if(param_est[(j-1)*t+i+1,7]>0.5){
        param_est[(j-1)*t+i+1,7] = param_est[(j-1)*t+i,7]/2
      }
      
      ## Ensure mu_y is in [-0.5,0.5] ##
      if(param_est[(j-1)*t+i+1,8]<(-0.5)){
        param_est[(j-1)*t+i+1,8] = param_est[(j-1)*t+i,8]/2
      }
      if(param_est[(j-1)*t+i+1,8]>0.5){
        param_est[(j-1)*t+i+1,8] = param_est[(j-1)*t+i,8]/2
      }
      
      
      ## Back Tracking Line Search ##
      if(back_track){
        
        ## Back Tracking Line Search For Each Parameter ##
        if(back_track_independent){
        
          checked = rep(FALSE,n_param)
          
          ## Loop over parameters ##
          for(l in back_track_param){
            
            while(!checked[l]){
              
              ## Parameter vector with only l^th parameter updated ##
              param_check = param_est[(j-1)*t+i,]
              param_check[l] = param_est[(j-1)*t+i+1,l]
              
              ## Observation, signal, bias parameters ##
              param_obs_new = param_obs_func(param_obs0$m,param_obs0$m_indices,
                                             param_check[param_obs_indices])
              param_sig_new = param_check[param_sig_indices]
            
              if(n_param_bias>0){
                param_bias_new = param_bias_func(param_bias0$m,
                                                 param_bias0$m_indices,
                                                 param_check[param_bias_indices])
              }
              else{
                param_bias_new = param_bias_func()
              }
            
              ## LGGSM matrices for this parameter vector ##
              Phi_new = Phi_f(wave,cosine_indices,n_cos,n,n_obs,coords)
              XBeta_new = XBeta_f(param_bias_new,n,t,n_obs)
              R_new = R_f(param_obs_new,n,n_obs)
              G_new = G_f(param_sig_new,wave,cosine_indices,dt=dt,n_cos=n_cos)
              Q_new = Q_f(param_sig_new,wave,n,n_cos=n_cos,nu=1,dt=dt,norm=TRUE,
                          spec_weights_mat = spec_weights_mat)
              
              ## Compute likelihood for this parameter vector ##
              ll_update = KalmanFilter(y = y_obs, Phi = Phi_new, R = R_new, 
                                       G = G_new, Q = Q_new, m0 = mu0, 
                                       R0 = Sigma0, linear_pred = XBeta_obs, 
                                       N = n_obs, t = 1, K = K, 
                                       log_lik = TRUE)$ll
              
              ## If new likelihood > old likelihood, accept new parameter ##
              if(ll_update > ll[(j-1)*t+i]){
                checked[l] = TRUE
              }
              
              checked[l] = TRUE
              
              ## If old likelihood >~ new likelihood ##
              if((ll[(j-1)*t+i]-ll_update)>back_track_tol){
                
                ## Half future step-sizes ##
                if(!back_track_time_independent){
                  step_size[((j-1)*t+i+1):(n_iterations*t),l] = 
                    back_track_scale*step_size[((j-1)*t+i+1):(n_iterations*t),l]
                }
                
                ## Re-update new parameter ##
                param_est[(j-1)*t+i+1,l] = 
                  back_track_scale*((1/back_track_scale-1)*param_est[(j-1)*t+i,l] + param_est[(j-1)*t+i+1,l])
                
                ## Check again ##
                checked[l] = FALSE
              }
            }
          }
        }
        
        ## Back Tracking Line Search For All Parameters ##
        if(!back_track_independent){
          
          checked = FALSE
          
          while(!checked){
            
            ## Parameter vector with all elements updated ##
            param_check = param_est[(j-1)*t+i+1,]
            
            ## Observation, signal, bias parameters ##
            param_obs_new = param_obs_func(param_obs0$m,param_obs0$m_indices,
                                           param_check[param_obs_indices])
            param_sig_new = param_check[param_sig_indices]
            
            if(n_param_bias>0){
              param_bias_new = param_bias_func(param_bias0$m,
                                               param_bias0$m_indices,
                                               param_check[param_bias_indices])
            }
            else{
              param_bias_new = param_bias_func()
            }
            
            ## LGGSM matrices for this parameter vector ##
            Phi_new = Phi_f(wave,cosine_indices,n_cos,n,n_obs,coords)
            XBeta_new = XBeta_f(param_bias_new,n,t,n_obs)
            R_new = R_f(param_obs_new,n,n_obs)
            G_new = G_f(param_sig_new,wave,cosine_indices,dt=dt,n_cos=n_cos)
            Q_new = Q_f(param_sig_new,wave,n,n_cos=n_cos,nu=1,dt=dt,norm=TRUE,
                        spec_weights_mat = spec_weights_mat)
            
            ## Compute likelihood for this parameter vector ##
            ll_update = KalmanFilter(y = y_obs, Phi = Phi_new, R = R_new, 
                                     G = G_new, Q = Q_new, m0 = mu0, 
                                     R0 = Sigma0, linear_pred = XBeta_obs, 
                                     N = n_obs, t = 1, K = K, log_lik = TRUE)$ll
            
            ## If new likelihood > old likelihood, accept new parameter ##
            if(ll_update > ll[(j-1)*t+i]){
              checked = TRUE
            }
            
            checked = TRUE
            
            ## If old likelihood >~ new likelihood ## 
            if((ll[(j-1)*t+i]-ll_update)>back_track_tol){
              
              ## Half step-size ##
              if(!back_track_time_independent){
                step_size[((j-1)*t+i+1):(n_iterations*t),] = 
                  back_track_scale*step_size[((j-1)*t+i+1):(n_iterations*t),]
              }
              
              ## Re-update new parameter ##
              param_est[(j-1)*t+i+1,] = 
                back_track_scale*((1/back_track_scale-1)*param_est[(j-1)*t+i,] + param_est[(j-1)*t+i+1,])
              
              ## Check again ##
              checked = FALSE
            }
            
            
          }
          
          
        }
      }
      
      
      ## Update 'initial' values for TKF ##
      mu0 = MTKF_param$mu_i_i[2,]
      Sigma0 = MTKF_param$Sigma_i_i[2,,]
      mu0_grad_param = lapply(MTKF_param$mu_i_i_grad, function(A) A[2,])
      Sigma0_grad_param = lapply(MTKF_param$Sigma_i_i_grad, function(A) A[2,,])
      mu0_grad_space = lapply(MTKF_space$mu_i_i_grad, function(A) A[2,])
      Sigma0_grad_space = lapply(MTKF_space$Sigma_i_i_grad, function(A) A[2,,])
      if(length(grad_params_hessian)>0){
        mu0_grad2_param = lapply(MTKF2_param$mu_i_i_grad2, function(A) A[2,])
        Sigma0_grad2_param = lapply(MTKF2_param$Sigma_i_i_grad2, function(A) A[2,,])
      }
      
      ## Print iteration ##
      if(!is.null(print_iter_freq)){
        if(iter_numbers[(j-1)*t+i]%%print_iter_freq==0){
          cat("Iteration ", iter_numbers[(j-1)*t+i], "\n")
        }
      }
      
      ## Mean parameter estimates ##
      if(mean_param){
        
        if(!is.null(param_est_starts)){
          
          for (k in 1:n_param){
            
            if(((j-1)*t+i)>(param_est_starts[[k]][1]+1)){
              
              ## Which start value to use to compute mean parameter estimate ##
              which_start[[k]] = max(which(((j-1)*t+i)>(param_est_starts[[k]]+1)))
              
              ## Which start value to use to plot mean parameter estimate ##
              if (k %in% param_sig_indices){
                mean_param_est_plot_start[[k]][which_start[[k]]] = param_sig_t_dep$t_vals[[max(which(unlist(param_sig_true$t_vals)<param_est_starts[[k]][which_start[[k]]]))]]
              }
              if (k %in% param_obs_indices){
                mean_param_est_plot_start[[k]][which_start[[k]]] = param_obs_t_dep$t_vals[[max(which(unlist(param_obs_true$t_vals)<param_est_starts[[k]][which_start[[k]]]))]]
              }
              if(n_param_bias>0){
                if (k %in% param_bias_indices){
                  mean_param_est_plot_start[[k]][which_start[[k]]] = param_bias_t_dep$t_vals[[max(which(unlist(param_bias_true$t_vals)<param_est_starts[[k]][which_start[[k]]]))]]
                }
              }
              
              ## Compute mean parameter estimate indices ##
              mean_param_est_indices[[k]] = 
                (param_est_starts[[k]][which_start[[k]]]+1):(min(((j-1)*t+i),param_est_ends[[k]][which_start[[k]]]))
              
              ## Compute mean parameter estimate ##
              mean_param_est[[k]][which_start[[k]]] = 
                mean(param_est[mean_param_est_indices[[k]],k])
            }
          }
        }
      }
      
      ## Plots ##
      if(plot){
        
        if(((j-1)*t+i)%%plot_freq==0){
          
          ## Open plotting device ##
          if(save_plot){
            if(((j-1)*t+i)==(n_iterations*t)){
              setwd(fig_wd)
              pdf(plot_filename,width=10,height=3.5)
              par(mar=c(2,4.5,1.5,1.5),mgp=c(1.5,.5,0))
            }
          }
          
          ## Number of Plots ##
          if(n_param_bias==0){par(mfrow=c(3,3))}
          if(n_param_bias>0){par(mfrow=c(3,4))}
          
          ## Graphics parameters ##
          par(mar=c(3,4.5,1.5,1.5),mgp=c(1.5,.5,0))
          
          ## Plot Points ##
          plot_points = seq(1,(j-1)*t+i,plot_point_freq)
          
          ## Signal Parameters ##
          for (k in param_sig_indices){
            
            y_lim = NULL
            y_lab = expression(theta)
            main = NULL
            
            if(!is.null(plot_limits)){
              y_lim = plot_limits[[k]]
            }

            if(!is.null(param_sig_names)){
              y_lab = parse(text=paste("hat(",toString(param_sig_names[[k]]),")",sep=""))
              main = param_sig_names[[k]]
            }
            
            plot(plot_points, param_est[plot_points,k], cex=0.4,
                 xlim = c(1,n_iterations*t+1), xlab = "t",ylab = y_lab,
                 ylim = y_lim, main = main)
            
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
              if(((j-1)*t+i)>(param_est_starts[[k]][1]+1)){
                for (m in 1:(which_start[[k]])){
                  if(plot_vertical){
                    abline(v=param_est_starts[[k]][m],col="black",lty=2)
                  }
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
                y_lab = parse(text=paste("hat(",toString(param_obs_names[[k-8]]),")",sep=""))
                main = param_obs_names[[k-8]]
              }
              
              plot(plot_points, param_est[plot_points,k], cex=0.4,
                   xlim = c(1,n_iterations*t+1) ,xlab="t", ylab=expression(theta),
                   ylim = y_lim, main = main, pch=19)
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
              
            if(mean_param){
              if(((j-1)*t+i)>(param_est_starts[[k]][1]+1)){
                for (m in 1:(which_start[[k]])){
                  if(plot_vertical){
                    abline(v=param_est_starts[[k]][m],col="black",lty=2)
                  }
                  if(length(param_obs_true$t_vals)==1){
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
          
          ## Bias Parameters ##
          if(n_param_bias>0){
            for (k in param_bias_indices){
              if(k==param_bias_indices[[1]]){
                
                y_lim = NULL
                y_lab = expression(theta)
                main = NULL
                
                if(!is.null(plot_limits)){
                  y_lim = plot_limits[[10]]
                }

                if(!is.null(param_bias_names)){
                  y_lab = parse(text=paste("hat(",toString(param_obs_names[[k-(n_param_sig+n_param_obs)]]),")",sep=""))
                  main = param_bias_names[[k-(n_param_sig+n_param_obs)]]
                }
                
                plot(plot_points, param_est[plot_points,k], cex=0.4,
                     xlim = c(1,n_iterations*t+1) ,xlab="t", ylab=expression(theta),
                     ylim = y_lim, main = main, pch=19)
              }
              
              if(k>param_bias_indices[[1]]){
                points(plot_points,param_est[plot_points,k],pch=19,cex=0.4)
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
                
              if(mean_param){
                if(((j-1)*t+i)>(param_est_starts[[k]][1]+1)){
                  for (m in 1:(which_start[[k]])){
                    if(plot_vertical){
                      abline(v=param_est_starts[[k]][m],col="black",lty=2)
                    }
                    if(length(param_bias_true$t_vals)==1){
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
          }
          
          ## Close plotting device ##
          if(save_plot){
            if(((j-1)*t+i)==(n_iterations*t)){
              dev.off()
            }
          }
        }
      }
    }
  }
  
  ## Ouput ##
  output = list(param_est = param_est, ll = ll, ll_grads_param = ll_grads_param,
                ll_grads_space = ll_grads_space, 
                obj_func = obj_func,
                obj_func_grads_param = obj_func_grads_param, 
                obj_func_grads_space = obj_func_grads_space,
                kf_means = kf_means)
  
  if(!is.null(grad_params_hessian)){
    output = c(output,list(ll_grads2_param=ll_grads2_param))
  }
  
  if(mean_param){
    output = c(output,list(mean_param_est = mean_param_est))
  }
  
  output = c(output,list(mu0 = mu0, Sigma0 = Sigma0, 
                         mu0_grad_param = mu0_grad_param,
                         Sigma0_grad_param = Sigma0_grad_param, 
                         mu0_grad_space = mu0_grad_space,
                         Sigma0_grad_space = Sigma0_grad_space))
  
  if(!is.null(grad_params_hessian)){
    output = c(output,list(mu0_grad2_param = mu0_grad2_param, 
                           Sigma0_grad2_param = Sigma0_grad2_param))
  }
  
  if(back_track){
    output = c(output,list(step_size = step_size))
  }
  
  ## Save ##
  if(save){
    saveRDS(output, file = filename)
  }
  
  ## Output ##
  return(output)
}


################################

#########################################
# SPECTRAL RECURSIVE MAXIMUM LIKELIHOOD #
#########################################

RML_spectral <- function(w, n, t, K = n^2, x = NULL, param_obs0, param_sig0, 
                         n_param, grad_params = 1:9, step_size,
                         log_scale = FALSE, log_indices,
                         hessian_indices = NULL, n_iterations = 1,
                         plot = FALSE, plot_limits = NULL, 
                         param_obs_true = NULL, param_sig_true = NULL, 
                         param_obs_names = NULL, param_sig_names = NULL, 
                         param_est_starts, param_est_ends, 
                         plot_point_freq = 1, plot_freq = 10, save_plot = FALSE,
                         plot_filename = "RML_plot.pdf", save = FALSE,
                         filename = "RML.pdf", mean_param = FALSE,
                         print_iter_freq = 100, 
                         iter_numbers = 1:(n_iterations*t),
                         m0 = NULL, R0 = NULL, m0_grad_param = NULL, 
                         R0_grad_param = NULL, m0_grad2_param = NULL, 
                         R0_grad2_param = NULL,
                         back_track = FALSE,back_track_tol = 1e-3, 
                         back_track_independent = TRUE,
                         back_track_time_independent = FALSE, 
                         back_track_scale = 0.5,
                         back_track_param = 1:n_param,
                         plot_vertical = TRUE,
                         W_fourier = NULL, dt = 1){
  
  ## Number of parameters ##
  n_param_sig = length(param_sig0)
  n_param_obs = param_obs0$m
  
  ## Parameter indices ##
  param_sig_indices = 1:n_param_sig
  param_obs_indices = (n_param_sig+1):(n_param_sig+n_param_obs)
  
  ## Matrix of parameter estimates ##
  param_est = matrix(0,ncol=n_param,nrow=t*n_iterations+1)
  
  ## Mean parameter estimates ##
  if(mean_param){
    
    ## Indices to compute mean parameter estimate(s) ##
    mean_param_est_indices = rep(list(c()),n_param)
    
    ## Mean parameter estimate(s) ##
    mean_param_est = rep(list(rep(0,100)),n_param)
    
    ## Start indices for mean parameter estimate(s) ##
    which_start = rep(list(c()),n_param)
    
    ## Start(s) for mean parameter estimate plot(s) ##
    mean_param_est_plot_start = rep(list(rep(0,100)),n_param)
    
  }
  
  ## Initial parameter estimate ##
  param0 = c(param_sig0,as.numeric(unlist(param_obs0$tau2_vals)))
  
  ## Rescale if param0 on log scale ##
  if(log_scale){param0[log_indices] <- exp(param0[log_indices])}
  
  ## Add initial parameter estimate to matrix of parameter estimates ##
  param_est[1,] = param0
  
  ## Matrices of conditional log-lik, log-lik gradients, 2nd derivatives ##
  ll = rep(0,t*n_iterations)
  ll_grads_param = matrix(0,ncol=n_param,nrow=t*n_iterations)
  ll_grads_space = rep(list(matrix(0,ncol=2,nrow=t*n_iterations)),n_obs)
  ll_grads2_param = matrix(0,ncol=n_param,nrow=t*n_iterations)
  
  ## Matrices of obj. func., conditional obj. func. grads ##
  obj_func = matrix(0,ncol=1,nrow=t*n_iterations)
  obj_func_grads_param = matrix(0,ncol=n_param,nrow=t*n_iterations)
  obj_func_grads_space = rep(list(matrix(0,ncol=2,nrow=t*n_iterations)),n_obs)
  
  ## Hessian & Non-Hessian parameters to estimate ##
  grad_params_hessian = intersect(hessian_indices,grad_params)
  grad_params_n_hessian = c(setdiff(hessian_indices,grad_params),
                            setdiff(grad_params,hessian_indices))
  
  ## SPDE_FT object ##
  SPDE_FT = spde_initialise(n,t,K)
  wave = SPDE_FT$wave
  cosine_indices = SPDE_FT$cosine_indices
  n_cos = SPDE_FT$n_cos
  n_cos_sin = length(cosine_indices)
  K = SPDE_FT$K
  
  ## Initial values ##
  if(is.null(m0)){
    m0 = rep(0,K)
  }
  if(is.null(R0)){
    R0 = Q_f_vec(param_sig0,wave=wave,n=n,n_cos=n_cos,nu=1,dt=dt,norm=TRUE)
  }
  if(is.null(m0_grad_param)){
    m0_grad_param = rep(list(rep(0,K)),n_param)
  }
  if(is.null(R0_grad_param)){
    R0_grad_param = Q_grad_f_vec(param_obs0,param_sig0,n=n,K=K,wave=wave,n_cos=n_cos,
                           nu=1,dt=dt,norm=TRUE)
  }
  if(is.null(m0_grad2_param)){
    m0_grad2_param = rep(list(rep(0,K)),n_param)
  }
  if(is.null(R0_grad2_param)){
    R0_grad2_param = Q_grad2_f_vec(param_obs0,param_sig0,n=n,K=K,wave=wave,
                             n_cos=n_cos,nu=1,dt=dt,norm=TRUE)
  }
  
  ### RML ###
  for (j in 1:n_iterations){
    
    for (i in 1:t){
      
      ## Current obs & sig parameters ##
      param_obs = param_obs_func(param_obs0$m,param_obs0$m_indices,
                                 param_est[(j-1)*t+i,param_obs_indices])
      param_sig = param_est[(j-1)*t+i,param_sig_indices]
      
      ## Current observation ##
      w_obs = matrix(w[i,],nrow=1)
      
      if(length(grad_params_hessian)>0){
        
        ## Run MTKF2 (parameters) ##
        MTKF2_param = Multi_TKF2_spectral(w = w_obs,param_obs = param_obs,
                                    param_sig = param_sig, n = n, t = 1, K = K,
                                    log_lik = TRUE, dt = dt, m0 = m0, R0 = R0,
                                    m0_grad = m0_grad_param, 
                                    R0_grad = R0_grad_param,
                                    m0_grad2 = m0_grad2_param, 
                                    R0_grad2 = R0_grad2_param,
                                    grad_params = grad_params_hessian,
                                    log_lik_grad = TRUE, log_lik_grad2 = TRUE)
        
        ## Log Likelihood ##
        ll[(j-1)*t+i] = unlist(MTKF2_param$ll)
        
        ## Log Likelihood (Parameter) Gradient ##
        ll_grads_param[(j-1)*t+i,grad_params_hessian] = 
          unlist(MTKF2_param$ll_grad)[grad_params_hessian]
        
        ## Log Likelihood 2nd (Parameter) Gradient ##
        ll_grads2_param[(j-1)*t+i,grad_params_hessian] = 
          unlist(MTKF2_param$ll_grad2)[grad_params_hessian]
        
        ## Objective Function ##
        obj_func[(j-1)*t+i] = 
          objective_function(t=1,diag(MTKF2_param$Sigma_i_i[2,]),W_fourier)
        
        ## Objective Function (Parameter) Gradient ##
        obj_func_grads_param[(j-1)*t+i,] = 
          objective_function_grad_param(n_param,t=1,
                                        lapply(MTKF2_param$Sigma_i_i_grad, function(A) diag(A[2,])),
                                        W_fourier)
        
        MTKF_param = MTKF2_param
      }
      
      if(length(grad_params_n_hessian)>0){
        
        ## Run MTKF (parameters) ##
        MTKF_param = Multi_TKF_spectral(w=matrix(w[i,],nrow=1),
                                  param_obs=param_obs,param_sig=param_sig,
                                  n=n,t=1,K=K,log_lik=TRUE,dt=dt,m0=m0,R0=R0,
                                  m0_grad=m0_grad_param,R0_grad=R0_grad_param,
                                  grad_params=grad_params_n_hessian,
                                  log_lik_grad=TRUE)
        
        ## Log Likelihood ##
        ll[(j-1)*t+i] = unlist(MTKF_param$ll)
        
        ## Log Likelihood Gradient ##
        ll_grads_param[(j-1)*t+i,grad_params_n_hessian] = 
          unlist(MTKF_param$ll_grad)[grad_params_n_hessian]
        
        ## Objective Function ##
        obj_func[(j-1)*t+i] = 
          objective_function(t=1,diag(MTKF_param$Sigma_i_i[2,]),W_fourier)
        
        ## Objective Function (Parameter) Gradient ##
        obj_func_grads_param[(j-1)*t+i,] = 
          objective_function_grad_param(n_param,t=1,
                                        lapply(MTKF_param$Sigma_i_i_grad, function(A) diag(A[2,])),
                                        W_fourier)
      }
      
      
      ## Parameter update ##
      
      ## Rescale to log-scale ##
      param_est[(j-1)*t+i+1,] = param_est[(j-1)*t+i,]
      param_est[(j-1)*t+i+1,log_indices] = 
        log(param_est[(j-1)*t+i+1,log_indices])
      
      ## Gradient ascent  ##
      param_est[(j-1)*t+i+1,grad_params_n_hessian] = 
        param_est[(j-1)*t+i+1,grad_params_n_hessian] + 
        step_size[(j-1)*t+i,grad_params_n_hessian]*
        ll_grads_param[(j-1)*t+i,grad_params_n_hessian]
      
      ## Gradient ascent (w/ Hessian) ##
      param_est[(j-1)*t+i+1,grad_params_hessian] =   
        param_est[(j-1)*t+i+1,grad_params_hessian] - 
        step_size[(j-1)*t+i,grad_params_hessian]*
        ll_grads2_param[(j-1)*t+i,grad_params_hessian]^(-1)*
        ll_grads_param[(j-1)*t+i,grad_params_hessian]
      
      ## Rescale to standard scale ##
      param_est[(j-1)*t+i+1,log_indices] = 
        exp(param_est[(j-1)*t+i+1,log_indices])
      
      ## Ensure alpha is in [0,pi/2] ##
      if(param_est[(j-1)*t+i+1,6]<0){
        param_est[(j-1)*t+i+1,6] = param_est[(j-1)*t+i,6]/2
      }
      if(param_est[(j-1)*t+i+1,6]>=pi/2){
        param_est[(j-1)*t+i+1,6] = (param_est[(j-1)*t+i,6]/2+pi/2)/2
      }
      
      ## Ensure mu_x is in [-0.5,0.5]
      if(param_est[(j-1)*t+i+1,7]<(-0.5)){
        param_est[(j-1)*t+i+1,7] = param_est[(j-1)*t+i,7]
      }
      if(param_est[(j-1)*t+i+1,7]>0.5){
        param_est[(j-1)*t+i+1,7] = param_est[(j-1)*t+i,7]
      }
      
      ## Ensure mu_y is in [-0.5,0.5]
      if(param_est[(j-1)*t+i+1,8]<(-0.5)){
        param_est[(j-1)*t+i+1,8] = param_est[(j-1)*t+i,8]
      }
      if(param_est[(j-1)*t+i+1,8]>0.5){
        param_est[(j-1)*t+i+1,8] = param_est[(j-1)*t+i,8]
      }
      
      ## Print iteration ##
      if(!is.null(print_iter_freq)){
        if(iter_numbers[(j-1)*t+i]%%print_iter_freq==0){
          cat("Iteration ", iter_numbers[(j-1)*t+i], "\n")
        }
      }
      
      
      
      
      ## Back Tracking Line Search ##
      if(back_track){
        
        ## Back Tracking Line Search For Each Parameter ##
        if(back_track_independent){
          
          checked = rep(FALSE,n_param)
          
          ## Loop over parameters ##
          for(l in back_track_param){
            
            while(!checked[l]){
              
              ## Parameter vector with only l^th parameter updated ##
              param_check = param_est[(j-1)*t+i,]
              param_check[l] = param_est[(j-1)*t+i+1,l]
              
              ## Observation, signal parameters ##
              param_obs_new = param_obs_func(param_obs0$m,param_obs0$m_indices,
                                             param_check[param_obs_indices])
              param_sig_new = param_check[param_sig_indices]
              
              ## Compute likelihood for this parameter vector ##
              ll_update = KF_spectral(w = w_obs, param_obs = param_obs_new,
                                      param_sig = param_sig_new, n = n, t = 1, 
                                      K = K, log_lik = TRUE, dt = dt, m0 = m0, 
                                      R0 = R0)$ll
              
              ## If new likelihood > old likelihood, accept new parameter ##
              if(ll_update > ll[(j-1)*t+i]){
                checked[l] = TRUE
              }
              
              checked[l] = TRUE
              
              ## If old likelihood >~ new likelihood ##
              if((ll[(j-1)*t+i]-ll_update)>back_track_tol){
                
                ## Half future step-sizes ##
                if(!back_track_time_independent){
                  step_size[((j-1)*t+i+1):(n_iterations*t),l] = 
                    back_track_scale*step_size[((j-1)*t+i+1):(n_iterations*t),l]
                }
                
                ## Re-update new parameter ##
                param_est[(j-1)*t+i+1,l] = 
                  back_track_scale*((1/back_track_scale-1)*param_est[(j-1)*t+i,l] + param_est[(j-1)*t+i+1,l])
                
                ## Check again ##
                checked[l] = FALSE
              }
            }
          }
        }
        
        ## Back Tracking Line Search For All Parameters ##
        if(!back_track_independent){
          
          checked = FALSE
          
          while(!checked){
            
            ## Parameter vector with all elements updated ##
            param_check = param_est[(j-1)*t+i+1,]
            
            ## Observation, signal parameters ##
            param_obs_new = param_obs_func(param_obs0$m,param_obs0$m_indices,
                                           param_check[param_obs_indices])
            param_sig_new = param_check[param_sig_indices]
            
            ## Compute likelihood for this parameter vector ##
            ll_update = KF_spectral(w = w_obs, param_obs = param_obs_new,
                                    param_sig = param_sig_new, n = n, t = 1, 
                                    K = K, log_lik = TRUE, dt = dt, m0 = m0, 
                                    R0 = R0)$ll
            
            ## If new likelihood > old likelihood, accept new parameter ##
            if(ll_update > ll[(j-1)*t+i]){
              checked = TRUE
            }
            
            checked = TRUE
            
            ## If old likelihood >~ new likelihood ## 
            if((ll[(j-1)*t+i]-ll_update)>back_track_tol){
              
              ## Half step-size ##
              if(!back_track_time_independent){
                step_size[((j-1)*t+i+1):(n_iterations*t),] = 
                  back_track_scale*step_size[((j-1)*t+i+1):(n_iterations*t),]
              }
              
              ## Re-update new parameter ##
              param_est[(j-1)*t+i+1,] = 
                back_track_scale*((1/back_track_scale-1)*param_est[(j-1)*t+i,] + param_est[(j-1)*t+i+1,])
              
              ## Check again ##
              checked = FALSE
            }
          }
        }
      }
      
      
      
      ## Update 'initial' values for TKF ##
      m0 = MTKF_param$mu_i_i[2,]
      R0 = MTKF_param$Sigma_i_i[2,]
      m0_grad_param = lapply(MTKF_param$mu_i_i_grad, function(A) A[2,])
      R0_grad_param = lapply(MTKF_param$Sigma_i_i_grad, function(A) A[2,])
      if(length(grad_params_hessian)>0){
        m0_grad2_param = lapply(MTKF2_param$mu_i_i_grad2, function(A) A[2,])
        R0_grad2_param = lapply(MTKF2_param$Sigma_i_i_grad2, function(A) A[2,])
      }
      
      ## Mean parameter estimates ##
      if(mean_param){
        
        if(!is.null(param_est_starts)){
          
          for (k in 1:n_param){
            
            if(((j-1)*t+i)>(param_est_starts[[k]][1]+1)){
              
              ## Which start value to use to compute mean parameter estimate ##
              which_start[[k]] = max(which(((j-1)*t+i)>(param_est_starts[[k]]+1)))
              
              ## Which start value to use to plot mean parameter estimate ##
              if (k %in% param_sig_indices){
                mean_param_est_plot_start[[k]][which_start[[k]]] = param_sig_t_dep$t_vals[[max(which(unlist(param_sig_true$t_vals)<param_est_starts[[k]][which_start[[k]]]))]]
              }
              if (k %in% param_obs_indices){
                mean_param_est_plot_start[[k]][which_start[[k]]] = param_obs_t_dep$t_vals[[max(which(unlist(param_obs_true$t_vals)<param_est_starts[[k]][which_start[[k]]]))]]
              }
              
              ## Compute mean parameter estimate indices ##
              mean_param_est_indices[[k]] = 
                (param_est_starts[[k]][which_start[[k]]]+1):(min(((j-1)*t+i),param_est_ends[[k]][which_start[[k]]]))
              
              ## Compute mean parameter estimate ##
              mean_param_est[[k]][which_start[[k]]] = 
                mean(param_est[mean_param_est_indices[[k]],k])
            }
          }
        }
      }
      
      ## Plots ##
      if(plot){
        
        if(((j-1)*t+i)%%plot_freq==0){
          
          ## Open plotting device ##
          if(save_plot){
            if(((j-1)*t+i)==(n_iterations*t)){
              setwd(fig_wd)
              pdf(plot_filename,width=10,height=6.5)
              par(mar=c(2,4.5,1.5,1.5),mgp=c(1.5,.5,0))
            }
          }
          
          ## Number of Plots ##
          par(mfrow=c(3,3))
          
          ## Graphics parameters ##
          par(mar=c(3,4.5,1.5,1.5),mgp=c(1.5,.5,0))
          
          ## Plot Points ##
          plot_points = seq(1,(j-1)*t+i,plot_point_freq)
          
          ## Signal Parameters ##
          for (k in param_sig_indices){
            
            y_lim = NULL
            y_lab = expression(theta)
            main = NULL
            
            if(!is.null(plot_limits)){
              y_lim = plot_limits[[k]]
            }
            
            if(!is.null(param_sig_names)){
              y_lab = parse(text=paste("hat(",toString(param_sig_names[[k]]),")",sep=""))
              main = param_sig_names[[k]]
            }
            
            plot(plot_points, param_est[plot_points,k], cex=0.4,
                 xlim = c(1,n_iterations*t+1), xlab = "t",ylab = y_lab,
                 ylim = y_lim, main = main)
            
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
              if(((j-1)*t+i)>(param_est_starts[[k]][1]+1)){
                for (m in 1:(which_start[[k]])){
                  if(plot_vertical){
                    abline(v=param_est_starts[[k]][m],col="black",lty=2)
                  }
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
                y_lab = parse(text=paste("hat(",toString(param_obs_names[[k-8]]),")",sep=""))
                main = param_obs_names[[k-8]]
              }
              
              plot(plot_points, param_est[plot_points,k], cex=0.4,
                   xlim = c(1,n_iterations*t+1) ,xlab="t", ylab=expression(theta),
                   ylim = y_lim, main = main)
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
            
            if(mean_param){
              if(((j-1)*t+i)>(param_est_starts[[k]][1]+1)){
                for (m in 1:(which_start[[k]])){
                  if(plot_vertical){
                    abline(v=param_est_starts[[k]][m],col="black",lty=2)
                  }
                  if(length(param_obs_true$t_vals)==1){
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
          
          ## Close plotting device ##
          if(save_plot){
            if(((j-1)*t+i)==(n_iterations*t)){
              dev.off()
            }
          }
        }
      }
      
    }
  }
  
  ## Ouput ##
  output = list(param_est=param_est,ll_grads_param=ll_grads_param)
  
  if(!is.null(grad_params_hessian)){
    output = c(output,list(ll_grads2_param = ll_grads2_param))
  }
  
  output = c(output,list(mu0 = m0, Sigma0 = R0, mu0_grad_param = m0_grad_param,
                         Sigma0_grad_param = R0_grad_param))
  
  if(!is.null(grad_params_hessian)){
    output = c(output,list(mu0_grad2_param = m0_grad2_param, 
                           Sigma0_grad2_param = R0_grad2_param))
  }
  
  if(back_track){
    output = c(output,list(step_size = step_size))
  }
  
  ## Save ##
  if(save){
    saveRDS(output, file = filename)
  }
  
  ## Output ##
  return(output)
  
}


#########################################

# ----

######################
# OBJECTIVE FUNCTION #
######################

objective_function <- function(t,Sigma_i_i,W_fourier){
  
  ## Sum of weighted trace of Sigma_i_i, i=1:t
  obj_func = sum(apply(Sigma_i_i,1,function(A) sum(diag(W_fourier%*%A)))[2:(t+1)])
  
  ## Output
  return(obj_func)
}

######################

###########################
# OBJECTIVE FUNCTION GRAD #
###########################

objective_function_grad_space <- function(n_obs,t,Sigma_i_i_grad,W_fourier){
  
  ## Gradient of objective function w.r.t (x,y) coords of all obs. locations ##
  obj_func_grad_list <- lapply(1:(2*n_obs),function(k) sum(apply(Sigma_i_i_grad[[k]],1,function(A) sum(diag(W_fourier%*%A)))[2:(t+1)]))

  ## Combine gradient w.r.t (x,y) for each observation location ##
  obj_func_grad = rep(list(c(0,0)),n_obs)
  for (i in 1:n_obs){
    obj_func_grad[[i]] = c(obj_func_grad_list[[2*i-1]],
                           obj_func_grad_list[[2*i]])
  }
  
  ## Output
  return(obj_func_grad)
}

objective_function_grad_param <- function(n_param,t,Sigma_i_i_grad,W_fourier){
  
  ## Gradient of objective function w.r.t (x,y) coords of all obs. locations ##
  obj_func_grad_list <- lapply(1:n_param,function(k) sum(apply(Sigma_i_i_grad[[k]],1,function(A) sum(diag(W_fourier%*%A)))[2:(t+1)]))
  
  ## Convert list into vector ##
  obj_func_grad = unlist(obj_func_grad_list)
  
  ## Output
  return(obj_func_grad)
}

###########################

####################################
# OFFLINE OPTIMAL SENSOR PLACEMENT #
####################################

BatchOptSens <- function(y, n_obs = dim(y)[2], t = dim(y)[1], K,
                         param_obs, param_sig, param_bias = param_bias_func(),
                         obs0, n_iterations, step_size, W_fourier, 
                         spec_weights_mat = NULL, plot = FALSE, W_coords = NULL, 
                         grad_obs = 1:n_obs, plot2d = FALSE,
                         plot_point_freq = 1, plot_freq = 50,
                         save_plot = FALSE, 
                         filename_plot = "BatchOptSens_plot.pdf",
                         save_plot2d = FALSE,
                         filename_plot2d = "BatchOptSens_2dplot.pdf",
                         grid = TRUE, print_iter_freq = 10,
                         save = FALSE, filename = "BatchOptSens",
                         iter_numbers = 1:n_iterations,
                         mu0 = NULL, Sigma0 = NULL, mu0_grad_param = NULL,
                         Sigma0_grad_param = NULL, mu0_grad_space = NULL,
                         Sigma0_grad_space = NULL, back_track = FALSE,
                         back_track_tol = 1e-3, back_track_independent = TRUE,
                         back_track_time_independent = FALSE, 
                         back_track_scale = 0.5,
                         back_track_obs = 1:n_obs,
                         plot_initial = FALSE,
                         leg = FALSE, dt = 1){
  
  ## Number of parameters ##
  n_param = length(param_sig0) + param_obs0$m + param_bias0$m
  
  ## Matrices of observation locations ##
  obs = rep(list(matrix(0,nrow=n_iterations+1,ncol=2)),n_obs)
  for (i in 1:n_obs){
    obs[[i]][1,] = obs0[[i]]
  }
  
  ## Matrices of log-lik,log-lik gradients, 2nd derivatives ##
  ll = rep(0,n_iterations)
  ll_grads_param = matrix(0,ncol=n_param,nrow=n_iterations)
  ll_grads_space = rep(list(matrix(0,ncol=2,nrow=n_iterations)),n_obs)
  ll_grads2_param = matrix(0,ncol=n_param,nrow=n_iterations)
  
  ## Matrices of obj. func., conditional obj. func. grads ##
  obj_func = matrix(0,ncol=1,nrow=n_iterations)
  obj_func_grads_param = matrix(0,ncol=n_param,nrow=n_iterations)
  obj_func_grads_space = rep(list(matrix(0,ncol=2,nrow=n_iterations)),n_obs)
  
  ## SPDE_FT object ##
  SPDE_FT = spde_initialise(n,t,K)
  wave = SPDE_FT$wave
  cosine_indices = SPDE_FT$cosine_indices
  n_cos = SPDE_FT$n_cos
  n_cos_sin = length(cosine_indices)
  K = SPDE_FT$K
  
  ## Compute obs0 indices ##
  obs0_indices = sapply(obs0, function(A) A[1]+A[2]*n+1)
  
  ## LGSSM matrices [indep of s, so outside loop]
  XBeta = XBeta_f(param_bias,n,t,n_obs)
  R = R_f(param_obs,n,n_obs)
  G = G_f(param_sig,wave,cosine_indices,dt=dt,n_cos)
  Q = Q_f(param_sig,wave,n,n_cos,nu,dt=dt,norm=TRUE,
          spec_weights_mat = spec_weights_mat)
  
  ## LGSSM (parameter) gradient arrays [indep of s, so outside loop] ##
  Phi_grad_param = Phi_grad_f(param_obs,param_sig,n,K,param_bias,n_obs)
  XBeta_grad_param = XBeta_grad_f(param_obs,param_sig,n,t,param_bias,n_obs)
  R_grad_param = R_grad_f(param_obs,param_sig,n,param_bias,n_obs)
  G_grad_param = G_grad_f(param_obs,param_sig,n,K,wave,cosine_indices,dt=dt,nu,n_cos,
                          param_bias)
  Q_grad_param = Q_grad_f(param_obs,param_sig,n,K,wave,n_cos,nu,dt=dt,norm=TRUE,
                          param_bias,spec_weights_mat = spec_weights_mat)
  
  ## LGSSM (space) gradient arrays [indep of s, so outside loop] ##
  XBeta_grad_space = XBeta_grad_space_f(n_obs,t)
  R_grad_space = R_grad_space_f(n_obs)
  G_grad_space = G_grad_space_f(n_obs,K)
  Q_grad_space = Q_grad_space_f(n_obs,K)
  
  ## Compute TKF 'grad_params' from 'grad_obs' ##  
  grad_params = as.numeric(sapply(grad_obs,function(x) c(2*x-1,2*x)))
  
  ## Initial values for MTKF ##
  if(is.null(mu0)){
    mu0 = rep(0,K)
  }
  if(is.null(Sigma0)){
    Sigma0 = Q_f(param_sig,wave,n,n_cos,nu=1,dt=dt,norm=TRUE,
                 spec_weights_mat = spec_weights_mat)
  }
  if(is.null(mu0_grad_param)){
    mu0_grad_param = rep(list(rep(0,K)),n_param)
  }
  if(is.null(Sigma0_grad_param)){
    Sigma0_grad_param = Q_grad_f(param_obs,param_sig,n,K,wave,n_cos,nu,dt=dt,norm=TRUE,
                                 param_bias,spec_weights_mat = spec_weights_mat)
  }
  if(is.null(mu0_grad_space)){
    mu0_grad_space = rep(list(rep(0,K)),2*n_obs)
  }
  if(is.null(Sigma0_grad_space)){
    Sigma0_grad_space = Q_grad_space_f(n_obs,K)
  }
  
  ## Run Batch Optimal Sensor Placement Algorithm ##
  for (i in 1:n_iterations){
    
    ## Current sensor locations ##
    obs_i = lapply(obs,function(A) A[i,])
    
    ## LGSSM matrices [dependent on s] ##
    Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs,coords=obs_i)
    
    ## LGSSM (space) gradient arrays [dependent on s] ##
    Phi_grad_space = Phi_grad_space_f(wave,cosine_indices,n_cos,n,n_obs,
                                            coords = obs_i,
                                            grad_coord_index = 1:n_obs)
    
    ## Run MTKF (parameters)
    MTKF_param = MultiTangentKalmanFilter(y = y, Phi = Phi, 
                                          R = R, G = G, Q = Q, 
                                          Phi_grad = Phi_grad_param, 
                                          R_grad = R_grad_param, 
                                          G_grad = G_grad_param, 
                                          Q_grad = Q_grad_param, 
                                          n = n_obs, t = t, K = K, log_lik = TRUE, 
                                          mu0 = mu0, Sigma0 = Sigma0, 
                                          mu0_grad = mu0_grad_param,
                                          Sigma0_grad = Sigma0_grad_param, 
                                          n_param = n_param,
                                          grad_params = 1:n_param,
                                          lin_pred = XBeta,
                                          lin_pred_grad = XBeta_grad_param)
    
    ## Log-likelihood ##
    ll[i] = unlist(MTKF_param$ll)
    
    ## Log-likelihood (parameter) gradient ##
    ll_grads_param[i,] = unlist(MTKF_param$ll_grad)
    
    # Objective function (parameter) gradient ##
    obj_func_grads_param[i,] = 
      objective_function_grad_param(n_param,t,MTKF_param$Sigma_i_i_grad,
                                    W_fourier)
    
    
    ## Run MTKF (space) ##
    MTKF_space = MultiTangentKalmanFilter(y = y, Phi = Phi, R = R, 
                                          G = G, Q = Q, Phi_grad = Phi_grad_space, 
                                          R_grad = R_grad_space, G_grad = G_grad_space, 
                                          Q_grad = Q_grad_space, n = n_obs, t = t, 
                                          K = K, log_lik = TRUE, mu0 = mu0, 
                                          Sigma0 = Sigma0, mu0_grad = mu0_grad_space,
                                          Sigma0_grad = Sigma0_grad_space,
                                          n_param = 2*n_obs, grad_params = grad_params, 
                                          lin_pred = XBeta, 
                                          lin_pred_grad = XBeta_grad_space)
    
    
    ## Log likelihood (space) gradient ##
    for (k in 1:n_obs){
      ll_grads_space[[k]][i,] = unlist(MTKF_space$ll_grad[(2*k-1):(2*k)])
    }
    
    ## Objective function ##
    obj_func[i] = objective_function(t,MTKF_space$Sigma_i_i,W_fourier)
    
    ## Objective function (space) gradient ##
    obj_func_grad = objective_function_grad_space(n_obs,t,MTKF_space$Sigma_i_i_grad,W_fourier)
    for (j in grad_obs){obj_func_grads_space[[j]][i,] = obj_func_grad[[j]]}
      
    ## Update observation locations ##
    for (j in grad_obs){
      obs[[j]][i+1,] = (obs[[j]][i,] - step_size[[j]][i,]*obj_func_grads_space[[j]][i,]) %% n
    }
    for (j in (1:n_obs)[-grad_obs]){
      obs[[j]][i+1,] = obs[[j]][i,]
    }
    
    
    
    ## Back Tracking Line Search ##
    if(back_track){
      
      ## Back Tracking Line Search For Each Sensor ##
      if(back_track_independent){
        
        checked = rep(FALSE,n_param)
        
        ## Loop over sensors ##
        for(l in back_track_obs){
          
          while(!checked[l]){
            
            ## Sensor locations with only l^th sensor updated ##
            obs_check = lapply(obs,function(A) A[i,])
            obs_check[[l]] = obs[[l]][i+1,]
            
            ## KF input(s) for this set of locations ##
            Phi_new = Phi_f(wave,cosine_indices,n_cos,n,n_obs,
                                coords=obs_check)
            
            ## Run KF for this set of locations ##
            KF = KalmanFilter(y=y,Phi=Phi_new,R=R,G=G,Q=Q,m0=mu0,
                              R0=Sigma0,log_lik=TRUE)
              
            ## Compute objective for this set of locations ##
            obj_func_update = objective_function(t,KF$R_i_i,W_fourier)
            
            ## If new objective < old objective, accept new locations ##
            
            if(obj_func_update < obj_func[i]){
              checked[l] = TRUE
            }
            
            checked[l] = TRUE
            
            ## If old objective <~ new objective ##
            if((obj_func_update-obj_func[i])>back_track_tol){
              
              print("a")
              
              ## Re-scale future step-sizes ##
              if(!back_track_time_independent){
                step_size[[l]][(i+1):(n_iterations),] = 
                  back_track_scale*step_size[[l]][(i+1):(n_iterations),]
              }
              
              ## Re-update sensor locations ##
              obs[[l]][i+1,] = 
                back_track_scale*((1/back_track_scale-1)*obs[[l]][i,] + obs[[l]][i+1,])
              
              ## Check again ##
              checked[l] = FALSE
            }
          }
        }
      }
      
      ## Back Tracking Line Search For All Parameters ##
      if(!back_track_independent){
        
        checked = FALSE
        
        while(!checked){
          
          ## Parameter vector with all elements updated ##
          obs_check = lapply(obs,function(A) A[i+1,])
          
          ## KF input(s) for this set of locations ##
          Phi_new = Phi_f(wave,cosine_indices,n_cos,n,n_obs,
                              coords=obs_check)
          
          ## Run KF for this set of locations ##
          KF = KalmanFilter(y=y,Phi=Phi_new,R=R,G=G,Q=Q,m0=mu0,
                            R0=Sigma0,log_lik=TRUE)
          
          ## Compute objective for this set of locations ##
          obj_func_update = objective_function(t,KF$R_i_i,W_fourier)
          
          ## If new objective < old objective, accept new locations ##
          if(obj_func_update < obj_func[i]){
            checked = TRUE
          }
          
          checked = TRUE
          
          ## If old objective <~ new objective ##
          if((obj_func_update-obj_func[i])>back_track_tol){
            
            ## Re-scale future step-sizes ##
            if(!back_track_time_independent){
              for (l in grad_obs){
                step_size[[l]][(i+1):(n_iterations),] = 
                  back_track_scale*step_size[[l]][(i+1):(n_iterations),]
              }
            }
            
            ## Re-update sensor locations ##
            for (l in grad_obs){
              obs[[l]][i+1,] = 
                back_track_scale*((1/back_track_scale-1)*obs[[l]][i,] + obs[[l]][i+1,])
            }
            
            ## Check again ##
            checked = FALSE
          }
        }
      }
    }
    
    
    
    ## Update initial values for TKF ##
    mu0 = MTKF_space$mu_i_i[t+1,]
    Sigma0 = MTKF_space$Sigma_i_i[t+1,,]
    mu0_grad_param = lapply(MTKF_param$mu_i_i_grad, function(A) A[t+1,])
    Sigma0_grad_param = lapply(MTKF_param$Sigma_i_i_grad, function(A) A[t+1,,])
    mu0_grad_space = lapply(MTKF_space$mu_i_i_grad, function(A) A[t+1,])
    Sigma0_grad_space = lapply(MTKF_space$Sigma_i_i_grad, function(A) A[t+1,,])
    
    ## Print iteration no. ##
    if(!is.null(print_iter_freq)){
      if(iter_numbers[i]%%print_iter_freq==0){
        cat("Iteration ", iter_numbers[i], "\n")
      }
    }
    
    ## Plot ##
    if(plot){
      
      if (i%%plot_freq==0){
        
        if(save_plot){
          if(i==n_iterations){
            setwd(fig_wd)
            pdf(filename_plot,width=10,height=5)
          }
        }
        
        ## Plot parameters ##
        par(mar=c(2.8,3.5,1.5,1.5),mgp=c(1.5,0.5,0))
        if(n_obs==1){par(mfrow=c(1,1))}
        if(n_obs==2){par(mfrow=c(1,2))}
        if(n_obs==3){par(mfrow=c(1,3))}
        if(n_obs==4){par(mfrow=c(2,2))}
        
        ## Plot Points ##
        plot_points = seq(1,i,plot_point_freq)
        
        ## Plot ##
        for(k in 1:n_obs){
          plot(plot_points,obs[[k]][plot_points,1],cex=.4,ylim=c(0,n),
               xlim=c(0,n_iterations),xlab="Iterations",
               ylab=expression(hat(s)),col=brewer.pal(8,"Greys")[6],
               main=paste("Sensor",k))
          points(plot_points,obs[[k]][plot_points,2],cex=.4,col=brewer.pal(8,"Greys")[8])
          legend(x=0.6*n_iterations,y=2.8,
                 legend=c(expression(hat(s[x])),expression(hat(s[y]))),
                 col=brewer.pal(8,"Greys")[c(6,8)],pch=19,ncol=2)
          
          ## Add dashed lines for weighted coordinates ##
          if(!is.null(W_coords)){
            abline(h=W_coords[[k]][1],lty=2,col=brewer.pal(8,"Reds")[3])
            abline(h=W_coords[[k]][2],lty=2,col=brewer.pal(8,"Reds")[7])
          }
        }
        
        if(save_plot){
          if(i==n_iterations){
            dev.off()
          }
        }
      }
    }
    
    ## 2D Plot ##
    if(plot2d){
      
      if (i%%plot_freq==0){
        
        if(save_plot2d){
          if(i==n_iterations){
            setwd(fig_wd)
            pdf(filename_plot2d,width=6,height=6)
          }
        }
        
        ## Plot parameters ##
        par(mar=c(4.0,3.5,2.5,1.5),mfrow=c(1,1),mgp=c(2,0.5,0))
        
        ## Plot Points ##
        plot_points = seq(1,i,plot_point_freq)
        
        ## Plot initial point ##
        if(!plot_initial){
          plot(obs[[1]][1,1],obs[[1]][1,2],xlim = c(0,n),ylim = c(0,n),
               xlab=expression(x),ylab=expression(y), 
               col=colorRampPalette(c("grey10", "grey100"))(100)[1])
        }
        
        if(plot_initial){
          plot(obs[[1]][1,1],obs[[1]][1,2],xlim = c(0,n),ylim = c(0,n),
               cex = 2, pch = 21, col=brewer.pal(8,"Greens")[7],
               bg = brewer.pal(8,"Greens")[2], xlab=expression(x),
               ylab=expression(y))
          
          for (k in 1:n_obs){
            points(obs[[k]][1,1],obs[[k]][1,2], cex = 2, pch = 21, 
                 col=brewer.pal(8,"Greens")[7],
                 bg = brewer.pal(8,"Greens")[2])
          }
        }
        
        ## Add legend ##
        if(leg){
          if(!is.null(W_coords)){
            legend(x="bottomright", inset = c(0.035,0.04), pt.cex = c(2,2,2), 
                   pch = c(21,21,21),
                   legend=c("Initial Sensor Locations",
                            "Weighted Observation Locations", 
                            "Final Sensor Locations"),
                   pt.bg = c(brewer.pal(8,"Greens")[2],
                             brewer.pal(8,"Reds")[2],
                             brewer.pal(8,"Blues")[2]),
                   col = c(brewer.pal(8,"Greens")[7],
                           brewer.pal(8,"Reds")[7],
                           brewer.pal(8,"Blues")[7]))
          }
          if(is.null(W_coords)){
            legend(x="bottomright", inset = c(0.035,0.04), pt.cex = c(2,2,2), 
                   pch = c(21,21,21),
                   legend=c("Initial Sensor Locations",
                            "Final Sensor Locations"),
                   pt.bg = c(brewer.pal(8,"Greens")[2],
                             brewer.pal(8,"Blues")[2]),
                   col = c(brewer.pal(8,"Greens")[7],
                           brewer.pal(8,"Blues")[7]))
          }
        }
        
        ## Add grid ##
        if(grid){
          abline(h=0:n,v=0:n,col="gray90",lty=2)
        }
        
        ## Highlight weighted coordinates ##
        if(!is.null(W_coords)){
          for(k in 1:n_obs){
            points(W_coords[[k]][1],W_coords[[k]][2],cex=2,
                   pch=21,col=brewer.pal(8,"Reds")[7],
                   bg = brewer.pal(8,"Reds")[2])
          }
        }
        
        ## Plot final point ##
        if(i==n_iterations){
          for (k in 1:n_obs){
            points(obs[[k]][n_iterations,1],obs[[k]][n_iterations,2],
                   cex = 2, pch = 21, col = brewer.pal(8,"Blues")[7],
                   bg = brewer.pal(8,"Blues")[2])
          }
        }
        
        ## Plot remaining points ##
        for (k in 1:n_obs){
          for (l in plot_points){
            points(obs[[k]][l,1],obs[[k]][l,2],
                   col=colorRampPalette(c("grey10", "grey100"))(100)[5*k-4],cex=.6)
          }
        }
        
        if(save_plot2d){
          if(i==n_iterations){
            dev.off()
          }
        }
      }
    }
  }
  
  ## Output ##
  output = list(obs = obs, obj_func = obj_func, 
                obj_func_grads_param = obj_func_grads_param,
                obj_func_grads_space = obj_func_grads_space,
                ll = ll, ll_grads_param = ll_grads_param,
                ll_grads_space = ll_grads_space, 
                mu0 = mu0, Sigma0 = Sigma0, 
                mu0_grad_param = mu0_grad_param,
                Sigma0_grad_param = Sigma0_grad_param,
                mu0_grad_space = mu0_grad_space, 
                Sigma0_grad_space = Sigma0_grad_space)
  
  if(back_track){
    output = c(output,list(step_size = step_size))
  }
  
  ## Save ##
  if(save){
    saveRDS(output, file = filename)
  }
  
  ## Output ##
  return(output)
}

####################################

###################################
# ONLINE OPTIMAL SENSOR PLACEMENT #
###################################

RecursOptSens <- function(y, n_obs = dim(y)[2], t = dim(y)[1], K, 
                          param_obs, param_sig, param_bias = param_bias_func(),
                          obs0, n_iterations, step_size, W_fourier,
                          spec_weights_mat = NULL, plot = FALSE, 
                          W_coords = NULL, grad_obs = 1:n_obs, plot2d = FALSE,
                          plot_point_freq = 1, plot_freq = 500,
                          save_plot = FALSE, 
                          filename_plot = "RecursOptSens_plot.pdf",
                          save_plot2d = FALSE,
                          filename_plot2d = "RecursOptSens_2dplot.pdf",
                          grid = TRUE, print_iter_freq = 100,
                          save = FALSE, filename = "RecursOptSens",
                          iter_numbers = 1:(n_iterations*t),
                          mu0 = NULL, Sigma0 = NULL, mu0_grad_param = NULL,
                          Sigma0_grad_param = NULL, mu0_grad_space = NULL,
                          Sigma0_grad_space = NULL, back_track = FALSE,
                          back_track_tol = 1e-3, back_track_independent = TRUE,
                          back_track_time_independent = FALSE, 
                          back_track_scale = 0.5,
                          back_track_obs = 1:n_obs,
                          plot_initial = FALSE,
                          leg = FALSE,
                          spde_sim = NULL, dt = 1){
  
  ## Number of parameters ##
  n_param = length(param_sig) + param_obs$m + param_bias$m
  
  ## Matrices of observation locations ##
  obs = rep(list(matrix(0,nrow=n_iterations*t+1,ncol=2)),n_obs)
  for (i in 1:n_obs){
    obs[[i]][1,] = obs0[[i]]
  }
  
  ## Matrices of log-lik, conditional log-lik gradients, 2nd derivatives ##
  ll = rep(0,t*n_iterations)
  ll_grads_param = matrix(0,ncol=n_param,nrow=t*n_iterations)
  ll_grads_space = rep(list(matrix(0,ncol=2,nrow=t*n_iterations)),n_obs)
  ll_grads2_param = matrix(0,ncol=n_param,nrow=t*n_iterations)
  
  ## Matrices of obj. func., conditional obj. func. grads ##
  obj_func = matrix(0,ncol=1,nrow=t*n_iterations)
  obj_func_grads_param = matrix(0,ncol=n_param,nrow=t*n_iterations)
  obj_func_grads_space = rep(list(matrix(0,ncol=2,nrow=t*n_iterations)),n_obs)
  
  ## Matrix of KF means ##
  kf_means = matrix(0,ncol=K,nrow=t*n_iterations)
  
  ## SPDE_FT object ##
  SPDE_FT = spde_initialise(n,t,K)
  wave = SPDE_FT$wave
  cosine_indices = SPDE_FT$cosine_indices
  n_cos = SPDE_FT$n_cos
  n_cos_sin = length(cosine_indices)
  K = SPDE_FT$K
  
  ## Compute obs0 indices ##
  obs0_indices = sapply(obs0, function(A) A[1]+A[2]*n+1)
  
  ## LGSSM matrices [indep of s, so outside loop]
  XBeta = XBeta_f(param_bias,n,t,n_obs)
  R = R_f(param_obs,n,n_obs)
  G = G_f(param_sig,wave,cosine_indices,dt=dt,n_cos)
  Q = Q_f(param_sig,wave,n,n_cos,nu,dt=dt,norm=TRUE,
          spec_weights_mat = spec_weights_mat)
  
  ## LGSSM (parameter) gradient arrays [indep of s, so outside loop] ##
  Phi_grad_param = Phi_grad_f(param_obs,param_sig,n,K,param_bias,n_obs)
  XBeta_grad_param = XBeta_grad_f(param_obs,param_sig,n,t,param_bias,n_obs)
  R_grad_param = R_grad_f(param_obs,param_sig,n,param_bias,n_obs)
  G_grad_param = G_grad_f(param_obs,param_sig,n,K,wave,cosine_indices,dt=dt,nu,n_cos,
                          param_bias)
  Q_grad_param = Q_grad_f(param_obs,param_sig,n,K,wave,n_cos,nu,dt=dt,norm=TRUE,
                          param_bias,spec_weights_mat = spec_weights_mat)
  
  ## LGSSM (space) gradient arrays [indep of s, so outside loop]
  XBeta_grad_space = XBeta_grad_space_f(n_obs,t)
  R_grad_space = R_grad_space_f(n_obs)
  G_grad_space = G_grad_space_f(n_obs,K)
  Q_grad_space = Q_grad_space_f(n_obs,K)
  
  ## Compute TKF 'grad_params' from 'grad_obs' ##  
  grad_params = as.numeric(sapply(grad_obs,function(x) c(2*x-1,2*x)))
  
  
  ## Run Recursive Sensor Placement Algorithm ##
  
  ## Initial values for MTKF ##
  if(is.null(mu0)){
    mu0 = rep(0,K)
  }
  if(is.null(Sigma0)){
    Sigma0 = Q_f(param_sig,wave,n,n_cos,nu=1,dt=dt,norm=TRUE,
                 spec_weights_mat = spec_weights_mat)
  }
  if(is.null(mu0_grad_param)){
    mu0_grad_param = rep(list(rep(0,K)),n_param)
  }
  if(is.null(Sigma0_grad_param)){
    Sigma0_grad_param = Q_grad_f(param_obs,param_sig,n,K,wave,n_cos,nu,dt=dt,norm=TRUE,
                                 param_bias,spec_weights_mat = spec_weights_mat)
  }
  if(is.null(mu0_grad_space)){
    mu0_grad_space = rep(list(rep(0,K)),2*n_obs)
  }
  if(is.null(Sigma0_grad_space)){
    Sigma0_grad_space = Q_grad_space_f(n_obs,K)
  }
  
  for (j in 1:n_iterations){
    for (i in 1:t){
      
      ## Current sensor locations ##
      obs_i = lapply(obs,function(A) A[(j-1)*t+i,])
      
      ## Inputs for KF ##
      Phi = Phi_f(wave,cosine_indices,n_cos,n,n_obs,coords=obs_i)
      
      ## Inputs for TKF ##
      Phi_grad_space = Phi_grad_space_f(wave,cosine_indices,n_cos,n,
                                              n_obs,coords = obs_i,
                                              grad_coord_index = 1:n_obs)
      
      ## Generate y-simulation ##
      
      ## The first is technically correct (i.e. generate new observations for
      ## the new sensor locations); but since the Ricatti equation is indep.
      ## of y, we can just use pre-simulated y as 'dummy' input ##
      
      #y_i = y_simulation(spde_sim=spde_sim,t_vals=i,n_obs=n_obs,
                           #obs_loc=obs_i)
      y_i = matrix(y[i,],nrow=1)
        
      
      ## Run MTKF (parameters) ##
      MTKF_param = MultiTangentKalmanFilter(y = y_i, 
                                            Phi = Phi, R = R, G = G, Q = Q, 
                                            Phi_grad = Phi_grad_param, 
                                            R_grad = R_grad_param, 
                                            G_grad = G_grad_param, 
                                            Q_grad = Q_grad_param, 
                                            n = n_obs, t = 1, K = K, log_lik = TRUE, 
                                            mu0 = mu0, Sigma0 = Sigma0, 
                                            mu0_grad = mu0_grad_param,
                                            Sigma0_grad = Sigma0_grad_param, 
                                            n_param = n_param,
                                            grad_params = 1:n_param,
                                            lin_pred = XBeta,
                                            lin_pred_grad = XBeta_grad_param)
      
      ## Log-likelihood ##
      ll[(j-1)*t+i] = unlist(MTKF_param$ll)
      
      ## Log-likelihood (parameter) gradient ##
      ll_grads_param[(j-1)*t+i,] = unlist(MTKF_param$ll_grad)
      
      # Objective function (parameter) gradient ##
      obj_func_grads_param[(j-1)*t+i,] = 
        objective_function_grad_param(n_param,t=1,MTKF_param$Sigma_i_i_grad,
                                      W_fourier)
      
      ## Run MTKF (space) ##
      MTKF_space = MultiTangentKalmanFilter(y = y_i, Phi = Phi, 
                                            R = R, G = G, Q = Q, Phi_grad = Phi_grad_space, 
                                            R_grad = R_grad_space, G_grad = G_grad_space, 
                                            Q_grad = Q_grad_space, n = n_obs, t = 1, 
                                            K = K, log_lik = TRUE, mu0 = mu0, 
                                            Sigma0 = Sigma0, mu0_grad = mu0_grad_space,
                                            Sigma0_grad = Sigma0_grad_space,
                                            n_param = 2*n_obs, 
                                            grad_params = grad_params, 
                                            lin_pred = XBeta, 
                                            lin_pred_grad = XBeta_grad_space)
      
      
      ## Log likelihood (space) gradient ##
      for (k in 1:n_obs){
        ll_grads_space[[k]][(j-1)*t+i,] = unlist(MTKF_space$ll_grad[(2*k-1):(2*k)])
      }
      
      ## Objective function ##
      obj_func[(j-1)*t+i] = objective_function(t=1,MTKF_space$Sigma_i_i,W_fourier)
      
      ## Objective function (space) gradient ##
      obj_func_grad = objective_function_grad_space(n_obs,t=1,MTKF_space$Sigma_i_i_grad,W_fourier)
      for (k in grad_obs){
        obj_func_grads_space[[k]][(j-1)*t+i,] = obj_func_grad[[k]]
      }
      
      ## KF means ##
      kf_means[(j-1)*t+i,] = MTKF_space$mu_i_i[2,]
      
      
      ## Update observation location ##
      for (k in grad_obs){
        obs[[k]][(j-1)*t+i+1,] = 
          (obs[[k]][(j-1)*t+i,] - step_size[[k]][(j-1)*t+i,]*obj_func_grads_space[[k]][(j-1)*t+i,]) %% n
      }
      for(k in (1:n_obs)[-grad_obs]){
        obs[[k]][(j-1)*t+i+1,] = obs[[k]][(j-1)*t+i,]
      }
      
      
      
      ## Back Tracking Line Search ##
      if(back_track){
        
        ## Back Tracking Line Search For Each Sensor ##
        if(back_track_independent){
          
          checked = rep(FALSE,n_param)
          
          ## Loop over sensors ##
          for(l in back_track_obs){
            
            while(!checked[l]){
              
              ## Sensor locations with only l^th sensor updated ##
              obs_check = lapply(obs,function(A) A[(j-1)*t+i,])
              obs_check[[l]] = obs[[l]][(j-1)*t+i+1,]
              
              ## KF input(s) for this set of locations ##
              Phi_new = Phi_f(wave,cosine_indices,n_cos,n,n_obs,
                                  coords=obs_check)
              
              ## Run KF for this set of locations ##
              KF = KalmanFilter(y=matrix(y[i,],nrow=1),t=1,
                                Phi=Phi_new,R=R,G=G,Q=Q,m0=mu0,
                                R0=Sigma0,log_lik=TRUE)
              
              ## Compute objective for this set of locations ##
              obj_func_update = objective_function(t=1,KF$R_i_i,W_fourier)
              
              ## If new objective < old objective, accept new locations ##
              if(obj_func_update < obj_func[(j-1)*t+i]){
                checked[l] = TRUE
              }
              
              checked[l] = TRUE
              
              ## If old objective <~ new objective ##
              if((obj_func_update-obj_func[(j-1)*t+i])>back_track_tol){
                
                ## Re-scale future step-sizes ##
                if(!back_track_time_independent){
                  step_size[[l]][((j-1)*t+i+1):(t*n_iterations),] = 
                    back_track_scale*step_size[[l]][((j-1)*t+i+1):(t*n_iterations),]
                }
                
                ## Re-update sensor locations ##
                obs[[l]][(j-1)*t+i+1,] = 
                  back_track_scale*((1/back_track_scale-1)*obs[[l]][(j-1)*t+i,] + obs[[l]][(j-1)*t+i+1,])
                
                ## Check again ##
                checked[l] = FALSE
              }
            }
          }
        }
        
        ## Back Tracking Line Search For All Parameters ##
        if(!back_track_independent){
          
          checked = FALSE
          
          while(!checked){
            
            ## Parameter vector with all elements updated ##
            obs_check = lapply(obs,function(A) A[(j-1)*t+i+1,])
            
            ## KF input(s) for this set of locations ##
            Phi_new = Phi_f(wave,cosine_indices,n_cos,n,n_obs,
                                coords=obs_check)
            
            ## Run KF for this set of locations ##
            KF = KalmanFilter(y=matrix(y[i,],nrow=1),t=1,Phi=Phi_new,
                              R=R,G=G,Q=Q,m0=mu0,
                              R0=Sigma0,log_lik=TRUE)
            
            ## Compute objective for this set of locations ##
            obj_func_update = objective_function(t=1,KF$R_i_i,W_fourier)
            
            ## If new objective < old objective, accept new locations ##
            if(obj_func_update < obj_func[(j-1)*t+i]){
              checked = TRUE
            }
            
            checked = TRUE
            
            ## If old objective <~ new objective ##
            if((obj_func_update-obj_func[(j-1)*t+i])>back_track_tol){
              
              ## Re-scale future step-sizes ##
              if(!back_track_time_independent){
                for (l in grad_obs){
                  step_size[[l]][((j-1)*t+i+1):(t*n_iterations),] = 
                    back_track_scale*step_size[[l]][((j-1)*t+i+1):(t*n_iterations),]
                }
              }
              
              ## Re-update sensor locations ##
              for (l in grad_obs){
                obs[[l]][(j-1)*t+i+1,] = 
                  back_track_scale*((1/back_track_scale-1)*obs[[l]][(j-1)*t+i,] + obs[[l]][(j-1)*t+i+1,])
              }
              
              ## Check again ##
              checked = FALSE
            }
          }
        }
      }
      
      ## Update initial values for TKF ##
      mu0 = MTKF_space$mu_i_i[2,]
      Sigma0 = MTKF_space$Sigma_i_i[2,,]
      mu0_grad_param = lapply(MTKF_param$mu_i_i_grad, function(A) A[2,])
      Sigma0_grad_param = lapply(MTKF_param$Sigma_i_i_grad, function(A) A[2,,])
      mu0_grad_space = lapply(MTKF_space$mu_i_i_grad, function(A) A[2,])
      Sigma0_grad_space = lapply(MTKF_space$Sigma_i_i_grad, function(A) A[2,,])
      
      ## Print iteration no. ##
      if(!is.null(print_iter_freq)){
        if(iter_numbers[(j-1)*t+i]%%print_iter_freq==0){
          cat("Iteration ", iter_numbers[(j-1)*t+i], "\n")
        }
      }
      
      ## Plot ##
      if(plot){
        
        if (((j-1)*t+i)%%plot_freq==0){
          
          if(save_plot){
            if(((j-1)*t+i)==(n_iterations*t)){
              setwd(fig_wd)
              pdf(filename_plot,width=9,height=4.5)
            }
          }
          
          ## Plot parameters ##
          par(mar=c(2.8,4.5,1.5,1))
          if(n_obs==1){par(mfrow=c(1,1))}
          if(n_obs==2){par(mfrow=c(1,2))}
          if(n_obs==3){par(mfrow=c(1,3))}
          if(n_obs==4){par(mfrow=c(2,2))}
          if(n_obs==8){par(mfrow=c(2,4))}
          
          ## Plot Points ##
          plot_points = seq(1,((j-1)*t+i),plot_point_freq)
          
          ## Plot ##
          for(k in 1:n_obs){
            plot(plot_points,obs[[k]][plot_points,1],cex=.4,ylim=c(0,n),
                 xlim=c(0,n_iterations*t),xlab="Time",
                 ylab=expression(hat(s)),col=brewer.pal(8,"Greys")[6],
                 main=paste("Sensor",k))
            points(plot_points,obs[[k]][plot_points,2],cex=.4,
                   col=brewer.pal(8,"Greys")[8])
            legend(x=0.5*(n_iterations*t),y=3.2,
                   legend=c(expression(hat(s[x])),expression(hat(s[y]))),
                   col=brewer.pal(8,"Greys")[c(6,8)],pch=19,ncol=2)
            
            ## Add dashed lines for weighted coordinates ##
            if(!is.null(W_coords)){
              abline(h=W_coords[[k]][1],lty=2,col=brewer.pal(8,"Reds")[3])
              abline(h=W_coords[[k]][2],lty=2,col=brewer.pal(8,"Reds")[7])
            }
          }
          
          if(save_plot){
            if(((j-1)*t+i)==(n_iterations*t)){
              dev.off()
            }
          }
        }
      }
      
      ## 2d Plot ##
      if(plot2d){
        
        if (((j-1)*t+i)%%plot_freq==0){
          
          if(save_plot2d){
            if(((j-1)*t+i)==(n_iterations*t)){
              setwd(fig_wd)
              pdf(filename_plot2d,width=6,height=6)
            }
          }
          
          ## Plot parameters ##
          par(mar=c(4.0,3.5,2.5,1.5),mfrow=c(1,1),mgp=c(2,0.5,0))
          
          ## Plot points ##
          plot_points = seq(1,((j-1)*t+i),plot_point_freq)
          
          ## Plot initial point ##
          if(!plot_initial){
            plot(obs[[1]][1,1],obs[[1]][1,2],xlim = c(0,n),ylim = c(0,n),
                 xlab=expression(x),ylab=expression(y), 
                 col=colorRampPalette(c("grey10", "grey100"))(100)[1])
          }
          
          if(plot_initial){
            plot(obs[[1]][1,1],obs[[1]][1,2],xlim = c(0,n),ylim = c(0,n),
                 cex = 2, pch = 21, col=brewer.pal(8,"Greens")[7],
                 bg = brewer.pal(8,"Greens")[2], xlab=expression(x),
                 ylab=expression(y))
            
            for (k in 1:n_obs){
              points(obs[[k]][1,1],obs[[k]][1,2], cex = 2, pch = 21, 
                     col=brewer.pal(8,"Greens")[7],
                     bg = brewer.pal(8,"Greens")[2])
            }
          }
          
          ## Add legend ##
          if(leg){
            if(!is.null(W_coords)){
              legend(x="bottomright", inset = c(0.035,0.04), pt.cex = c(2,2,2), 
                     pch = c(21,21,21),
                     legend=c("Initial Sensor Locations",
                              "Weighted Observation Locations", 
                              "Final Sensor Locations"),
                     pt.bg = c(brewer.pal(8,"Greens")[2],
                              brewer.pal(8,"Reds")[2],
                              brewer.pal(8,"Blues")[2]),
                     col = c(brewer.pal(8,"Greens")[7],
                             brewer.pal(8,"Reds")[7],
                             brewer.pal(8,"Blues")[7]))
            }
            if(is.null(W_coords)){
              legend(x="bottomright", inset = c(0.035,0.04), pt.cex = c(2,2,2), 
                     pch = c(21,21,21),
                     legend=c("Initial Sensor Locations",
                              "Final Sensor Locations"),
                     pt.bg = c(brewer.pal(8,"Greens")[2],
                               brewer.pal(8,"Blues")[2]),
                     col = c(brewer.pal(8,"Greens")[7],
                             brewer.pal(8,"Blues")[7]))
            }
          }
          
          ## Add grid ##
          if(grid){
            abline(h=0:n,v=0:n,col="gray90",lty=2)
          }
          
          ## Highlight weighted coordinates ##
          if(!is.null(W_coords)){
            for(k in 1:n_obs){
              points(W_coords[[k]][1],W_coords[[k]][2],cex=2,
                     pch=21,col=brewer.pal(8,"Reds")[7],
                     bg = brewer.pal(8,"Reds")[2])
            }
          }
          
          ## Plot final point ##
          if(((j-1)*t+i)==(n_iterations*t)){
            for (k in 1:n_obs){
            points(obs[[k]][n_iterations*t,1],obs[[k]][n_iterations*t,2],
                   cex = 2, pch = 21, col = brewer.pal(8,"Blues")[7],
                   bg = brewer.pal(8,"Blues")[2])
            }
          }
          
          ## Plot remaining points ##
          for (k in 1:n_obs){
            for (l in plot_points){
              points(obs[[k]][l,1],obs[[k]][l,2],
                     col=colorRampPalette(c("grey10", "grey100"))(100)[5*k-4],cex=.6)
            }
          }
          
        
          
          if(save_plot2d){
            if(((j-1)*t+i)==(n_iterations*t)){
              dev.off()
            }
          }
        }
      }
      
    }
  }
  
  ## Output ##
  output = list(obs = obs, obj_func = obj_func, 
                obj_func_grads_param = obj_func_grads_param,
                obj_func_grads_space = obj_func_grads_space,
                ll = ll, ll_grads_param = ll_grads_param,
                ll_grads_space = ll_grads_space,
                mu0 = mu0, Sigma0 = Sigma0, 
                mu0_grad_param = mu0_grad_param,
                Sigma0_grad_param = Sigma0_grad_param,
                mu0_grad_space = mu0_grad_space, 
                Sigma0_grad_space = Sigma0_grad_space,
                kf_means = kf_means)
  
  if(back_track){
    output = c(output,list(step_size = step_size))
  }
  
  ## Save ##
  if(save){
    saveRDS(output, file = filename)
  }
  
  ## Output ##
  return(output)
}

###################################

# ----

####################
# JOINT RML & ROSP #
####################

RML_ROSP <- function(spde_sim, n, n_obs, t, t_RML, t_ROSP, K, n_param, 
                     param_obs0,param_sig0, param_bias0, log_scale = FALSE, 
                     log_indices,grad_params, hessian_indices, obs0, grad_obs, 
                     n_iterations, step_size_RML, step_size_ROSP, W_fourier,
                     spec_weights_mat = NULL, RML_plot = FALSE, 
                     plot_limits = NULL, param_obs_true = NULL, 
                     param_sig_true = NULL, param_bias_true = NULL, 
                     param_obs_names = NULL, param_sig_names = NULL, 
                     param_bias_names = NULL, mean_param = FALSE, 
                     param_est_starts = NULL, param_est_ends = NULL, 
                     RML_plot_point_freq = 1, RML_plot_freq = 100, 
                     save_RML_plot = FALSE, RML_plot_filename = "RML_plot.pdf", 
                     ROSP_plot = FALSE, ROSP_plot2d = FALSE, W_coords = NULL, 
                     ROSP_plot_point_freq = 1, ROSP_plot_freq = 100,
                     save_ROSP_plot = FALSE, save_ROSP_plot2d = FALSE, 
                     ROSP_plot_filename = "ROSP_plot.pdf",
                     ROSP_plot2d_filename = "ROSP_plot2d.pdf",
                     ROSP_grid = TRUE, save_RML = FALSE, RML_filename = "RML",
                     save_ROSP = FALSE, ROSP_filename = "ROSP", 
                     RML_print_iter_freq = NULL, ROSP_print_iter_freq = NULL,
                     print_iter_freq = NULL, 
                     iter_numbers = 1:(n_iterations*t),
                     mu0 = NULL, Sigma0 = NULL, mu0_grad_param = NULL,
                     mu0_grad_space = NULL, Sigma0_grad_param = NULL, 
                     Sigma0_grad_space = NULL,mu0_grad2_param = NULL, 
                     Sigma0_grad2_param = NULL,
                     RML_back_track = FALSE,
                     RML_back_track_tol = 1e-3, 
                     RML_back_track_independent = TRUE,
                     RML_back_track_time_independent = FALSE, 
                     RML_back_track_scale = 0.5,
                     RML_back_track_param = 1:n_param,
                     ROSP_back_track = FALSE,
                     ROSP_back_track_tol = 1e-3, 
                     ROSP_back_track_independent = TRUE,
                     ROSP_back_track_time_independent = FALSE, 
                     ROSP_back_track_scale = 0.5,
                     ROSP_back_track_obs = 1:n_obs,
                     ROSP_plot_initial = FALSE,
                     ROSP_leg = FALSE,
                     y_sim_t_dep = TRUE,
                     dt = 1){
  
  ## ------------------------------------------------
  
  ## Number of parameters ##
  n_param_sig = length(param_sig0)
  n_param_obs = param_obs0$m
  n_param_bias = param_bias0$m
  
  ## Parameter indices ##
  param_sig_indices = 1:n_param_sig
  param_obs_indices = (n_param_sig+1):(n_param_sig+n_param_obs)
  if(n_param_bias>0){param_bias_indices = (n_param_sig+n_param_obs+1):n_param}
  
  ## ------------------------------------------------
  
  ## Matrix of parameter estimates ##
  param_est = matrix(0,ncol=n_param,nrow=t*n_iterations+1)
  
  ## Initial parameter estimate ##
  param0 = c(param_sig0,as.numeric(unlist(param_obs0$tau2_vals)))
  if(param_bias0$m>0){
    param0 = c(param0,as.numeric(unlist(param_bias0$bias_vals)))
  }
  
  ## Rescale if param0 on log scale
  if(log_scale){param0[log_indices] <- exp(param0[log_indices])}
  
  ## Initial parameter estimate ##
  param_est[1,] = param0
  
  ## ------------------------------------------------
  
  ## Vector of log-likelihoods ##
  log_liks = rep(0,t*n_iterations+1)
  
  ## Vector of objective functions ##
  obj_funcs = rep(0,t*n_iterations+1)
  
  ## Matrix of KF means ##
  kf_means = matrix(0,nrow=t*n_iterations+1,ncol=K)
  
  ## ------------------------------------------------
  
  ## Matrices of observation locations ##
  obs = rep(list(matrix(0,ncol=2,nrow=n_iterations*t+1)),n_obs)
  
  ## Initial observation locations ##
  for (i in 1:n_obs){obs[[i]][1,] = obs0[[i]]}
  
  ## ------------------------------------------------
  
  ## Observation simulation ##
  y_sim = matrix(0,nrow=t,ncol=n_obs)
  
  if(!y_sim_t_dep){
    y_sim = y_simulation(spde_sim = spde_sim, t_vals = 1:t, 
                         n_obs = n_obs, obs_loc = obs0, 
                         obs_loc_t_dep = FALSE)
  }
  
  ## ------------------------------------------------
  
  ## Number of time epochs ##
  M = t/(t_RML+t_ROSP)
  
  ## ------------------------------------------------
  
  ## SPDE_FT object ##
  SPDE_FT = spde_initialise(n,t,K)
  wave = SPDE_FT$wave
  cosine_indices = SPDE_FT$cosine_indices
  n_cos = SPDE_FT$n_cos
  n_cos_sin = length(cosine_indices)
  K = SPDE_FT$K
  
  ## ------------------------------------------------
  
  ## Initial values ##
  
  if(is.null(mu0)){
    mu0 = rep(0,K)
  }
  if(is.null(Sigma0)){
    Sigma0 = Q_f(param_sig0,wave,n,n_cos=n_cos,nu=1,dt=dt,norm=TRUE,
                 spec_weights_mat = spec_weights_mat)
  }
  if(is.null(mu0_grad_param)){
    mu0_grad_param = rep(list(rep(0,K)),n_param)
  }
  if(is.null(mu0_grad_space)){
    mu0_grad_space = rep(list(rep(0,K)),2*n_obs)
  }
  if(is.null(Sigma0_grad_param)){
    Sigma0_grad_param = Q_grad_f(param_obs0,param_sig0,n,K,wave,n_cos,nu,dt=dt,
                           norm=TRUE,param_bias0,
                           spec_weights_mat = spec_weights_mat)
  }
  if(is.null(Sigma0_grad_space)){
    Sigma0_grad_space = Q_grad_space_f(n_obs,K)
  }
  if(is.null(mu0_grad2_param)){
    mu0_grad2_param = rep(list(rep(0,K)),n_param)
  }
  if(is.null(Sigma0_grad2_param)){
    Sigma0_grad2_param = Q_grad2_f(param_obs0,param_sig0,n,K,wave,n_cos,nu,dt=dt,
                             norm=TRUE,param_bias0,
                             spec_weights_mat = spec_weights_mat)
  }
  
  ## ------------------------------------------------
  
  ## Joint RML and ROSP ##
  
  for (i in 1:M){
    
    ## ----------------------------------
     
    ## Times for current epoch ##
    t_start = (i-1)*(t_RML+t_ROSP)+1
    t_end = i*(t_RML+t_ROSP)
    t_vals = t_start:t_end
    
    ## RML times in current epoch ##
    t_start_RML = t_start
    t_end_RML = (i-1)*(t_RML+t_ROSP)+t_RML
    t_vals_RML = t_start_RML:t_end_RML
    
    ## ROSP times in current epoch ##
    t_start_ROSP = t_end_RML + 1
    t_end_ROSP = t_end
    t_vals_ROSP = t_start_ROSP:t_end_ROSP
    
    ## ----------------------------------
    
    ## Current parameters ##
    param_obs = param_obs_func(param_obs0$m,param_obs0$m_indices,
                               param_est[t_start,param_obs_indices])
    param_sig = param_est[t_start,param_sig_indices]
    
    if(n_param_bias>0){
      param_bias = param_bias_func(param_bias0$m,param_bias0$m_indices,
                                   param_est[t_start,param_bias_indices])
    }
    else{
      param_bias = param_bias_func()
    }
    
    ## ----------------------------------
    
    ## Current observation locations ##
    obs_current = lapply(1:n_obs,function(j) obs[[j]][t_start,])

    if(y_sim_t_dep){
      coords_RML = obs_current
    } 
    if(!y_sim_t_dep){
      coords_RML = obs0
    }
    
    ## ----------------------------------
    
    ## Simulate observations for t = t_vals ##
    if(y_sim_t_dep){
      y_sim_current = y_simulation(spde_sim = spde_sim, t_vals = t_vals, 
                                   n_obs = n_obs, obs_loc = obs_current, 
                                   obs_loc_t_dep = FALSE)
      
      ## Update 'y_sim' ##
      y_sim[t_vals,] = y_sim_current
    }
    
    ## RML 'y_sim' ##
    y_sim_RML = matrix(y_sim[t_vals_RML,],nrow=length(t_vals_RML))
    
    ## ROSP 'y_sim' & 'spde_sim' ##
    y_sim_ROSP = matrix(y_sim[t_vals_ROSP,],nrow=length(t_vals_ROSP))
    
    ## NB: technically speaking, we should generate new observations as the 
    ## inputs for ROSP, as the sensors move. However, since Ricatti equation is 
    ## independent of y, it is cheaper (and more convenient) to use this pre-
    ## simulated 'dummy' y input. Note that the sensor placement does now
    ## depend on y explicitly, as the parameter estimates (hence the model)
    ## depend on the y, and the sensor placement depends on the the model.
    ## ----------------------------------
    
    ## Run RML for t=t_vals_RML ## 
    RML_sim = RML(y = y_sim_RML, n = n, t = t_RML, K = K, 
                  param_obs0 = param_obs, param_sig0 = param_sig, 
                  param_bias0 = param_bias, n_param = n_param, 
                  grad_params = grad_params, 
                  step_size = matrix(step_size_RML[t_vals_RML,],nrow=length(t_vals_RML)), 
                  log_scale = log_scale, log_indices = log_indices, 
                  hessian_indices = hessian_indices, n_iterations = 1, 
                  n_obs = n_obs, coords = coords_RML, plot = FALSE, 
                  plot_limits = NULL, param_obs_true = NULL, 
                  param_sig_true = NULL, param_bias_true = NULL,
                  param_obs_names = NULL, param_sig_names = NULL,
                  param_bias_names = NULL, param_est_starts = NULL, 
                  param_est_ends = NULL, plot_point_freq = NULL, 
                  plot_freq = NULL, save_plot = FALSE, plot_filename = NULL, 
                  save = FALSE, filename = NULL, mean_param = FALSE, 
                  print_iter_freq = RML_print_iter_freq,
                  iter_numbers = iter_numbers[t_vals_RML],
                  mu0 = mu0, Sigma0 = Sigma0, mu0_grad_param = mu0_grad_param,
                  Sigma0_grad_param = Sigma0_grad_param, 
                  mu0_grad_space = mu0_grad_space, 
                  Sigma0_grad_space = Sigma0_grad_space,
                  mu0_grad2_param = mu0_grad2_param,
                  Sigma0_grad2_param = Sigma0_grad2_param,
                  back_track = RML_back_track,
                  back_track_tol = RML_back_track_tol, 
                  back_track_independent = RML_back_track_independent,
                  back_track_time_independent = RML_back_track_time_independent, 
                  back_track_scale = RML_back_track_scale,
                  back_track_param = RML_back_track_param,
                  W_fourier = W_fourier,spec_weights_mat = spec_weights_mat,
                  dt = dt)
    
    ## ----------------------------------
    
    ## Update step sizes ##
    if(RML_back_track){
      step_size_RML = matrix(rep(RML_sim$step_size[t_RML,],t),nrow=t,byrow=T)
    }
    
    ## ----------------------------------
    
    ## Update 'initial values' for KF ##
    mu0 = RML_sim$mu0
    Sigma0 = RML_sim$Sigma0
    mu0_grad_param = RML_sim$mu0_grad_param
    Sigma0_grad_param = RML_sim$Sigma0_grad_param
    mu0_grad_space = RML_sim$mu0_grad_space
    Sigma0_grad_space = RML_sim$Sigma0_grad_space
    if(!is.null(hessian_indices)){
      mu0_grad2_param = RML_sim$mu0_grad2_param
      Sigma0_grad2_param = RML_sim$Sigma0_grad2_param
    }
    
    ## ----------------------------------
    
    ## Update 'param_est' ##
    param_est[t_vals_RML+1,] = RML_sim$param_est[2:(dim(RML_sim$param_est)[1]),]
    
    ## Update 'obs' ##
    for (j in 1:n_obs){
      for (k in (t_vals_RML+1)){
        obs[[j]][k,] = obs_current[[j]]
      }
    }
    
    ## ----------------------------------
    
    ## Update 'log_liks' ##
    log_liks[t_vals_RML] = RML_sim$ll
    
    ## Update 'obj_funcs' ##
    obj_funcs[t_vals_RML] = RML_sim$obj_func
    
    ## Update 'kf_means' ##
    kf_means[t_vals_RML,] = RML_sim$kf_means
    
    ## ----------------------------------
    
    
    ## Plot ##
    if(RML_plot){
      
      if((0%in%(t_vals%%RML_plot_freq)) | (i==M)){
        
        ## Open plot ##
        if(save_RML_plot){
          if(i==M){
            setwd(fig_wd)
            pdf(RML_plot_filename,width=7,height=7)
          }
        }
      
        ## Number of Plots ##
        par(mar=c(2.8,4.0,1.5,1.5),mgp=c(2.5,1,0))
        if(n_param_bias==0){par(mfrow=c(3,3))}
        if(n_param_bias>0){par(mfrow=c(3,4))}
        
        ## Plot Points ##
        plot_points = seq(1,t_end_RML,RML_plot_point_freq)
        
        ## Signal Parameters ##
        for (k in param_sig_indices){
          
          y_lim = NULL
          y_lab = expression(theta)
          main = NULL
          
          if(!is.null(plot_limits)){
            y_lim = plot_limits[[k]]
          }
          
          if(!is.null(param_sig_names)){
            y_lab = 
              parse(text=paste("hat(",toString(param_sig_names[[k]]),")",sep=""))
            main = param_sig_names[[k]]
          }
          
          plot(plot_points, param_est[plot_points,k], cex=0.4,
               xlim = c(1,n_iterations*t+1), xlab = "t",ylab = y_lab,
               ylim = y_lim, main = main)
          
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
              y_lab = parse(text=paste("hat(",toString(param_obs_names[[k-8]]),")",sep=""))
              main = param_obs_names[[k-8]]
            }
            
            plot(plot_points, param_est[plot_points,k], cex=0.4,
                 xlim = c(1,n_iterations*t+1) ,xlab="t", ylab=y_lab,
                 ylim = y_lim, main = main)
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
          
          if(mean_param){
            if(t_end_RML>(param_est_starts[[k]][1]+1)){
              for (m in 1:(which_start[[k]])){
                abline(v=param_est_starts[[k]][m],col="black",lty=2)
                if(length(param_obs_true$t_vals)==1){
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
        
        ## Bias Parameters ##
        if(n_param_bias>0){
          for (k in param_bias_indices){
            if(k==param_bias_indices[[1]]){
              
              y_lim = NULL
              y_lab = expression(theta)
              main = NULL
              
              if(!is.null(plot_limits)){
                y_lim = plot_limits[[10]]
              }
              
              if(!is.null(param_bias_names)){
                y_lab = parse(text=paste("hat(",toString(param_obs_names[[k-(n_param_sig+n_param_obs)]]),")",sep=""))
                main = param_bias_names[[k-(n_param_sig+n_param_obs)]]
              }
              
              plot(plot_points, param_est[plot_points,k], cex=0.4,
                   xlim = c(1,n_iterations*t+1) ,xlab="t", ylab=expression(theta),
                   ylim = y_lim, main = main)
            }
            
            if(k>param_bias_indices[[1]]){
              points(plot_points,param_est[plot_points,k],pch=19,cex=0.4)
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
            
            if(mean_param){
              if(t_end_RML>(param_est_starts[[k]][1]+1)){
                for (m in 1:(which_start[[k]])){
                  abline(v=param_est_starts[[k]][m],col="black",lty=2)
                  if(length(param_bias_true$t_vals)==1){
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
        }
        
        ## Close Plot ##
        if(save_RML_plot){
          if(i==M){
            dev.off()
           }
        }
      }
    }
    
    ## ----------------------------------
    
    ## Current parameters ##
    param_obs = param_obs_func(param_obs$m, param_obs$m_indices,
                               param_est[t_start_ROSP,param_obs_indices])
    param_sig = param_est[t_start_ROSP,param_sig_indices]
    if(n_param_bias>0){
      param_bias = param_bias_func(param_bias$m,param_bias$m_indices,
                                   param_est[t_start_ROSP,param_bias_indices])
    }
    
    ## Current observation locations ##
    for (j in 1:n_obs){obs_current[[j]] = obs[[j]][t_start_ROSP,]}
    
    ## ----------------------------------
    
    ## Run ROSP for t=t_vals_ROSP ##
    ROSP_sim = RecursOptSens(y = y_sim_ROSP, 
                             n_obs = n_obs, t = t_ROSP, K = K, 
                             param_obs = param_obs, 
                             param_sig = param_sig,  param_bias = param_bias, 
                             obs0 = obs_current, n_iterations = 1, 
                             step_size = step_size_ROSP, W_fourier = W_fourier,
                             spec_weights_mat = spec_weights_mat,
                             plot = FALSE, W_coords = NULL, grad_obs = grad_obs, 
                             plot2d = FALSE, plot_point_freq = NULL, 
                             plot_freq = NULL,  save_plot = FALSE,  
                             filename_plot = NULL, save_plot2d = FALSE,
                             filename_plot2d = NULL, grid = FALSE, 
                             print_iter_freq = ROSP_print_iter_freq,
                             iter_numbers = iter_numbers[t_vals_ROSP],
                             mu0 = mu0, Sigma0 = Sigma0, 
                             mu0_grad_param = mu0_grad_param,
                             Sigma0_grad_param = Sigma0_grad_param,
                             mu0_grad_space = mu0_grad_space,
                             Sigma0_grad_space = Sigma0_grad_space,
                             back_track = ROSP_back_track,
                             back_track_tol = ROSP_back_track_tol, 
                             back_track_independent = ROSP_back_track_independent,
                             back_track_time_independent = ROSP_back_track_time_independent, 
                             back_track_scale = ROSP_back_track_scale,
                             back_track_obs = ROSP_back_track_obs, 
                             dt = dt)
    
    ## ----------------------------------
    
    ## Update step sizes ##
    if(ROSP_back_track){
      step_size_ROSP = rep(list(matrix(0,nrow=t_ROSP,ncol=2)),2)
      for (p in 1:n_obs){
        step_size_ROSP[[p]] = matrix(rep(ROSP_sim$step_size[[p]][t_ROSP,],t),
                                     nrow=t,byrow=T)
      }
    }
    
    ## ----------------------------------
    
    ## Update 'initial values' for KF ##
    mu0 = ROSP_sim$mu0
    Sigma0 = ROSP_sim$Sigma0
    mu0_grad_param = ROSP_sim$mu0_grad_param
    Sigma0_grad_param = ROSP_sim$Sigma0_grad_param
    mu0_grad_space = ROSP_sim$mu0_grad_space
    Sigma0_grad_space = ROSP_sim$Sigma0_grad_space
    
    ## ----------------------------------
    
    ## Update 'param_est' ##
    for (j in (t_vals_ROSP+1)){
      param_est[j,] = param_est[t_start_ROSP,]
    }
    
    ## Update 'obs' ##
    for (j in 1:n_obs){
      obs[[j]][t_vals_ROSP+1,] = ROSP_sim$obs[[j]][2:(dim(ROSP_sim$obs[[j]])[1]),]
    }
    
    ## ----------------------------------
    
    ## Update 'log_liks' ##
    log_liks[t_vals_ROSP] = ROSP_sim$ll
    
    ## Update 'obj_funcs' ##
    obj_funcs[t_vals_ROSP] = ROSP_sim$obj_func
    
    ## Update 'kf_means' ##
    kf_means[t_vals_ROSP,] = ROSP_sim$kf_means
    
    ## ----------------------------------
    
    ## Plot ##
    if(ROSP_plot){
      
      if(0%in%(t_vals%%ROSP_plot_freq) | (i==M)){
        
        if(save_ROSP_plot){
          if(i==M){
            setwd(fig_wd)
            pdf(ROSP_plot_filename,width=7,height=7)
          }
        }
        
        ## Plot parameters ##
        par(mar=c(2.8,3.5,1.5,1.5))
        if(n_obs==1){par(mfrow=c(1,1))}
        if(n_obs==2){par(mfrow=c(1,2))}
        if(n_obs==3){par(mfrow=c(1,3))}
        if(n_obs==4){par(mfrow=c(2,2))}
        if(n_obs==8){par(mfrow=c(2,4))}
        
        ## Plot Points ##
        plot_points = seq(1,t_end_ROSP,ROSP_plot_point_freq)
        
        
        
        ## Plot ##
        for(k in 1:n_obs){
          plot(plot_points,obs[[k]][plot_points,1],cex=.4,ylim=c(0,n),
               xlim=c(0,n_iterations*t),xlab="Time",
               ylab=expression(hat(s)),col=brewer.pal(8,"Greys")[6],
               main=paste("Sensor",k))
          points(plot_points,obs[[k]][plot_points,2],cex=.4,
                 col=brewer.pal(8,"Greys")[8])
          legend(x=0.3*(n_iterations*t),y=3.2,
                 legend=c(expression(hat(s[x])),expression(hat(s[y]))),
                 col=brewer.pal(8,"Greys")[c(6,8)],pch=19,ncol=2)
          
          ## Add dashed lines for weighted coordinates ##
          if(!is.null(W_coords)){
            abline(h=W_coords[[k]][1],lty=2,col=brewer.pal(8,"Reds")[3])
            abline(h=W_coords[[k]][2],lty=2,col=brewer.pal(8,"Reds")[7])
          }
        }
        
        ## Close Plot ##
        if(save_ROSP_plot){
          if(i==M){
            dev.off()
          }
        }
        
      }
    }
    
    ## 2D Plot ##
    if(ROSP_plot2d){
      
      if(0%in%(t_vals%%ROSP_plot_freq) | (i==M)){
        
        if(save_ROSP_plot2d){
          if(i==M){
            setwd(fig_wd)
            pdf(ROSP_plot2d_filename,width=7,height=7)
          }
        }
        
        ## Plot parameters ##
        par(mar=c(4.0,3.5,2.5,1.5),mfrow=c(1,1),mgp=c(2,0.5,0))
        
        ## Plot points ##
        plot_points = seq(1,t_end_ROSP,ROSP_plot_point_freq)
        
        ## Plot initial point ##
        if(!ROSP_plot_initial){
          plot(obs[[1]][1,1],obs[[1]][1,2],xlim = c(0,n),ylim = c(0,n),
               xlab=expression(x),ylab=expression(y), 
               col=colorRampPalette(c("grey10", "grey100"))(100)[1])
        }
        
        if(ROSP_plot_initial){
          plot(obs[[1]][1,1],obs[[1]][1,2],xlim = c(0,n),ylim = c(0,n),
               cex = 2, pch = 21, col=brewer.pal(8,"Greens")[7],
               bg = brewer.pal(8,"Greens")[2], xlab=expression(x),
               ylab=expression(y))
          
          for (k in grad_obs){
            points(obs[[k]][1,1],obs[[k]][1,2], cex = 2, pch = 21, 
                   col=brewer.pal(8,"Greens")[7],
                   bg = brewer.pal(8,"Greens")[2])
          }
        }
        
        ## Add legend ##
        if(ROSP_leg){
          if(length(grad_obs)==n_obs){
            if(!is.null(W_coords)){
              legend(x="bottomright", inset = c(0.035,0.04), pt.cex = c(2,2,2), 
                     pch = c(21,21,21),
                     legend=c("Initial Sensor Locations",
                              "Weighted Observation Locations", 
                              "Final Sensor Locations"),
                     pt.bg = c(brewer.pal(8,"Greens")[2],
                               brewer.pal(8,"Reds")[2],
                               brewer.pal(8,"Blues")[2]),
                     col = c(brewer.pal(8,"Greens")[7],
                             brewer.pal(8,"Reds")[7],
                             brewer.pal(8,"Blues")[7]))
            }
            if(is.null(W_coords)){
              legend(x="bottomright", inset = c(0.035,0.04), pt.cex = c(2,2,2), 
                     pch = c(21,21,21),
                     legend=c("Initial Sensor Locations",
                              "Final Sensor Locations"),
                     pt.bg = c(brewer.pal(8,"Greens")[2],
                               brewer.pal(8,"Blues")[2]),
                     col = c(brewer.pal(8,"Greens")[7],
                             brewer.pal(8,"Blues")[7]))
            }
          }
          if(length(grad_obs)!=n_obs){
            legend(x="bottomright", inset = c(0.035,0.04), pt.cex = c(2,2,2), 
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
          }
        }
        
        ## Add grid ##
        if(ROSP_grid){
          abline(h=0:n,v=0:n,col="gray90",lty=2)
        }
        
        
        ## Highlight weighted coordinates ##
        if(!is.null(W_coords)){
          for(k in 1:length(W_coords)){
            points(W_coords[[k]][1],W_coords[[k]][2],cex=2,
                   pch=21,col=brewer.pal(8,"Purples")[7],
                   bg = brewer.pal(8,"Purples")[2])
          }
        }
        
        ## Plot final point ##
        if(t_end_ROSP==(n_iterations*t)){
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
        }
        
        ## Plot remaining points ##
        for (k in grad_obs){
          for (l in plot_points){
            points(obs[[k]][l,1],obs[[k]][l,2],
                   col=colorRampPalette(c("grey10", "grey100"))(100)[5*k-4],cex=.6)
          }
        }
        
        ## Close Plot ##
        if(save_ROSP_plot2d){
          if(i==M){
            dev.off()
          }
        }
      }
    }
    
    ## ----------------------------------
    
    
  }
  
  ## Output ##
  output = list(param_est = param_est, obs = obs, log_liks = log_liks, 
                obj_funcs = obj_funcs, kf_means = kf_means)
  
  ## Save ##
  if(save_RML && save_ROSP){
    filename = paste(RML_filename,"_",ROSP_filename,sep="")
    save(output,file = filename)
  }
  if(save_RML && !(save_ROSP)){
    save(list(param_est=param_est),file = RML_filename)
  }
  if(save_ROSP && !(save_RML)){
    save(list(obs=obs),file = RML_filename)
  }

  ## Output ##
  return(output)
}

####################

#########################
# JOINT RML & ROSP (v2) #
#########################

## Slightly different implementation than RML_ROSP ##

## -> v1 and v2 are essentially equivalent, and yield the same numerical
##    results; we generally use v1 as it is faster

## -> v2 is a slightly more direct implementation of two-timescale gradient 
##    descent; it updates parameters & sensors using the same observations

## -> v2 still allows one iterate (i.e. the slow timescale) to be held constant 
##    for a number of iterations, but does so by adjusting the appropriate step
##    sizes to zero

## -> In v2, t_RML (resp. t_ROSP) are the proportion of iterations in 
##    max(t_RML,t_ROSP) for which the parameters (resp. sensors) will be 
##    updated. 
##    - E.g. (1): t_RML = 1, t_ROSP = 1 implies parameters and sensors updated 
##      each timestep
##    - E.g. (2): t_RML = 1, t_ROSP = 2 implies parameters updated 1 out of 
##      every 2 timesteps, sensors updated every timestep
##    - E.g. (3): t_RML = 3, t_ROSP = 2 implies parameters updated every
##      timestep, sensors updated 2 out of every 3 timesteps

## -> In v2, RML_ROSP_print_iter_freq (= NULL, integer) is an additional 
##    argument, which represents the interval at which the iterations will be 
##    printed
##      - NB: RML_print_iter_freq and ROSP_print_iter_freq should always be set
##        to NULL. These are maintained as arguments in the function RML_ROSP_v2
##        to keep consistency with the function RML_ROSP

## -> In v2, y_sim_t_dep (= TRUE, FALSE) is an additional argument, which 
##    represents whether or not to simulate new observations as the sensor 
##    locations are updated
##      - NB: y_sim_t_dep = TRUE corresponds to the assumption that the 
##        sensors can actually be moved in real time according to the gradient 
##        updates, and thus we need to simulate new observations as they move.
##      - NB: y_sim_t_dep = FALSE corresponds to the assumption that the
##        sensors can't actually be moved in real time according to the gradient
##        updates, and thus we can simulate all observations at the outset.
## -> NB (1): this has now been added to the function RML_ROSP
## -> NB (2): we assume throughout our numerics that the sensors do move in real
##    time, and set y_sim_t_dep = TRUE. This can actually benefit parameter
##    estimation in some cases (maybe explore this as an additional simulation)?

RML_ROSP_v2 <- function(spde_sim, n, n_obs, t, t_RML, t_ROSP, K, n_param, 
                     param_obs0,param_sig0, param_bias0, log_scale = FALSE, 
                     log_indices,grad_params, hessian_indices, obs0, grad_obs, 
                     n_iterations, step_size_RML, step_size_ROSP, W_fourier,
                     spec_weights_mat = NULL, RML_plot = FALSE, 
                     plot_limits = NULL, param_obs_true = NULL, 
                     param_sig_true = NULL, param_bias_true = NULL, 
                     param_obs_names = NULL, param_sig_names = NULL, 
                     param_bias_names = NULL, mean_param = FALSE, 
                     param_est_starts = NULL, param_est_ends = NULL, 
                     RML_plot_point_freq = 1, RML_plot_freq = 100, 
                     save_RML_plot = FALSE, RML_plot_filename = "RML_plot.pdf", 
                     ROSP_plot = FALSE, ROSP_plot2d = FALSE, W_coords = NULL, 
                     ROSP_plot_point_freq = 1, ROSP_plot_freq = 100,
                     save_ROSP_plot = FALSE, save_ROSP_plot2d = FALSE, 
                     ROSP_plot_filename = "ROSP_plot.pdf",
                     ROSP_plot2d_filename = "ROSP_plot2d.pdf",
                     ROSP_grid = TRUE, save_RML = FALSE, RML_filename = "RML",
                     save_ROSP = FALSE, ROSP_filename = "ROSP", 
                     RML_print_iter_freq = NULL, ROSP_print_iter_freq = NULL,
                     print_iter_freq = NULL, 
                     iter_numbers = 1:(n_iterations*t),
                     mu0 = NULL, Sigma0 = NULL, mu0_grad_param = NULL,
                     mu0_grad_space = NULL, Sigma0_grad_param = NULL, 
                     Sigma0_grad_space = NULL,mu0_grad2_param = NULL, 
                     Sigma0_grad2_param = NULL,
                     RML_back_track = FALSE,
                     RML_back_track_tol = 1e-3, 
                     RML_back_track_independent = TRUE,
                     RML_back_track_time_independent = FALSE, 
                     RML_back_track_scale = 0.5,
                     RML_back_track_param = 1:n_param,
                     ROSP_back_track = FALSE,
                     ROSP_back_track_tol = 1e-3, 
                     ROSP_back_track_independent = TRUE,
                     ROSP_back_track_time_independent = FALSE, 
                     ROSP_back_track_scale = 0.5,
                     ROSP_back_track_obs = 1:n_obs,
                     ROSP_plot_initial = FALSE,
                     ROSP_leg = FALSE,
                     RML_ROSP_print_iter_freq = NULL,
                     y_sim_t_dep = TRUE,
                     dt = 1){
  
  ## ------------------------------------------------
  
  ## Number of parameters ##
  n_param_sig = length(param_sig0)
  n_param_obs = param_obs0$m
  n_param_bias = param_bias0$m
  
  ## Parameter indices ##
  param_sig_indices = 1:n_param_sig
  param_obs_indices = (n_param_sig+1):(n_param_sig+n_param_obs)
  if(n_param_bias>0){param_bias_indices = (n_param_sig+n_param_obs+1):n_param}
  
  ## ------------------------------------------------
  
  ## Matrix of parameter estimates ##
  param_est = matrix(0,ncol=n_param,nrow=t*n_iterations+1)
  
  ## Initial parameter estimate ##
  param0 = c(param_sig0,as.numeric(unlist(param_obs0$tau2_vals)))
  if(param_bias0$m>0){
    param0 = c(param0,as.numeric(unlist(param_bias0$bias_vals)))
  }
  
  ## Rescale if param0 on log scale
  if(log_scale){param0[log_indices] <- exp(param0[log_indices])}
  
  ## Initial parameter estimate ##
  param_est[1,] = param0
  
  ## ------------------------------------------------
  
  ## Vector of log-likelihoods ##
  log_liks = rep(0,t*n_iterations+1)
  
  ## Vector of objective functions ##
  obj_funcs = rep(0,t*n_iterations+1)
  
  ## Matrix of KF means ##
  kf_means = matrix(0,nrow=t*n_iterations+1,ncol=K)
  
  ## ------------------------------------------------
  
  ## Matrices of observation locations ##
  obs = rep(list(matrix(0,ncol=2,nrow=n_iterations*t+1)),n_obs)
  
  ## Initial observation locations ##
  for (i in 1:n_obs){obs[[i]][1,] = obs0[[i]]}
  
  ## ------------------------------------------------
  
  ## Observation simulation ##
  y_sim = matrix(0,nrow=t,ncol=n_obs)
  
  if(!y_sim_t_dep){
    y_sim = y_simulation(spde_sim = spde_sim, t_vals = 1:t, 
                         n_obs = n_obs, obs_loc = obs0, 
                         obs_loc_t_dep = FALSE)
  }
  
  ## ------------------------------------------------
  
  ## Step-size adjustment ##

  ## Time Epochs ##
  M_length = max(t_RML,t_ROSP)
  M = t/M_length
  
  ## Step-size indicator matrix ##
  step_size_indicator = matrix(0,nrow=t,ncol=2)
  for (i in 1:M){
    
    ## Times for current epoch ##
    t_start = (i-1)*M_length+1
    t_end = i*M_length
    t_vals = t_start:t_end
    
    ## RML times in current epoch ##
    t_start_RML = t_start
    t_end_RML = t_start + t_RML - 1
    t_vals_RML = t_start_RML:t_end_RML
    step_size_indicator[t_vals_RML,1] = 1
    
    ## ROSP times in current epoch ##
    t_start_ROSP = t_start
    t_end_ROSP = t_start + t_ROSP - 1
    t_vals_ROSP = t_start_ROSP:t_end_ROSP
    step_size_indicator[t_vals_ROSP,2] = 1
  }
  
  ## Adjust step sizes ##
  if(length(which(step_size_indicator[1:t,1]==0))>0){
    step_size_RML[which(step_size_indicator[1:t,1]==0),] = rep(0,ncol(step_size_RML))
  }
  
  if(length(which(step_size_indicator[1:t,2]==0))>0){
    for (i in 1:n_obs){
      step_size_ROSP[[i]][which(step_size_indicator[1:t,2]==0),] = rep(0,2)
    }
  }
  
  ## ------------------------------------------------
  
  ## SPDE_FT object ##
  SPDE_FT = spde_initialise(n,t,K)
  wave = SPDE_FT$wave
  cosine_indices = SPDE_FT$cosine_indices
  n_cos = SPDE_FT$n_cos
  n_cos_sin = length(cosine_indices)
  K = SPDE_FT$K
  
  ## ------------------------------------------------
  
  ## Initial values ##
  
  if(is.null(mu0)){
    mu0 = rep(0,K)
  }
  if(is.null(Sigma0)){
    Sigma0 = Q_f(param_sig0,wave,n,n_cos=n_cos,nu=1,dt=dt,norm=TRUE,
                 spec_weights_mat = spec_weights_mat)
  }
  if(is.null(mu0_grad_param)){
    mu0_grad_param = rep(list(rep(0,K)),n_param)
  }
  if(is.null(mu0_grad_space)){
    mu0_grad_space = rep(list(rep(0,K)),2*n_obs)
  }
  if(is.null(Sigma0_grad_param)){
    Sigma0_grad_param = Q_grad_f(param_obs0,param_sig0,n,K,wave,n_cos,nu,dt=dt,
                                 norm=TRUE,param_bias0,
                                 spec_weights_mat = spec_weights_mat)
  }
  if(is.null(Sigma0_grad_space)){
    Sigma0_grad_space = Q_grad_space_f(n_obs,K)
  }
  if(is.null(mu0_grad2_param)){
    mu0_grad2_param = rep(list(rep(0,K)),n_param)
  }
  if(is.null(Sigma0_grad2_param)){
    Sigma0_grad2_param = Q_grad2_f(param_obs0,param_sig0,n,K,wave,n_cos,nu,dt=dt,
                                   norm=TRUE,param_bias0,
                                   spec_weights_mat = spec_weights_mat)
  }
  
  ## ------------------------------------------------
  
  ## Joint RML and ROSP ##
  
  for (i in 1:t){
    
    ## ----------------------------------
    
    ## Current parameters ##
    param_obs = param_obs_func(param_obs0$m,param_obs0$m_indices,
                               param_est[i,param_obs_indices])
    param_sig = param_est[i,param_sig_indices]
    
    if(n_param_bias>0){
      param_bias = param_bias_func(param_bias0$m,param_bias0$m_indices,
                                   param_est[i,param_bias_indices])
    }
    else{
      param_bias = param_bias_func()
    }
    
    ## ----------------------------------
    
    ## Current observation locations ##
    obs_current = lapply(1:n_obs,function(j) obs[[j]][i,])
    
    if(y_sim_t_dep){
      coords_RML = obs_current
    }  
    if(!y_sim_t_dep){
      coords_RML = obs0
    }
    
    ## ----------------------------------
    
    ## Current observations ##
    if(y_sim_t_dep){
      y_sim_current = y_simulation(spde_sim = spde_sim, t_vals = i, 
                                   n_obs = n_obs, obs_loc = obs_current, 
                                   obs_loc_t_dep = FALSE)
      y_sim[i,] = y_sim_current
    }
    
    ## RML 'y_sim' ##
    y_sim_RML = matrix(y_sim[i,],nrow=1)
    
    ## ROSP 'y_sim' ##
    y_sim_ROSP = matrix(y_sim[i,],nrow=1)
    
    ## ----------------------------------
    
    ## Run RML ##
    RML_sim = RML(y = y_sim_RML, n = n, t = 1, K = K, 
                  param_obs0 = param_obs, param_sig0 = param_sig, 
                  param_bias0 = param_bias, n_param = n_param, 
                  grad_params = grad_params, 
                  step_size = matrix(step_size_RML[i,],nrow=1), 
                  log_scale = log_scale, log_indices = log_indices, 
                  hessian_indices = hessian_indices, n_iterations = 1, 
                  n_obs = n_obs, coords = coords_RML, plot = FALSE, 
                  plot_limits = NULL, param_obs_true = NULL, 
                  param_sig_true = NULL, param_bias_true = NULL,
                  param_obs_names = NULL, param_sig_names = NULL,
                  param_bias_names = NULL, param_est_starts = NULL, 
                  param_est_ends = NULL, plot_point_freq = NULL, 
                  plot_freq = NULL, save_plot = FALSE, plot_filename = NULL, 
                  save = FALSE, filename = NULL, mean_param = FALSE, 
                  print_iter_freq = RML_print_iter_freq,
                  iter_numbers = iter_numbers[i],
                  mu0 = mu0, Sigma0 = Sigma0, mu0_grad_param = mu0_grad_param,
                  Sigma0_grad_param = Sigma0_grad_param, 
                  mu0_grad_space = mu0_grad_space, 
                  Sigma0_grad_space = Sigma0_grad_space,
                  mu0_grad2_param = mu0_grad2_param,
                  Sigma0_grad2_param = Sigma0_grad2_param,
                  back_track = RML_back_track,
                  back_track_tol = RML_back_track_tol, 
                  back_track_independent = RML_back_track_independent,
                  back_track_time_independent = RML_back_track_time_independent, 
                  back_track_scale = RML_back_track_scale,
                  back_track_param = RML_back_track_param,
                  W_fourier = W_fourier,spec_weights_mat = spec_weights_mat, 
                  dt = dt)
    
    ## ----------------------------------
    
    ## Run ROSP ##
    ROSP_sim = RecursOptSens(y = y_sim_ROSP, 
                             n_obs = n_obs, t = 1, K = K, 
                             param_obs = param_obs, 
                             param_sig = param_sig,  param_bias = param_bias, 
                             obs0 = obs_current, 
                             n_iterations = 1, 
                             step_size = step_size_ROSP,
                             lapply(step_size_ROSP,function(A) matrix(A[i,],nrow=1)), 
                             W_fourier = W_fourier,
                             spec_weights_mat = spec_weights_mat,
                             plot = FALSE, W_coords = NULL, grad_obs = grad_obs, 
                             plot2d = FALSE, plot_point_freq = NULL, 
                             plot_freq = NULL,  save_plot = FALSE,  
                             filename_plot = NULL, save_plot2d = FALSE,
                             filename_plot2d = NULL, save = FALSE,
                             grid = FALSE, 
                             print_iter_freq = ROSP_print_iter_freq,
                             iter_numbers = iter_numbers[i],
                             mu0 = mu0, Sigma0 = Sigma0, 
                             mu0_grad_param = mu0_grad_param,
                             Sigma0_grad_param = Sigma0_grad_param,
                             mu0_grad_space = mu0_grad_space,
                             Sigma0_grad_space = Sigma0_grad_space,
                             back_track = ROSP_back_track,
                             back_track_tol = ROSP_back_track_tol, 
                             back_track_independent = ROSP_back_track_independent,
                             back_track_time_independent = ROSP_back_track_time_independent, 
                             back_track_scale = ROSP_back_track_scale,
                             back_track_obs = ROSP_back_track_obs, 
                             dt = dt)
  
    ## ----------------------------------
    
    ## Update 'initial values' for KF and tangent KFs ##
    mu0 = RML_sim$mu0
    Sigma0 = RML_sim$Sigma0
    mu0_grad_param = RML_sim$mu0_grad_param
    Sigma0_grad_param = RML_sim$Sigma0_grad_param
    mu0_grad_space = ROSP_sim$mu0_grad_space
    Sigma0_grad_space = ROSP_sim$Sigma0_grad_space
    if(!is.null(hessian_indices)){
      mu0_grad2_param = RML_sim$mu0_grad2_param
      Sigma0_grad2_param = RML_sim$Sigma0_grad2_param
    }
    
    ## ----------------------------------
    
    ## Update 'param_est' ##
    param_est[i+1,] = RML_sim$param_est[2:(dim(RML_sim$param_est)[1]),]
    
    ## Update 'obs' ##
    for (j in 1:n_obs){
      obs[[j]][i+1,] = ROSP_sim$obs[[j]][2:(dim(ROSP_sim$obs[[j]])[1]),]
    }
    
    ## ----------------------------------
    
    ## Update 'log_liks' ##
    log_liks[i] = RML_sim$ll
    
    ## Update 'obj_funcs' ##
    obj_funcs[i] = ROSP_sim$obj_func
    
    ## Update 'kf_means' ##
    kf_means[i,] = RML_sim$kf_means
    
    ## ----------------------------------
    
    ## Plots ##
    if(RML_plot){
      
      if((0%in%(i%%RML_plot_freq)) | (i==t)){
        
        ## Open plot ##
        if(save_RML_plot){
          if(i==t){
            setwd(fig_wd)
            pdf(RML_plot_filename,width=7,height=7)
          }
        }
        
        ## Number of Plots ##
        par(mar=c(2.8,4.0,1.5,1.5),mgp=c(2.5,1,0))
        if(n_param_bias==0){par(mfrow=c(3,3))}
        if(n_param_bias>0){par(mfrow=c(3,4))}
        
        ## Plot Points ##
        plot_points = seq(1,i,RML_plot_point_freq)
        
        ## Signal Parameters ##
        for (k in param_sig_indices){
          
          y_lim = NULL
          y_lab = expression(theta)
          main = NULL
          
          if(!is.null(plot_limits)){
            y_lim = plot_limits[[k]]
          }
          
          if(!is.null(param_sig_names)){
            y_lab = 
              parse(text=paste("hat(",toString(param_sig_names[[k]]),")",sep=""))
            main = param_sig_names[[k]]
          }
          
          plot(plot_points, param_est[plot_points,k], cex=0.4,
               xlim = c(1,n_iterations*t+1), xlab = "t",ylab = y_lab,
               ylim = y_lim, main = main)
          
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
            if(i>(param_est_starts[[k]][1]+1)){
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
              y_lab = parse(text=paste("hat(",toString(param_obs_names[[k-8]]),")",sep=""))
              main = param_obs_names[[k-8]]
            }
            
            plot(plot_points, param_est[plot_points,k], cex=0.4,
                 xlim = c(1,n_iterations*t+1) ,xlab="t", ylab=y_lab,
                 ylim = y_lim, main = main)
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
          
          if(mean_param){
            if(i>(param_est_starts[[k]][1]+1)){
              for (m in 1:(which_start[[k]])){
                abline(v=param_est_starts[[k]][m],col="black",lty=2)
                if(length(param_obs_true$t_vals)==1){
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
        
        ## Bias Parameters ##
        if(n_param_bias>0){
          for (k in param_bias_indices){
            if(k==param_bias_indices[[1]]){
              
              y_lim = NULL
              y_lab = expression(theta)
              main = NULL
              
              if(!is.null(plot_limits)){
                y_lim = plot_limits[[10]]
              }
              
              if(!is.null(param_bias_names)){
                y_lab = parse(text=paste("hat(",toString(param_obs_names[[k-(n_param_sig+n_param_obs)]]),")",sep=""))
                main = param_bias_names[[k-(n_param_sig+n_param_obs)]]
              }
              
              plot(plot_points, param_est[plot_points,k], cex=0.4,
                   xlim = c(1,n_iterations*t+1) ,xlab="t", ylab=expression(theta),
                   ylim = y_lim, main = main)
            }
            
            if(k>param_bias_indices[[1]]){
              points(plot_points,param_est[plot_points,k],pch=19,cex=0.4)
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
            
            if(mean_param){
              if(i>(param_est_starts[[k]][1]+1)){
                for (m in 1:(which_start[[k]])){
                  abline(v=param_est_starts[[k]][m],col="black",lty=2)
                  if(length(param_bias_true$t_vals)==1){
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
        }
        
        ## Close Plot ##
        if(save_RML_plot){
          if(i==t){
            dev.off()
          }
        }
      }
    }
    
    if(ROSP_plot){
      
      if(0%in%(i%%ROSP_plot_freq) | (i==t)){
        
        if(save_ROSP_plot){
          if(i==t){
            setwd(fig_wd)
            pdf(ROSP_plot_filename,width=7,height=7)
          }
        }
        
        ## Plot parameters ##
        par(mar=c(2.8,3.5,1.5,1.5))
        if(n_obs==1){par(mfrow=c(1,1))}
        if(n_obs==2){par(mfrow=c(1,2))}
        if(n_obs==3){par(mfrow=c(1,3))}
        if(n_obs==4){par(mfrow=c(2,2))}
        if(n_obs==8){par(mfrow=c(2,4))}
        
        ## Plot Points ##
        plot_points = seq(1,i,ROSP_plot_point_freq)
        
        
        
        ## Plot ##
        for(k in 1:n_obs){
          plot(plot_points,obs[[k]][plot_points,1],cex=.4,ylim=c(0,n),
               xlim=c(0,n_iterations*t),xlab="Time",
               ylab=expression(hat(s)),col=brewer.pal(8,"Greys")[6],
               main=paste("Sensor",k))
          points(plot_points,obs[[k]][plot_points,2],cex=.4,
                 col=brewer.pal(8,"Greys")[8])
          legend(x=0.3*(n_iterations*t),y=3.2,
                 legend=c(expression(hat(s[x])),expression(hat(s[y]))),
                 col=brewer.pal(8,"Greys")[c(6,8)],pch=19,ncol=2)
          
          ## Add dashed lines for weighted coordinates ##
          if(!is.null(W_coords)){
            abline(h=W_coords[[k]][1],lty=2,col=brewer.pal(8,"Reds")[3])
            abline(h=W_coords[[k]][2],lty=2,col=brewer.pal(8,"Reds")[7])
          }
        }
        
        ## Close Plot ##
        if(save_ROSP_plot){
          if(i==t){
            dev.off()
          }
        }
        
      }
    }
    
    if(ROSP_plot2d){
      
      if(0%in%(i%%ROSP_plot_freq) | (i==t)){
        
        if(save_ROSP_plot2d){
          if(i==t){
            setwd(fig_wd)
            pdf(ROSP_plot2d_filename,width=7,height=7)
          }
        }
        
        ## Plot parameters ##
        par(mar=c(4.0,3.5,2.5,1.5),mfrow=c(1,1),mgp=c(2,0.5,0))
        
        ## Plot points ##
        plot_points = seq(1,i,ROSP_plot_point_freq)
        
        ## Plot initial point ##
        if(!ROSP_plot_initial){
          plot(obs[[1]][1,1],obs[[1]][1,2],xlim = c(0,n),ylim = c(0,n),
               xlab=expression(x),ylab=expression(y), 
               col=colorRampPalette(c("grey10", "grey100"))(100)[1])
        }
        
        if(ROSP_plot_initial){
          plot(obs[[1]][1,1],obs[[1]][1,2],xlim = c(0,n),ylim = c(0,n),
               cex = 2, pch = 21, col=brewer.pal(8,"Greens")[7],
               bg = brewer.pal(8,"Greens")[2], xlab=expression(x),
               ylab=expression(y))
          
          for (k in grad_obs){
            points(obs[[k]][1,1],obs[[k]][1,2], cex = 2, pch = 21, 
                   col=brewer.pal(8,"Greens")[7],
                   bg = brewer.pal(8,"Greens")[2])
          }
        }
        
        ## Add legend ##
        if(ROSP_leg){
          if(length(grad_obs)==n_obs){
            if(!is.null(W_coords)){
              legend(x="bottomright", inset = c(0.035,0.04), pt.cex = c(2,2,2), 
                     pch = c(21,21,21),
                     legend=c("Initial Sensor Locations",
                              "Weighted Observation Locations", 
                              "Final Sensor Locations"),
                     pt.bg = c(brewer.pal(8,"Greens")[2],
                               brewer.pal(8,"Reds")[2],
                               brewer.pal(8,"Blues")[2]),
                     col = c(brewer.pal(8,"Greens")[7],
                             brewer.pal(8,"Reds")[7],
                             brewer.pal(8,"Blues")[7]))
            }
            if(is.null(W_coords)){
              legend(x="bottomright", inset = c(0.035,0.04), pt.cex = c(2,2,2), 
                     pch = c(21,21,21),
                     legend=c("Initial Sensor Locations",
                              "Final Sensor Locations"),
                     pt.bg = c(brewer.pal(8,"Greens")[2],
                               brewer.pal(8,"Blues")[2]),
                     col = c(brewer.pal(8,"Greens")[7],
                             brewer.pal(8,"Blues")[7]))
            }
          }
          if(length(grad_obs)!=n_obs){
            legend(x="bottomright", inset = c(0.035,0.04), pt.cex = c(2,2,2), 
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
          }
        }
        
        ## Add grid ##
        if(ROSP_grid){
          abline(h=0:n,v=0:n,col="gray90",lty=2)
        }
        
        
        ## Highlight weighted coordinates ##
        if(!is.null(W_coords)){
          for(k in 1:length(W_coords)){
            points(W_coords[[k]][1],W_coords[[k]][2],cex=2,
                   pch=21,col=brewer.pal(8,"Purples")[7],
                   bg = brewer.pal(8,"Purples")[2])
          }
        }
        
        ## Plot final point ##
        if(i==(n_iterations*t)){
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
        }
        
        ## Plot remaining points ##
        for (k in grad_obs){
          for (l in plot_points){
            points(obs[[k]][l,1],obs[[k]][l,2],
                   col=colorRampPalette(c("grey10", "grey100"))(100)[5*k-4],cex=.6)
          }
        }
        
        ## Close Plot ##
        if(save_ROSP_plot2d){
          if(i==t){
            dev.off()
          }
        }
      }
    }
    
    ## ----------------------------------
    
    ## Print current iteration ##
    if(!is.null(RML_ROSP_print_iter_freq)){
      if(i%%RML_ROSP_print_iter_freq==0){
        cat("Iteration ", i, "\n")
      }
    }
    
    ## ----------------------------------
    
  }
  
  ## Output ##
  output = list(param_est = param_est, obs = obs, log_liks = log_liks, 
                obj_funcs = obj_funcs, kf_means = kf_means)
  
  ## Save ##
  if(save_RML && save_ROSP){
    filename = paste(RML_filename,"_",ROSP_filename,sep="")
    save(output,file = filename)
  }
  if(save_RML && !(save_ROSP)){
    save(list(param_est=param_est),file = RML_filename)
  }
  if(save_ROSP && !(save_RML)){
    save(list(obs=obs),file = RML_filename)
  }
  
  ## Output ##
  return(output)
}

#########################

