#include <R.h>
#include <Rmath.h>
#include <stdlib.h>

void RML_spectral(double y_FFT[], double mu_i_i1[], double mu_i_i[], double Sigma_i_i1[], double Sigma_i_i[], double mu_i_i1_grad[],double mu_i_i_grad[], double Sigma_i_i1_grad[],double Sigma_i_i_grad[], double R[], double G_cos[], double G_1[], double G_2[], double Q[], double R_grad[], double G_grad_cos[], double G_grad_1[], double G_grad_2[], double Q_grad[], int *i, int *t, int *n_cos, int *n_cos_sin, double *tau2, double theta[], int est_indices_index[], int *n_est, int fixed_indices_index[], int *n_theta, double ll_conditional_grad_auxiliary[], int *r, double step_size_vec[]){
  int j,k,K;
  double ll_conditional_grad;
  K = 2 * *n_cos_sin + *n_cos;
  for (k=0; k<*n_cos;k++){
    mu_i_i1[k]=G_cos[k]*mu_i_i[k];
    Sigma_i_i1[k]=Sigma_i_i[k]*(pow(G_cos[k],2))+Q[k];
  }
  for (k=0; k<*n_cos_sin;k++){
    mu_i_i1[*n_cos+2*k]=G_1[k]*mu_i_i[*n_cos+2*k]+G_2[k]*mu_i_i[*n_cos+2*k+1];
    mu_i_i1[*n_cos+2*k+1]=G_1[k]*mu_i_i[*n_cos+2*k+1]-G_2[k]*mu_i_i[*n_cos+2*k];
    Sigma_i_i1[*n_cos+2*k]=Sigma_i_i[*n_cos+2*k]*(pow(G_1[k],2)+pow(G_2[k],2))+Q[*n_cos+2*k];
    Sigma_i_i1[*n_cos+2*k+1]=Sigma_i_i[*n_cos+2*k+1]*(pow(G_1[k],2)+pow(G_2[k],2))+Q[*n_cos+2*k+1];
  }
  for(k=0; k<K;k++){
    mu_i_i[k+K] = mu_i_i1[k] + Sigma_i_i1[k]/(Sigma_i_i1[k]+*tau2)*(y_FFT[k+(*i-1)*K]-mu_i_i1[k]);
    Sigma_i_i[k+K] = Sigma_i_i1[k] - Sigma_i_i1[k]/(Sigma_i_i1[k]+*tau2)*Sigma_i_i1[k];
  }
  for (j=0;j<*n_theta;j++){
    if(est_indices_index[j]==1){
      for (k=0; k<*n_cos;k++){
        mu_i_i1_grad[k+j*K] = G_cos[k]*mu_i_i_grad[k+j*K] + G_grad_cos[k+j**n_cos]*mu_i_i[k];
        Sigma_i_i1_grad[k+j*K] = Sigma_i_i_grad[k+j*K]*(pow(G_cos[k],2)) + Sigma_i_i[k]*G_grad_cos[k+j**n_cos]*(pow(G_cos[k],2))*2 + Q_grad[k+j*K];
      }
      for (k=0; k<*n_cos_sin;k++){
        mu_i_i1_grad[*n_cos+2*k+j*K] = G_1[k]*mu_i_i_grad[*n_cos+2*k+j*K] + G_2[k]*mu_i_i_grad[*n_cos+2*k+1+j*K] + G_grad_1[k+j**n_cos_sin]*mu_i_i[*n_cos+2*k] + G_grad_2[k+j**n_cos_sin]*mu_i_i[*n_cos+2*k+1];
        mu_i_i1_grad[*n_cos+2*k+1+j*K] = G_1[k]*mu_i_i_grad[*n_cos+2*k+1+j*K] - G_2[k]*mu_i_i_grad[*n_cos+2*k+j*K] + G_grad_1[k+j**n_cos_sin]*mu_i_i[*n_cos+2*k+1] - G_grad_2[k+j**n_cos_sin]*mu_i_i[*n_cos+2*k];
        Sigma_i_i1_grad[*n_cos+2*k+j*K] = Sigma_i_i_grad[*n_cos+2*k+j*K]*(pow(G_1[k],2)+pow(G_2[k],2)) + Sigma_i_i[*n_cos+2*k]*2*(G_1[k]*G_grad_1[k+j**n_cos_sin]+G_2[k]*G_grad_2[k+j**n_cos_sin]) + Q_grad[*n_cos+2*k+j*K];
        Sigma_i_i1_grad[*n_cos+2*k+1+j*K] = Sigma_i_i_grad[*n_cos+2*k+1+j*K]*(pow(G_1[k],2)+pow(G_2[k],2)) + Sigma_i_i[*n_cos+2*k+1]*2*(G_1[k]*G_grad_1[k+j**n_cos_sin]+G_2[k]*G_grad_2[k+j**n_cos_sin]) + Q_grad[*n_cos+2*k+1+j*K];
      }
      ll_conditional_grad[j]=0;
      for (k=0; k<K;k++){
        mu_i_i_grad[k+*n_theta*K+j*K] = mu_i_i1_grad[k+j*K] + Sigma_i_i1_grad[k+j*K]/(Sigma_i_i1[k]+*tau2)*(y_FFT[k+(*i-1)*K]-mu_i_i1[k]) - Sigma_i_i1[k]/(Sigma_i_i1[k]+*tau2)*(Sigma_i_i1_grad[k+j*K]+R_grad[j])/(Sigma_i_i1[k]+*tau2)*(y_FFT[k+(*i-1)*K]-mu_i_i1[k]) - Sigma_i_i1[k]/(Sigma_i_i1[k]+*tau2)*mu_i_i1_grad[k+j*K];
        Sigma_i_i_grad[k+*n_theta*K+j*K] = Sigma_i_i1_grad[k+j*K] - Sigma_i_i1_grad[k+j*K]/(Sigma_i_i1[k]+*tau2)*Sigma_i_i1[k] + Sigma_i_i1[k]/(Sigma_i_i1[k]+*tau2)*(Sigma_i_i1_grad[k+j*K]+R_grad[j])/(Sigma_i_i1[k]+*tau2)*Sigma_i_i1[k] - Sigma_i_i1[k]/(Sigma_i_i1[k]+*tau2)*Sigma_i_i1_grad[k+j*K];
        ll_conditional_grad_auxiliary[k] = -0.5*(Sigma_i_i1_grad[k+j*K]+R_grad[j])/(Sigma_i_i1[k]+*tau2) + 0.5*mu_i_i1_grad[k+j*K]/(Sigma_i_i1[k]+*tau2)*(y_FFT[k+(*i-1)*K]-mu_i_i1[k]) + 0.5*(y_FFT[k+(*i-1)*K]-mu_i_i1[k])/(Sigma_i_i1[k]+*tau2)*(Sigma_i_i1_grad[k+j*K]+R_grad[j])/(Sigma_i_i1[k]+*tau2)*(y_FFT[k+(*i-1)*K]-mu_i_i1[k]) + 0.5*(y_FFT[k+(*i-1)*K]-mu_i_i1[k])/(Sigma_i_i1[k]+*tau2)*mu_i_i1_grad[k+j*K];
        ll_conditional_grad = ll_conditional_grad + ll_conditional_grad_auxiliary[k];
      }
      theta[*n_theta*((*r-1)**t+*i)+j] = theta[*n_theta*((*r-1)**t+(*i-1))+j] + step_size_vec[*n_theta*((*r-1)**t+(*i-1))+j]*ll_conditional_grad;
    }
    if (est_indices_index[j]==0){
      theta[*n_theta*((*r-1)**t+*i)+j] = theta[*n_theta*((*r-1)**t+(*i-1))+j];
    }
  }
}


  