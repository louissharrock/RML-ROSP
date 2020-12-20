#include <R.h>
#include <Rmath.h>
#include <stdlib.h>


void kf_spectral_2(double y_FFT[], double mu_i_i1[], double mu_i_i[]
                     , double Sigma_i_i1[], double Sigma_i_i[]
                     , double G_cos[], double G_1[], double G_2[]
                     , double Q[], int *T, int *n_cos, int *n_cos_sin
                     , double *tau2){
  int t,j,k,K;
  K = 2 * *n_cos_sin + *n_cos;
  for (t=0; t<*T;t++){
    for (k=0; k<*n_cos;k++){
      mu_i_i1[k + t*K]=G_cos[k]*mu_i_i[k+t*K];
      Sigma_i_i1[k + t*K]=Sigma_i_i[k+t*K]*(pow(G_cos[k],2))+Q[k];
    }
    for (k=0; k<*n_cos_sin;k++){
      mu_i_i1[*n_cos+2*k+t*K]=G_1[k]*mu_i_i[*n_cos+2*k+t*K]+G_2[k]*mu_i_i[*n_cos+2*k+1+t*K];
      mu_i_i1[*n_cos+2*k+1+t*K]=G_1[k]*mu_i_i[*n_cos+2*k+1+t*K]-G_2[k]*mu_i_i[*n_cos+2*k+t*K];
      Sigma_i_i1[*n_cos+2*k+t*K]=Sigma_i_i[*n_cos+2*k+t*K]*(pow(G_1[k],2)+pow(G_2[k],2))+Q[*n_cos+2*k];
      Sigma_i_i1[*n_cos+2*k+1+t*K]=Sigma_i_i[*n_cos+2*k+1+t*K]*(pow(G_1[k],2)+pow(G_2[k],2))+Q[*n_cos+2*k+1];
    }
    for(k=0; k<K;k++){
      mu_i_i[k+(t+1)*K] = mu_i_i1[k+t*K] + (Sigma_i_i1[k+t*K]*(y_FFT[k+t*K]-mu_i_i1[k+t*K])/(Sigma_i_i1[k+t*K]+*tau2));
      Sigma_i_i[k+(t+1)*K] = Sigma_i_i1[k+t*K]*(1-Sigma_i_i1[k+t*K]/(Sigma_i_i1[k+t*K]+*tau2));
    }
  }
}

void tkf_spectral(double y_FFT[], double mu_i_i1[], double mu_i_i[]
                    , double Sigma_i_i1[], double Sigma_i_i[]
                    , double mu_i_i1_grad[], double mu_i_i_grad[]
                    , double Sigma_i_i1_grad[], double Sigma_i_i_grad[]
                    , double G_cos[], double G_1[], double G_2[]
                    , double G_grad_cos[], double G_grad_1[], double G_grad_2[]
                    , double Q[], double Q_grad[], int *T, int *n_cos
                    , int *n_cos_sin, double *tau2, double R_grad[]){
  int t,j,k,K;
  K = 2 * *n_cos_sin + *n_cos;
  for (t=0; t<*T;t++){
    for (k=0; k<*n_cos;k++){
      mu_i_i1_grad[k+t*K] = G_cos[k]*mu_i_i_grad[k+t*K] + G_grad_cos[k]*mu_i_i[k+t*K];
      Sigma_i_i1_grad[k+t*K] = Sigma_i_i_grad[k+t*K]*(pow(G_cos[k],2)) + 2*Sigma_i_i[k+t*K]*G_cos[k]*G_grad_cos[k] + Q_grad[k];
    }
    for (k=0; k<*n_cos_sin;k++){
      mu_i_i1_grad[*n_cos+2*k+t*K] = G_1[k]*mu_i_i_grad[*n_cos+2*k+t*K] + G_2[k]*mu_i_i_grad[*n_cos+2*k+1+t*K] + G_grad_1[k]*mu_i_i[*n_cos+2*k+t*K] + G_grad_2[k]*mu_i_i[*n_cos+2*k+1+t*K];
      mu_i_i1_grad[*n_cos+2*k+1+t*K] = G_1[k]*mu_i_i_grad[*n_cos+2*k+1+t*K] - G_2[k]*mu_i_i_grad[*n_cos+2*k+t*K] + G_grad_1[k]*mu_i_i[*n_cos+2*k+1+t*K] - G_grad_2[k]*mu_i_i[*n_cos+2*k+t*K];
      Sigma_i_i1_grad[*n_cos+2*k+t*K] = Sigma_i_i_grad[*n_cos+2*k+t*K]*(pow(G_1[k],2)+pow(G_2[k],2)) + 2*Sigma_i_i[*n_cos+2*k+t*K]*(G_1[k]*G_grad_1[k]+G_2[k]*G_grad_2[k]) + Q_grad[*n_cos+2*k];
      Sigma_i_i1_grad[*n_cos+2*k+1+t*K] = Sigma_i_i_grad[*n_cos+2*k+1+t*K]*(pow(G_1[k],2)+pow(G_2[k],2)) + 2*Sigma_i_i[*n_cos+2*k+1+t*K]*(G_1[k]*G_grad_1[k]+G_2[k]*G_grad_2[k]) + Q_grad[*n_cos+2*k+1];
    }
    for (k=0; k<K;k++){
      mu_i_i_grad[k+(t+1)*K] = mu_i_i1_grad[k+t*K] + (Sigma_i_i1_grad[k+t*K]*(y_FFT[k+t*K]-mu_i_i1[k+t*K]))/(Sigma_i_i1[k+t*K]+*tau2) - (Sigma_i_i1[k+t*K]*(Sigma_i_i1_grad[k+t*K]+R_grad[k])*(y_FFT[k+t*K]-mu_i_i1[k+t*K]))/((Sigma_i_i1[k+t*K]+*tau2)*(Sigma_i_i1[k+t*K]+*tau2)) - (Sigma_i_i1[k+t*K]*mu_i_i1_grad[k+t*K])/(Sigma_i_i1[k+t*K]+*tau2);
      Sigma_i_i_grad[k+(t+1)*K] = Sigma_i_i1_grad[k+t*K] - 2*(Sigma_i_i1_grad[k+t*K]*Sigma_i_i1[k+t*K])/(Sigma_i_i1[k+t*K]+*tau2) + (Sigma_i_i1[k+t*K]*Sigma_i_i1[k+t*K]*(Sigma_i_i1_grad[k+t*K]+R_grad[k]))/((Sigma_i_i1[k+t*K]+*tau2)*(Sigma_i_i1[k+t*K]+*tau2)) ;
    }
  }
}

void ll_grad_spectral(double *ll_grad, double conditional_ll_grad[]
                        , double y_FFT[], double mu_i_i1[]
                        , double Sigma_i_i1[], double mu_i_i1_grad[]
                        , double Sigma_i_i1_grad[], int *T, int *K
                        , double *tau2, double *R_grad){
  int t,j,k;
  *ll_grad=0;
  for(t=0; t<*T;t++){
    conditional_ll_grad[t]=0;
  }
  for(t=0; t<*T;t++){
    for(k=0; k<*K;k++){
      conditional_ll_grad[t] = conditional_ll_grad[t] - 0.5*(Sigma_i_i1_grad[k+t**K]+*R_grad)/(Sigma_i_i1[k+t**K]+*tau2) + (mu_i_i1_grad[k+t**K]*(y_FFT[k+t**K]-mu_i_i1[k+t**K]))/(Sigma_i_i1[k+t**K]+*tau2) + 0.5*((y_FFT[k+t**K]-mu_i_i1[k+t**K])*(y_FFT[k+t**K]-mu_i_i1[k+t**K])*(Sigma_i_i1_grad[k+t**K]+*R_grad))/((Sigma_i_i1[k+t**K]+*tau2)*(Sigma_i_i1[k+t**K]+*tau2));
    }
    *ll_grad = *ll_grad + conditional_ll_grad[t];
  }
}
