#include <R.h>
#include <Rmath.h>
#include <stdlib.h>

void propagate_spectral(double xtp1[], double xt[], double G11C[], double G11[], double G12[], int *NFc, int *ns){
  int j;
  for(j=0; j<*ns;j++){
    xtp1[j]=G11C[j] * xt[j];
  }
  for(j=0; j<*NFc;j++){
    xtp1[2 * j+*ns]=G11[j] * xt[2 * j+*ns] + G12[j] * xt[2 * j+*ns+1];
    xtp1[2 * j+*ns + 1]=G11[j] * xt[2 * j+*ns + 1] - G12[j] * xt[2 * j+*ns];
  }
}