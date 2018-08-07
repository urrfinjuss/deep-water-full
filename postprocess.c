#include "header.h"

void minmax_ongrid(long double *in, long *imax, long *imin){
  long double 	mmax = in[0], mmin = in[0]; 
	
  *imax = 0;
  *imin = 0;
  for (long j = 1; j < state.number_modes; j++) { 
    if (mmax < in[j]) {
      mmax = in[j];
      *imax = j;	
    } else if (mmin > in[j]) {
      mmin = in[j];
      *imin = j;      
    }
  }
}

void evaluate_anywhere(fftwl_complex *in, fftwl_complex w, fftwl_complex *f) {
  long int 	N = state.number_modes;
  fftwl_complex eye_p = -1.L*cexpl(1.IL*w);
  fftwl_complex fp, fm, dfp, dfm;

  fm = in[N/2-1];
  fp = in[N/2+1];
  dfm = -1.IL*(N/2-1)*fm;
  dfp =  1.IL*(N/2-1)*fp; 
  for (long j = N/2-2; j > 0; j--) {
    /* get the value at w by summing Fourier series with Horner's method */
    fm = in[  j] + conjl(eye_p)*fm;
    fp = in[N-j] +       eye_p *fp;
    /* get the derivative at w by summing Fourier series with Horner's method */
    dfm = -1.IL*j*in[  j] + conjl(eye_p)*dfm;
    dfp =  1.IL*j*in[N-j] +       eye_p *dfp; 
  }
  fm = conjl(eye_p)*fm + 0.5IL*cimagl(in[0]);
  fp =       eye_p *fp + 0.5IL*cimagl(in[0]);

  dfm = conjl(eye_p)*dfm;
  dfp =        eye_p*dfp;
  /* f[0] and f[1] contains the value of the function and the derivative at w */
  f[0] = fm + fp;
  f[1] = dfm + dfp;
}


