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

long double newton_method(fftwl_complex *in, long double u0, long double tol, long *iter_count) {
  long double 	u = u0;
  fftwl_complex f[2];
  *iter_count = 0;

  evaluate_anywhere(in, u, f);
  if (cabsl(f[0]) < tol ) return u;
  else while (1) {
    u += - f[0]/f[1];
    evaluate_anywhere(in, u, f);
    (*iter_count)++;
    if (cabsl(f[0]) < tol ) break;
  }
  return u;
}

void find_peak(fftwl_complex *inZ, fftwl_complex **z_xtr) {
	long		imin, imax, N = state.number_modes;
	long		itmin, itmax;
	long double 	qmin, qmax;

	for (long j = 0; j < state.number_modes; j++) {
	  tmpc[0][j] = cimagl(1.L/(data[0][j]*data[0][j]*state.number_modes));
	  tmpr[0][j] = cimagl(inZ[j]);
	}
	fftwl_execute(ift0);
	minmax_ongrid(tmpr[0], &imin, &imax);
	
	qmin = PI*(2.L*imin/N - 1.L);
	qmax = PI*(2.L*imax/N - 1.L);
	
        qmin = newton_method(tmpc[0], qmin, 1e-10L, &itmin);
        qmax = newton_method(tmpc[0], qmax, 1e-10L, &itmax);
        printf("Find minimum converged in %4ld iterations\n", itmin);
        printf("Find maximum converged in %4ld iterations\n", itmax);

	evaluate_anywhere(inZ, qmax, z_xtr[0]);
	evaluate_anywhere(inZ, qmin, z_xtr[1]);
	printf("Location of minimum (%23.16LE,%23.16LE)\n", creall(z_xtr[1][0]), cimagl(z_xtr[1][0]));
	printf("Location of maximum (%23.16LE,%23.16LE)\n", creall(z_xtr[0][0]), cimagl(z_xtr[0][0]));
}
