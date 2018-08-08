#include "header.h"

static fftwl_complex z_xtr[4];
static long double init_qmin, init_qmax;

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

void init_find_peaks(fftwl_complex *inZ) {
  long 	imin, imax, N = state.number_modes;
  
  for (long j = 0; j < N; j++) {
    tmpr[0][j] = cimagl(inZ[j]);
  }
  minmax_ongrid(tmpr[0], &imax, &imin);
  init_qmin = PI*(2.L*imin/N - 1.L);
  init_qmax = PI*(2.L*imax/N - 1.L);
}

void find_peak(fftwl_complex *inZ) {
  //long		imin, imax;
  long	 	N = state.number_modes;
  long		itmin, itmax;
  long double 	qmin, qmax;

  for (long j = 0; j < state.number_modes; j++) {
    tmpc[0][j] = cimagl(1.L/(data[0][j]*data[0][j]*N));
    tmpc[1][j] = inZ[j]/N;
    tmpr[0][j] = cimagl(inZ[j]);
  }
  fftwl_execute(ift0);
  fftwl_execute(ift1);
  /* ------  init qmin and qmax once and use from the previous step   ------ */
  //printf("Estimate for max is %23.16Le\t and min is %23.16Le\n", qmax, qmin);
  qmin = newton_method(tmpc[0], init_qmin, 1e-10L, &itmin);
  qmax = newton_method(tmpc[0], init_qmax, 1e-10L, &itmax);
  //printf("Initial guess %23.16Le and converged at %23.16Le\n", init_qmax, qmax);
  init_qmin = qmin;
  init_qmax = qmax; 
  /* ------  end continuous setting of initial guess -----  */


  /* ------  init qmin and qmax with on grid values every n-th step   ------ */
  /*
  minmax_ongrid(tmpr[0], &imax, &imin);
  qmin = PI*(2.L*imin/N - 1.L);
  qmax = PI*(2.L*imax/N - 1.L);
  //printf("Estimate for max is %23.16Le\t and min is %23.16Le\n", qmax, qmin);
  qmin = newton_method(tmpc[0], qmin, 1e-10L, &itmin);
  qmax = newton_method(tmpc[0], qmax, 1e-10L, &itmax);
  */
  /* ------  end ongrid setting of initial guess -----  */ 

  //printf("Minimum %23.16Le converged in %4ld iterations\n", qmin, itmin);
  //printf("Maximum %23.16Le converged in %4ld iterations\n", qmax, itmax);
  
  evaluate_anywhere(tmpc[1], qmax, &z_xtr[0]);
  evaluate_anywhere(tmpc[1], qmin, &z_xtr[2]);
  state.qmax = qmax;
  state.xmax = qmax + creall(z_xtr[0]);
  state.ymax = cimagl(z_xtr[0]);
  //printf("Value at maximum (%23.16LE,%23.16LE)\n", qmax+creall(z_xtr[0]), cimagl(z_xtr[0]));
  //printf("Value at minimum (%23.16LE,%23.16LE)\n", qmin+creall(z_xtr[2]), cimagl(z_xtr[2]));
}
