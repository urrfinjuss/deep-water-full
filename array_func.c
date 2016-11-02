#include "header.h"

static long double overN;

void init_arrayf() {
  overN = 1.L/state.number_modes;
}

void inverse(fftwl_complex *a, fftwl_complex *x) {
  // inverts a*x = 1 to find x by series inversion
  // a,b -- input inverse Fourier coefficients,
  // x   -- output inverse Fourier coefficients
  x[0] = 1.0L/a[0];
  for (long int j = 1; j < state.number_modes/2; j++) {
    fftwl_complex sum = 0.L;
    for (long int l = j-1; l > -1; l--) {
      sum += a[j-l]*x[l];
    }
    x[j] = - sum/a[0];
  }
  memset(x + state.number_modes/2, 0, (state.number_modes/2)*sizeof(fftwl_complex));
}

void square_ft(fftwl_complex *Z, fftwl_complex *x){
  // inverts x^2 = Z to find x by series inversion
  // Z  --  input inverse Fourier coefficients
  // x  --  output inverse Fourier coefficients
  x[0] = csqrtl(Z[0]);
  x[1] = 0.5L*Z[1]/x[0];
  for (long int j = 2; j < state.number_modes/2; j++) {
    fftwl_complex sum = 0.L;
    for (long int m = 0; m < j-1; m++) {
      sum += x[m+1]*x[j-m-1];
    }
    x[j] = 0.5L*(Z[j] - sum)/x[0];
  }
  memset(x + state.number_modes/2, 0, (state.number_modes/2)*sizeof(fftwl_complex));
}
