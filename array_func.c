#include "header.h"

static fftwl_plan *pb, *pf;
static long double overN;

void init_arrayf() {
  overN = 1.L/state.N;
}

void inverse(fftwl_complex *a, fftwl_complex *b, fftwl_complex *x) {
  // inverts a*x = b to find x by series inversion
  // a,b -- input inverse Fourier coefficients,
  // x   -- output inverse Fourier coefficients
  x[0] = b[0]/a[0];
  for (long int j = 1; j < state.N/2; j++) {
    fftwl_complex sum = 0.L;
    for (long int l = j-1; l > -1; l--) {
      sum += a[j-l]*x[l];
    }
    x[j] = (b[j] - sum/a[0]);
  }
  memset(x + state.N/2, 0, (state.N/2)*sizeof(fftwl_complex));
}
