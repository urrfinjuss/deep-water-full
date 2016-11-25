#include "header.h"

static fftwl_complex *coeffQ, *coeffP;
static fftwl_complex *arrayQ, *arrayP;

void set_Q0(unsigned long d) {
  long double 		overN = 1.L/state.number_modes;
  long double 		q  = 0.L;
  fftwl_complex 	xi = 0.L;  
  arrayQ = fftwl_malloc(state.number_modes*sizeof(fftwl_complex));
  arrayP = fftwl_malloc(state.number_modes*sizeof(fftwl_complex));
 

  for (long int j = 0; j < state.number_modes; j++) arrayQ[j] = 1.L; 
  for (long int l = 0; l < d; l++) {
    xi = 2.IL + cexp(2.IL*PI*l/d);
    for (long int j = 0; j < state.number_modes; j++) {
      q = PI*(2.L*j*overN - 1.L);
      arrayQ[j] = arrayQ[j]*(tanl(0.5L*q) - xi);
    }
  }
  arrayQ[0] = 1.L; // irrelevant value on the side
  complex_array_out("padeQ.txt",arrayQ);  
}
