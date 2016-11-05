#include "header.h"

#define COMPLEX_ARRAYS 	6	
#define REAL_ARRAYS	2

static const unsigned long int n_complex_arrays = COMPLEX_ARRAYS;
static const unsigned long int n_real_arrays = REAL_ARRAYS;
static fftwl_complex	*aux_array;

fftwl_plan	ft0, ift0;
fftwl_plan	ft1, ift1;
long double 	**tmpr;
fftwl_complex 	**tmpc;
fftwl_complex 	**data;

void init_memory() {
  printf("Requesting %lu complex arrays of length %lu\n", n_complex_arrays, state.number_modes);
  printf("Requesting %lu real arrays of length %lu\n", n_real_arrays, state.number_modes);
  tmpc = fftwl_malloc(n_complex_arrays*sizeof(fftwl_complex *)); 
  tmpr = fftwl_malloc(n_real_arrays*sizeof(long double *));   
  data = fftwl_malloc(2*sizeof(fftwl_complex *)); 
}

void allocate_memory() {
  for (long int j = 0; j < n_complex_arrays; j++) {
    tmpc[j] = fftwl_malloc(state.number_modes*sizeof(fftwl_complex));
    memset(tmpc[j], 0, state.number_modes*sizeof(fftwl_complex));
  }
  for (long int j = 0; j < n_real_arrays; j++) {
    tmpr[j] = fftwl_malloc(state.number_modes*sizeof(long double));
    memset(tmpr[j], 0, state.number_modes*sizeof(long double));
  }
  conf.dq =  fftwl_malloc(state.number_modes*sizeof(long double));
  conf.w =  fftwl_malloc((state.number_modes/2-1)*sizeof(fftwl_complex));
  aux_array = fftwl_malloc(state.number_modes*sizeof(fftwl_complex));
  data[0] = fftwl_malloc(state.number_modes*sizeof(fftwl_complex));
  data[1] = fftwl_malloc(state.number_modes*sizeof(fftwl_complex));
  memset(data[0], 0, state.number_modes*sizeof(fftwl_complex));
  memset(data[1], 0, state.number_modes*sizeof(fftwl_complex));

  ft0  = fftwl_plan_dft_1d(state.number_modes, tmpc[0], tmpc[0], FFTW_FORWARD, FMODE);
  ft1  = fftwl_plan_dft_1d(state.number_modes, tmpc[1], tmpc[1], FFTW_FORWARD, FMODE);
  ift0 = fftwl_plan_dft_1d(state.number_modes, tmpc[0], tmpc[0], FFTW_BACKWARD, FMODE);
  ift1 = fftwl_plan_dft_1d(state.number_modes, tmpc[1], tmpc[1], FFTW_BACKWARD, FMODE);
}


void deallocate_memory() {
  for (long int j = 0; j < n_complex_arrays; j++) fftwl_free(tmpc[j]);
  for (long int j = 0; j < n_real_arrays; j++) fftwl_free(tmpr[j]);	 
  fftwl_free(conf.dq);   fftwl_free(conf.w); 
  fftwl_free(aux_array); fftwl_free(data[0]);  fftwl_free(data[1]);
  fftwl_destroy_plan(ft0);
  fftwl_destroy_plan(ft1);
  fftwl_destroy_plan(ift0);
  fftwl_destroy_plan(ift1);
}

void fft_shift(fftwl_complex *in){
  for (long int j = 0; j < state.number_modes/2; j++) {
    aux_array[j] 		= in[state.number_modes/2+j];
    aux_array[j+state.number_modes/2]	= in[j]; 
  }
  memcpy(in, aux_array, state.number_modes*sizeof(fftwl_complex));
}



