#include "header.h"

#define COMPLEX_ARRAYS 	4	
#define REAL_ARRAYS	2

static unsigned long int n_local_size;
static const unsigned long int n_complex_arrays = COMPLEX_ARRAYS;
static const unsigned long int n_real_arrays = REAL_ARRAYS;
static fftwl_complex 	*bQ, *bV;

fftwl_plan	ft0, ift0;
fftwl_plan	ft1, ift1;
long double 	**tmpr;
fftwl_complex 	**tmpc;

void init_memory() {
  n_local_size =	state.number_modes;
  printf("Requesting %lu complex arrays of length %lu\n", n_complex_arrays, n_local_size);
  printf("Requesting %lu real arrays of length %lu\n", n_real_arrays, n_local_size);
  tmpc = fftwl_malloc(n_complex_arrays*sizeof(fftwl_complex *)); 
  tmpr = fftwl_malloc(n_real_arrays*sizeof(long double *)); 
}

void allocate_memory() {
  for (long int j = 0; j < n_complex_arrays; j++) {
    tmpc[j] = fftwl_malloc(n_local_size*sizeof(fftwl_complex));
    memset(tmpc[j], 0, n_local_size*sizeof(fftwl_complex));
  }
  for (long int j = 0; j < n_real_arrays; j++) {
    tmpr[j] = fftwl_malloc(n_local_size*sizeof(long double));
    memset(tmpr[j], 0, n_local_size*sizeof(long double));
  }
  conf.dq =  fftwl_malloc(n_local_size*sizeof(long double));
  conf.w =  fftwl_malloc((n_local_size/2-1)*sizeof(fftwl_complex));
  ft0  = fftwl_plan_dft_1d(n_local_size, tmpc[0], tmpc[0], FFTW_FORWARD, FMODE);
  ft1  = fftwl_plan_dft_1d(n_local_size, tmpc[1], tmpc[1], FFTW_FORWARD, FMODE);
  ift0 = fftwl_plan_dft_1d(n_local_size, tmpc[0], tmpc[0], FFTW_BACKWARD, FMODE);
  ift1 = fftwl_plan_dft_1d(n_local_size, tmpc[1], tmpc[1], FFTW_BACKWARD, FMODE);
}

void backup_arrays() {
  unsigned long int n_local_previous = n_local_size;
  bQ = fftwl_malloc(n_local_size*sizeof(fftwl_complex));
  bV = fftwl_malloc(n_local_size*sizeof(fftwl_complex));
  memcpy(bQ, tmpc[0], n_local_size*sizeof(fftwl_complex));
  memcpy(bV, tmpc[1], n_local_size*sizeof(fftwl_complex));
  deallocate_memory();
  n_local_size = state.number_modes;
  allocate_memory();
  if (n_local_previous < n_local_size) {
	unsigned long int mem_offset = n_local_size-n_local_previous/2;
	memcpy(tmpc[0], bQ, (n_local_previous/2)*sizeof(fftwl_complex));
	memcpy(tmpc[1], bV, (n_local_previous/2)*sizeof(fftwl_complex));
	memcpy(tmpc[0]+mem_offset, bQ, (n_local_previous/2)*sizeof(fftwl_complex));
	memcpy(tmpc[1]+mem_offset, bV, (n_local_previous/2)*sizeof(fftwl_complex));
  } else if (n_local_previous > n_local_size) {
	unsigned long int mem_offset = n_local_previous - n_local_size/2;
	memcpy(tmpc[0], bQ, (n_local_size/2)*sizeof(fftwl_complex));
        memcpy(tmpc[1], bV, (n_local_size/2)*sizeof(fftwl_complex));
        memcpy(tmpc[0]+mem_offset, bQ, (n_local_size/2)*sizeof(fftwl_complex));
        memcpy(tmpc[1]+mem_offset, bV, (n_local_size/2)*sizeof(fftwl_complex));
  } else {
	memcpy(tmpc[0], bQ, n_local_size*sizeof(fftwl_complex));
	memcpy(tmpc[1], bV, n_local_size*sizeof(fftwl_complex));
  }
}

void deallocate_memory() {
  for (long int j = 0; j < n_complex_arrays; j++) fftwl_free(tmpc[j]);
  for (long int j = 0; j < n_real_arrays; j++) fftwl_free(tmpr[j]);	 
  fftwl_free(conf.dq); fftwl_free(conf.w); 
  fftwl_destroy_plan(ft0);
  fftwl_destroy_plan(ft1);
  fftwl_destroy_plan(ift0);
  fftwl_destroy_plan(ift1);
}

