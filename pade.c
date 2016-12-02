#include "header.h"

//static fftwl_complex *coeffQ, *coeffP;
static long double *M;
static fftwl_complex *arrayQ, *arrayP, *W, *tmp;
static fftwl_complex **Gees, **Cees;

void pade_real_out(char *fname, long double *in) {
  unsigned long 	N = state.number_modes;
  long double		q, overN = 1.L/N;
  FILE *fh = fopen(fname, "w");
  fprintf(fh, "# 1. q 2. Array\n\n");
  for (long int j = 0; j < N-1; j++) {
    q = PI*(2.0L*(j+1)*overN - 1.0L);
    fprintf(fh, "%.18LE\t%.18LE\n", q, in[j]);
  }
  fclose(fh);
}

void pade_complex_out(char *fname, fftwl_complex *in) {
  unsigned long 	N = state.number_modes;
  long double		q, overN = 1.L/N;
  FILE *fh = fopen(fname, "w");
  fprintf(fh, "# 1. q 2. re 3. im\n\n");
  for (long int j = 0; j < N-1; j++) {
    q = PI*(2.0L*(j+1)*overN - 1.0L);
    fprintf(fh, "%.18LE\t%.18LE\t%.18LE\n", q, creall(in[j]), cimagl(in[j]));
  }
  fclose(fh);
}

void set_weight() {
  unsigned long 	N = state.number_modes;
  long double		q, overN = 1.L/N;

  for (long int j = 0; j < N-1; j++) {
    q = PI*(2.0L*(j+1)*overN - 1.0L);
    M[j] = 0.5L/powl(cosl(0.5*q), 2)/creall(arrayQ[j]*conjl(arrayQ[j]));
  }
}

void prepare_array(fftwl_complex *in) {
  unsigned long N = state.number_modes;
  for (unsigned long j = 0; j < N-1; j++) {
    W[j] = in[j+1] - in[0];
  }
  set_weight();
}

fftwl_complex dot(fftwl_complex *in1, fftwl_complex *in2) {
   unsigned long N = state.number_modes;
   fftwl_complex  rvalue = 0.L;
   for (unsigned long j = 0; j < N-1; j++) {
      rvalue += in1[j]*conjl(in2[j])*M[j];
   }
   return rvalue;
}

void allocate_pade(unsigned long D) {
  unsigned long 	N = state.number_modes;
  long double 		q, overN = 1.L/N;
  fftwl_complex		xi;
 
  Gees = fftwl_malloc(4*sizeof(fftwl_complex *));
  Cees = fftwl_malloc(4*sizeof(fftwl_complex *));
  M = fftwl_malloc((N-1)*sizeof(long double));
  W = fftwl_malloc((N-1)*sizeof(fftwl_complex));
  tmp = fftwl_malloc((N-1)*sizeof(fftwl_complex));
  arrayQ = fftwl_malloc((N - 1)*sizeof(fftwl_complex));
  arrayP = fftwl_malloc((N - 1)*sizeof(fftwl_complex));
  for (int k = 0; k < 4; k++) {
    Gees[k] = fftwl_malloc((N - 1)*sizeof(fftwl_complex));
    Cees[k] = fftwl_malloc((2*D + 1)*sizeof(fftwl_complex));
  }
  memset(Cees, 0, 4*(2*D+1)*sizeof(fftwl_complex));
  for (unsigned long j = 0; j < N-1; j++) arrayQ[j] = 1.0L;
  for (int k = 0; k < D; k++) {
    xi = -1.0L + 1.0L*cexpl(0.5L*PI*(k + 0.5L)/D); 
    for (unsigned long j = 0; j < N-1; j++) {
      q = PI*(2.0L*(j + 1)*overN - 1.0L);
      arrayQ[j] = arrayQ[j]*(tanl(0.5L*q) + 1.IL*xi);
    }
  }
  prepare_array(data[0]);
  pade_complex_out("Q0.txt", W);
  pade_real_out("Weight.txt", M);
  printf("Dot Prod: Re = %.12LE\tIm = %.12LE\n", creall(dot(W,W)), cimagl(dot(W,W)));
}

void deallocate_pade() {
   fftwl_free(arrayQ); 
   fftwl_free(arrayP);
   fftwl_free(M);
   fftwl_free(W);
   fftwl_free(tmp);
   for (int k = 0; k < 4; k++) {
     fftwl_free(Gees[k]);
     fftwl_free(Cees[k]);
   }
   fftwl_free(Gees);
   fftwl_free(Cees);
}

void gram_schmidt(fftwl_complex *in, unsigned long D){
  unsigned long N = state.number_modes;
  long double 	s, overN = 1.L/N;

  // step 0
  memset(Cees, 0, 4*(2*D + 1)*sizeof(fftwl_complex));
  memcpy(Gees[0], W, (N-1)*sizeof(fftwl_complex));
 
  // step 1
  for (unsigned long j = 0; j < N-1; j++) tmp[j] = 1.L;
  Cees[0][1] = dot(tmp, Gees[0])/dot(Gees[0], Gees[0]);
  for (unsigned long j = 0; j < N-1; j++) Gees[1][j] = 1.L - Cees[0][1]*Gees[0][j];

  // step 2
  for (unsigned long j = 0; j < N-1; j++) {
    s = tanl(0.5L*PI*(2.L*(j+1)*overN - 1.L));
    tmp[j] = s*Gees[0][j]; 
  }
  Cees[0][2] = dot(tmp, Gees[1])/dot(Gees[1],Gees[1]);
  for (unsigned long j = 0; j < N-1; j++) tmp[j] = tmp[j] - Cees[0][2]*Gees[1][j]; 
  
  Cees[1][2] = dot(tmp, Gees[0])/dot(Gees[0],Gees[0]);
  for (unsigned long j = 0; j < N-1; j++) tmp[j] = tmp[j] - Cees[1][2]*Gees[0][j]; 
  
  // step 3+
  for ()



  for (unsigned long j = 0; j < N-1; j++) {
    tmp = 1.L - Cees[]
    s = tanl(0.5L*PI*(2.L*(j+1)*overN - 1.L));
  }
}






