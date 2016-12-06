#include "header.h"

//static fftwl_complex *coeffQ, *coeffP;
static long double *M;
static fftwl_complex *arrayQ,  *arrayP;
static fftwl_complex **A, **B, **C, **D;
static fftwl_complex *W, *tmp;
static fftwl_complex **Gees, **Cees;
static long double *Ens;

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

void allocate_pade(unsigned long nD) {
  unsigned long 	N = state.number_modes;
  long double 		q, overN = 1.L/N;
  fftwl_complex		xi;
 
  Gees = fftwl_malloc(4*sizeof(fftwl_complex *));
  Cees = fftwl_malloc(4*sizeof(fftwl_complex *));
  Ens = fftwl_malloc(4*sizeof(long double));
  M = fftwl_malloc((N-1)*sizeof(long double));
  W = fftwl_malloc((N-1)*sizeof(fftwl_complex));
  tmp = fftwl_malloc((N-1)*sizeof(fftwl_complex));
  arrayQ = fftwl_malloc((N - 1)*sizeof(fftwl_complex));
  arrayP = fftwl_malloc((N - 1)*sizeof(fftwl_complex));
  A = fftwl_malloc(4*sizeof(fftwl_complex *));
  B = fftwl_malloc(4*sizeof(fftwl_complex *));
  C = fftwl_malloc(4*sizeof(fftwl_complex *));
  D = fftwl_malloc(4*sizeof(fftwl_complex *));
  for (int k = 0; k < 4; k++) {
    Gees[k] = fftwl_malloc((N - 1)*sizeof(fftwl_complex));
    Cees[k] = fftwl_malloc((2*nD + 1)*sizeof(fftwl_complex));
    A[k] = fftwl_malloc((N-1)*sizeof(fftwl_complex));
    B[k] = fftwl_malloc((N-1)*sizeof(fftwl_complex));
    C[k] = fftwl_malloc((N-1)*sizeof(fftwl_complex));
    D[k] = fftwl_malloc((N-1)*sizeof(fftwl_complex));
  }
  for (unsigned long j = 0; j < N-1; j++) arrayQ[j] = 1.0L;
  for (int k = 0; k < nD; k++) {
    xi = -1.0L + 1.0L*cexpl(0.5L*PI*(k + 0.5L)/nD); 
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
     fftwl_free(A[k]);
     fftwl_free(B[k]);
     fftwl_free(C[k]);
     fftwl_free(D[k]);
   }
   fftwl_free(Gees);
   fftwl_free(Cees);
}

void evaluate_poly_array(unsigned long nD) {
  unsigned long N = state.number_modes;
  long double	s, overN = 1.L/N;
  
  for (unsigned long j = 0; j < N-1; j++) {
    s = tanl(0.5L*PI*(2.L*(j+1)*overN - 1.L));
    // --------   step 0
    A[0][j] = 0.L;
    B[0][j] = 1.L;
    // --------   step 1
    A[1][j] = -1.L;
    B[1][j] = -Cees[0][1];
    // --------   step 2
    A[2][j] = s*A[0][j] - Cees[0][2]*A[1][j] - Cees[1][2]*A[0][j];
    B[2][j] = s*B[0][j] - Cees[0][2]*B[1][j] - Cees[1][2]*B[0][j];
    // --------   step 3
    A[3][j] = s*A[1][j] - Cees[0][3]*A[2][j] - Cees[1][3]*A[1][j] - Cees[2][3]*A[0][j];
    B[3][j] = s*B[1][j] - Cees[0][3]*B[2][j] - Cees[1][3]*B[1][j] - Cees[2][3]*B[0][j];
    // --------   step 4+
    for (unsigned long k = 4; k < 2*nD + 1; k++) {
      tmp[j] = s*A[2][j] - Cees[0][k]*A[3][j] - Cees[1][k]*A[2][j] - Cees[2][k]*A[1][j] - Cees[3][k]*A[0][j];	
      A[0][j] = A[1][j];
      A[1][j] = A[2][j];
      A[2][j] = A[3][j];
      A[3][j] = tmp[j]; 

      tmp[j] = s*B[2][j] - Cees[0][k]*B[3][j] - Cees[1][k]*B[2][j] - Cees[2][k]*B[1][j] - Cees[3][k]*B[0][j];
      B[0][j] = B[1][j];
      B[1][j] = B[2][j];
      B[2][j] = B[3][j];
      B[3][j] = tmp[j]; 
    }  
  }
  memcpy(arrayQ, B[3], (N-1)*sizeof(fftwl_complex));
  memcpy(arrayP, A[3], (N-1)*sizeof(fftwl_complex));
}

void fmul_sub(fftwl_complex *in1, fftwl_complex *in2, fftwl_complex cnum, fftwl_complex *out) {
  unsigned long N = state.number_modes;
  for (unsigned long j = 0; j < N-1; j++) {
    out[j] = in1[j] - cnum*in2[j];
  }
}

void sigma_mul(fftwl_complex *in, fftwl_complex *out) {
  unsigned long N = state.number_modes;
  long double	s, overN = 1.L/N;
  for (unsigned long j = 0; j < N-1; j++) {
    s = tanl(0.5L*PI*(2.L*(j+1)*overN - 1.L));
    out[j] = s*in[j];
  }
}

void gram_schmidt(unsigned long nD){
  unsigned long N = state.number_modes;

  memset(Ens, 0, 4*sizeof(long double));
  for (unsigned k = 0; k < 4; k++) {
    memset(Cees[k], 0, (2*nD + 1)*sizeof(fftwl_complex));
  }
  // step 0
  memcpy(Gees[0], W, (N-1)*sizeof(fftwl_complex));
  Ens[0] = 1.0L/dot(Gees[0],Gees[0]);
  // step 1
  for (unsigned long j = 0; j < N-1; j++) tmp[j] = 1.L;
  Cees[0][1] = dot(tmp, Gees[0])*Ens[0];
  fmul_sub(tmp, Gees[0], Cees[0][1], Gees[1]);
  Ens[1] = 1.0L/dot(Gees[1],Gees[1]);
  // step 2
  sigma_mul(Gees[0], tmp);
  Cees[0][2] = dot(tmp, Gees[1])*Ens[1];
  fmul_sub(tmp, Gees[1], Cees[0][2], tmp);
  Cees[1][2] = dot(tmp, Gees[0])*Ens[0];
  fmul_sub(tmp, Gees[0], Cees[1][2], Gees[2]);
  Ens[2] = 1.0L/dot(Gees[2],Gees[2]);
  // step 3
  sigma_mul(Gees[1], tmp);
  Cees[0][3] = dot(tmp, Gees[2])*Ens[2];
  fmul_sub(tmp, Gees[2], Cees[0][3], tmp);
  Cees[1][3] = dot(tmp, Gees[1])*Ens[1];
  fmul_sub(tmp, Gees[1], Cees[1][3], tmp);
  Cees[2][3] = dot(tmp, Gees[0])*Ens[0];
  fmul_sub(tmp, Gees[0], Cees[2][3], Gees[3]);
  Ens[3] = 1.0L/dot(Gees[3],Gees[3]);
  // step 4+
  for (long int j = 4; j < 2*nD + 1; j++) {
      sigma_mul(Gees[2], tmp);
      Cees[0][j] = dot(tmp, Gees[3])*Ens[3];
      fmul_sub(tmp, Gees[3], Cees[0][j], tmp);
      Cees[1][j] = dot(tmp, Gees[2])*Ens[2];
      fmul_sub(tmp, Gees[2], Cees[1][j], tmp);
      Cees[2][j] = dot(tmp, Gees[1])*Ens[1];
      fmul_sub(tmp, Gees[1], Cees[2][j], tmp);
      Cees[3][j] = dot(tmp, Gees[0])*Ens[0];
      fmul_sub(tmp, Gees[0], Cees[3][j], tmp);
   
      memcpy(Gees[0], Gees[1], (N-1)*sizeof(fftwl_complex));
      memcpy(Gees[1], Gees[2], (N-1)*sizeof(fftwl_complex));
      memcpy(Gees[2], Gees[3], (N-1)*sizeof(fftwl_complex));
      memcpy(Gees[3], tmp, (N-1)*sizeof(fftwl_complex));

      Ens[0] = Ens[1];
      Ens[1] = Ens[2];
      Ens[2] = Ens[3];	
      Ens[3] = 1.0L/dot(tmp, tmp);
  }
}

void compute_rational(unsigned long nD) {
  unsigned long N = state.number_modes; 

  allocate_pade(nD);
  pade_real_out("m0.txt", M);
  gram_schmidt(nD);
  evaluate_poly_array(nD);
  set_weight();
  pade_real_out("m1.txt", M);
  gram_schmidt(nD);
  evaluate_poly_array(nD);
  set_weight();
  pade_real_out("m2.txt", M);
  gram_schmidt(nD);
  evaluate_poly_array(nD);
  set_weight();
  pade_real_out("m3.txt", M);
  gram_schmidt(nD);
  evaluate_poly_array(nD);
  set_weight();
  pade_real_out("m4.txt", M);
  gram_schmidt(nD);
  evaluate_poly_array(nD);
  set_weight();
  pade_real_out("m5.txt", M);
  for (unsigned int j = 0; j < N-1; j++){
    tmp[j] = arrayP[j]/arrayQ[j];
  }
  pade_complex_out("rational.txt", tmp);
  deallocate_pade();
  exit(1);
}




