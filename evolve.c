#include "header.h"

void compute_rhs(fftwl_complex *inQ, fftwl_complex *inV, fftwl_complex *outQ, fftwl_complex *outV) {
  unsigned long 	N = state.number_modes;
  long double 		overN = 1.L/N;
  long double 		sigma = state.surface_tension;
  long double 		g = state.gravity;
  fftwl_complex 	w2 = cexpl(-1.IL*conf.origin_offset - 2.L*atanhl(conf.scaling));
  fftwl_complex		w1 = conjl(w2);
  fftwl_complex 	b1U = 0.L, b2U = 0.L;
  //fftwl_complex 	b1B = 0.L, b2B = 0.L;

  complex_array_out("Q.ph.txt", data[0]);
  complex_array_out("V.ph.txt", data[1]);
  memcpy(tmpc[0], inQ, N*sizeof(fftwl_complex));
  memcpy(tmpc[1], inV, N*sizeof(fftwl_complex));
  fftwl_execute(ift0); 
  fftwl_execute(ift1); 
  memset(tmpc[0]+N/2, 0, N/2*sizeof(fftwl_complex));
  memset(tmpc[1]+N/2, 0, N/2*sizeof(fftwl_complex));
  for (long int j = 0; j < N/2; j++) {
    tmpc[0][j] = -1.IL*j*tmpc[0][j]*overN;
    tmpc[1][j] = -1.IL*j*tmpc[1][j]*overN;
  }
  fftwl_execute(ft0);
  fftwl_execute(ft1);
  for (long int j = 0; j < N; j++) {
    tmpc[2][j]= 2.L*creall(inV[j]*conjl(inQ[j]*inQ[j]))*overN;
    tmpc[3][j]= inV[j]*conjl(inV[j])+4.L*sigma*conf.dq[j]*cimagl(tmpc[0][j]*conjl(inQ[j]));
    tmpc[3][j]= tmpc[3][j]*overN;
  }
  //project(tmpc[3], tmpc[4]); 
  //complex_array_out("B.ph.old.txt", tmpc[4]);
  //for (long int j = 0; j < N; j++) {
  //  tmpc[3][j]= tmpc[3][j]*overN;
  //}
  fftwl_execute(ift2);
  fftwl_execute(ift3);
  tmpc[4][0] = 0.L;
  //b2B = tmpc[3][N/2-1];	b1B = tmpc[3][N/2+1];
  b2U = tmpc[2][N/2-1];	b1U = tmpc[2][N/2+1];
  for (long int j = 1; j < N/2 - 1; j++) {
    //b1B = b1B*w1 + tmpc[3][N/2+1+j];
    //b2B = b2B*w2 + tmpc[3][N/2-1-j];
    b1U = b1U*w1 + tmpc[2][N/2+1+j];
    b2U = b2U*w2 + tmpc[2][N/2-1-j];
  }
  memset(tmpc[2]+N/2, 0, N/2*sizeof(fftwl_complex));
  memset(tmpc[3]+N/2, 0, N/2*sizeof(fftwl_complex));
  memset(tmpc[4]+N/2, 0, N/2*sizeof(fftwl_complex));
  tmpc[2][0] = 0.5L*tmpc[2][0];
  //tmpc[3][0] = 0.5L*tmpc[3][0];  
  for (long int j = 0; j < N/2; j++) {
    tmpc[3][j] = -1.IL*j*tmpc[3][j];
    tmpc[4][j] = -1.IL*j*tmpc[2][j];
  }
  fftwl_execute(ft2);
  fftwl_execute(ft3);
  fftwl_execute(ft4);
  for (long int j = 0; j < N; j++) {
    tmpc[2][j] += 0.5L*(b1U*w1 - b2U*w2);
    //tmpc[3][j] += 0.5L*(b1B*w1 - b2B*w2); 	// thats irrelevant for B'
  }
  // Summary:
  // tmpc[0] <--  stores Q'
  // tmpc[1] <--  stores V'
  // tmpc[2] <--  stores U 
  // tmpc[3] <--  stores B'
  // tmpc[4] <--  stores U'
  for (long int j = 0; j < N; j++) {
    outQ[j] = 0.5IL*conf.dq[j]*(2.L*tmpc[0][j]*tmpc[2][j]-tmpc[4][j]*inQ[j]);
    outV[j] = g*(inQ[j]*inQ[j]-1.L)+1.IL*conf.dq[j]*(tmpc[2][j]*tmpc[1][j]-inQ[j]*inQ[j]*tmpc[3][j]);
  }
}
