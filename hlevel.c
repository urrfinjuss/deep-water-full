#include "header.h"

/*void compute_rhs(fftwl_complex *inQ, fftwl_complex *inV, fftwl_complex *outQ, fftwl_complex *outV) {
  long double overN = 1.L/state.number_modes;
}*/

void project(fftwl_complex *in, fftwl_complex *out) {
  long double 		overN = 1.L/state.number_modes;
  fftwl_complex 	b0 = 0.L;	

  memcpy(tmpc[0], in, state.number_modes*sizeof(fftwl_complex));
  fftwl_execute(ift0);
  // first we need to compute the proper value of the constant:  
  for (long int j = state.number_modes/2 - 2; j > -1; j--) {
    b0 += 0.5L*(tmpc[0][state.number_modes-j-1]*conjl(conf.w[j]) - tmpc[0][j+1]*conf.w[j]);
  }
  memset(tmpc[0]+state.number_modes/2, 0, state.number_modes/2*sizeof(fftwl_complex));
  for (long int j = 1; j < state.number_modes/2; j++) {
    tmpc[0][j] = tmpc[0][j]*overN;
  }
  tmpc[0][0] = 0.5L*tmpc[0][0]*overN;
  fftwl_execute(ft0);
  for (long int j = 0; j < state.number_modes; j++) {
    out[j] = tmpc[0][j] + b0;
  } 
  //printf("shift_re = %.16LE\nshift_im = %.16LE\n", creall(b0), cimagl(b0));
}

void restore_potential(fftwl_complex *inQ, fftwl_complex *inV, fftwl_complex *out){
  long double overN = 1.L/state.number_modes;
  long double K = 0.L, S = 0.L;
  fftwl_complex z0 = 0.L;
  long int k = 0;

  convertQtoZ(inQ, tmpc[5]); // tmpc[5] <- z-tilde
  for (long int j = 0; j < state.number_modes; j++) {
    tmpc[0][j] = inQ[j]*inQ[j]*overN;
    tmpc[1][j] = -1.IL*inV[j]*overN;
    tmpc[2][j] = tmpc[5][j]*overN;
  }  
  fftwl_execute(ift0);
  fftwl_execute(ift1);
  fftwl_execute(ift2);
  memcpy(tmpc[5], tmpc[2], state.number_modes*sizeof(fftwl_complex)); 
  linear_solve(tmpc[0], tmpc[1], tmpc[2]);
  memcpy(tmpc[1], tmpc[2], state.number_modes*sizeof(fftwl_complex));
  div_jacobian(tmpc[1], tmpc[0]);
  for (long int j = 1; j < state.number_modes/2; j++) {
    tmpc[0][j] = 1.0IL*tmpc[0][j]/j;
  }
  memset(tmpc[0]+state.number_modes/2, 0, state.number_modes*sizeof(fftwl_complex)/2);
  compute_zero_mode_complex(tmpc[0], 0.L, &z0);
  tmpc[0][0] = z0;
  for (long int j = state.number_modes/2-1; j > 0; j--) {
    k = j + state.number_modes/2;
    K += j*tmpc[0][j]*conjl(tmpc[0][j]);
    S += 1.0L/cabsl(inQ[j]*inQ[j])/conf.dq[j];
    S += 1.0L/cabsl(inQ[k]*inQ[k])/conf.dq[k];
  }
  S = 0.L;
  for (long int j = 0; j < state.number_modes; j++) {
    S += 1.0L/cabsl(inQ[j]*inQ[j])/conf.dq[j];
  }
  //state.surfaceE = state.surface_tension*S*overN;
  fftwl_execute(ft0);
  memcpy(out, tmpc[0], state.number_modes*sizeof(fftwl_complex));
}

void get_hamiltonian(fftwl_complex *inQ, fftwl_complex *inV) {
  long double 	overN = 1.L/state.number_modes;
 
  /*  -------  Get Kinetic Energy --------- */ 
  state.kineticE = 0.L;
  for (long int j = 0; j < state.number_modes; j++) {
    tmpc[1][j] = -1.IL*inV[j]*overN/(inQ[j]*inQ[j]*conf.dq[j]);
  }  
  /* --------------------------------------  */
  /* tmpc[0] has z_u, and tmpc[1] has \Phi_u */
  /* --------------------------------------  */
  fftwl_execute(ift1);
  for (long int j = state.number_modes-1; j > 0; j--) {
    state.kineticE += tmpc[1][j]*conjl(tmpc[1][j])/j;
  }
  state.kineticE = 0.5L*PI*state.kineticE;
  /* -------- Get Potential Gravity Energy --------- */
  state.potentialE = 0.L;
  convertQtoZ(inQ, tmpc[5]); // tmpc[5] <- z-tilde
  for (long int j = 0; j < state.number_modes; j++) {
    tmpc[0][j] = cimagl(tmpc[5][j])*cimagl(tmpc[5][j])*overN;
    tmpc[1][j] = cimagl(tmpc[5][j])*overN;
    tmpc[2][j] = 1.L*overN/conf.dq[j];
  }
  fftwl_execute(ift0);
  fftwl_execute(ift1);
  fftwl_execute(ift2);
  for (long int j = state.number_modes/2 - 1; j > 0; j--) {
    state.potentialE += 2.L*creall(tmpc[0][j]*conjl(tmpc[2][j]));
    state.potentialE += 2.L*creall(j*tmpc[0][j]*conjl(tmpc[1][j]));
  }
  state.potentialE += creall(tmpc[0][0]*tmpc[2][0]);
  state.potentialE = PI*state.gravity*state.potentialE;
  /* --------- Get Potential Surface Energy --------- */
  state.surfaceE = 0.L;
  for (long int j = 0; j < state.number_modes; j++) {
    state.surfaceE += overN/cabsl(inQ[j]*inQ[j])/conf.dq[j];
  }
  state.surfaceE = 2.L*PI*state.surface_tension*state.surfaceE;

}

void get_momentum(fftwl_complex *inQ, fftwl_complex *inV) {
  long double 	overN = 1.L/state.number_modes;
  long int 	Np = state.number_modes;
  fftwl_complex	P = 0.L;
  
  for (long int j = 0; j < state.number_modes; j++) {
    tmpc[0][j] =            1.L/(inQ[j]*inQ[j]*conf.dq[j])*overN;
    tmpc[1][j] =  cimagl(inV[j]/(inQ[j]*inQ[j]*conf.dq[j]))*overN;
    //tmpc[2][j] =   -1.IL*conjl(inV[j])/(inQ[j]*inQ[j]*conf.dq[j])*overN;
  }  
  fftwl_execute(ift0);
  fftwl_execute(ift1);
  //fftwl_execute(ift2);
  for (long int j = state.number_modes/2-1; j > 0; j--) {
    P +=  tmpc[0][j]*conjl(tmpc[1][j])/j;
    P += -tmpc[0][Np-j]*conjl(tmpc[1][Np-j])/j;
  }
  state.momentum = 2.L*PI*P;
}

/*
void convertQtoZ(fftwl_complex *in, fftwl_complex *out) {
  // computes Z from Q by series inversion of equation:	   //
  // 		   z_u Q^2 = 1				   //
  unsigned long N = state.number_modes;
  long double overN = 1.L/N;
  long double S0 = 0.L, T0 = 0.L;
  long double P = 0.L, mean_level = 0.L;
  fftwl_complex z0 = 0.L;

  for (long int j = 0; j < N; j++) {
    tmpc[0][j] = in[j]*in[j]*overN;
  }
  fftwl_execute(ift0);
  inverse(tmpc[0],tmpc[1]);
  fftwl_execute(ft1);
  for (long int j = 0; j < N; j++) {
    tmpc[1][j] = (tmpc[1][j] - 1.L)*overN;
  }
  fftwl_execute(ift1);
  div_jacobian(tmpc[1], tmpc[0]);
  for (long int j = N/2-1; j > 0; j--) {
    tmpc[0][j] =  1.0IL*tmpc[0][j]/j;	// this stores z_k, k = 1, ..., N/2
    S0 += -0.5L*j*creall(tmpc[0][j]*conjl(tmpc[0][j])); 
    T0 += -creall(tmpc[0][j]);
  }
  compute_zero_mode_complex(tmpc[0], 1.IL*S0, &z0);
  tmpc[0][0] = T0 - creall(z0) + z0;
  memset(tmpc[0]+N/2, 0, N/2*sizeof(fftwl_complex));
  tmpc[1][0] = cimagl(z0);
  for (long int j = 1; j < N/2; j++) {
    tmpc[1][j] = -0.5IL*tmpc[0][j];
    tmpc[1][N-j] = conjl(tmpc[1][j]);
  }
  div_jacobian(tmpc[1], tmpc[4]); // now S
  fftwl_execute(ft0);
  memcpy(out, tmpc[0], N*sizeof(fftwl_complex));
  for (long int j = 0; j < N; j++) {
    tmpc[2][j] = cpowl(cimagl(tmpc[0][j]),2)*overN;
    tmpc[0][j] = cimagl(tmpc[0][j])*overN;
  }
  fftwl_execute(ift0);
  fftwl_execute(ift2);
  for (long int j = N/2-1; j > 0; j--) {
    mean_level += 2.0L*j*creall(tmpc[0][j]*conjl(tmpc[0][j]));
    P += 2.L*creall((tmpc[4][j] + j*tmpc[2][j])*conjl(tmpc[1][j]));
  }
  mean_level = 2.L*PI*(S0 + mean_level);
  P += creall(tmpc[4][0]*conjl(tmpc[1][0]));
  //state.potentialE = 2.L*PI*state.gravity*P;
  state.potentialE = 0.5L*state.gravity*P;
  state.mean_level = mean_level;
} */


void convertQtoZ(fftwl_complex *in, fftwl_complex *out) {
  // computes Z from Q by series inversion of equation:	   //
  // 		   z_u Q^2 = 1				   //
  unsigned long N = state.number_modes;
  long double overN = 1.L/N;
  long double S0 = 0.L, T0 = 0.L;
  fftwl_complex z0 = 0.L;

  for (long int j = 0; j < N; j++) {
    tmpc[0][j] = in[j]*in[j]*overN;
  }
  fftwl_execute(ift0);
  inverse(tmpc[0],tmpc[1]);
  fftwl_execute(ft1);
  for (long int j = 0; j < N; j++) {
    tmpc[1][j] = (tmpc[1][j] - 1.L)*overN;
  }
  fftwl_execute(ift1);
  div_jacobian(tmpc[1], tmpc[0]);
  for (long int j = N/2-1; j > 0; j--) {
    tmpc[0][j] =  1.0IL*tmpc[0][j]/j;	// this stores z_k, k = 1, ..., N/2
    S0 += -0.5L*j*creall(tmpc[0][j]*conjl(tmpc[0][j])); 
    T0 += -creall(tmpc[0][j]);
  }
  compute_zero_mode_complex(tmpc[0], 1.IL*S0, &z0);
  tmpc[0][0] = T0 + 1.IL*cimagl(z0);
  memset(tmpc[0]+N/2, 0, N/2*sizeof(fftwl_complex));
  fftwl_execute(ft0);
  memcpy(out, tmpc[0], N*sizeof(fftwl_complex));
}

void convertZtoQ(fftwl_complex *in, fftwl_complex *out) {
  // computes Q from tilde-Z (z) by series inversion of:   //
  // 		   (z_q q_u + 1) Q^2 = 1		   //
  // 							   //
  // in		-- array with Z-tilde (z)		   //
  // out	-- array with Q				   //
  unsigned long N = state.number_modes;
  long double overN = 1.L/N;

  memcpy(tmpc[0], in, N*sizeof(fftwl_complex));
  fftwl_execute(ift0);
  for (long int j = 0; j < N/2; j++) {
    tmpc[0][j] =  -1.0IL*((fftwl_complex)(tmpc[0][j]*j))*overN;
    tmpc[0][N-1-j] = (fftwl_complex) 0.0L;
  }
  memset(tmpc[0]+N/2, 0, N/2*sizeof(fftwl_complex));
  fftwl_execute(ft0); 
  for (long int j = 0; j < N; j++) {
    tmpc[0][j] = (tmpc[0][j]*conf.dq[j] + 1.L)*overN; 
  }
  fftwl_execute(ift0);
  inverse(tmpc[0], tmpc[1]);
  square_ft(tmpc[1], tmpc[0]);
  fftwl_execute(ft0);
  memcpy(out, tmpc[0], N*sizeof(fftwl_complex));
}
