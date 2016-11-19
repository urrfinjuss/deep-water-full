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
}

void restore_potential(fftwl_complex *inQ, fftwl_complex *inV, fftwl_complex *out){
  long double overN = 1.L/state.number_modes;
  long double x0 = 0.L, y0 = 0.L;
  long double K = 0.L;

  for (long int j = 0; j < state.number_modes; j++) {
    tmpc[0][j] = inQ[j]*inQ[j]*overN;
    tmpc[1][j] = -1.IL*inV[j]*overN;
  }  
  fftwl_execute(ift0);
  fftwl_execute(ift1);
  linear_solve(tmpc[0], tmpc[1], tmpc[2]);
  memcpy(tmpc[1], tmpc[2], state.number_modes*sizeof(fftwl_complex));
  div_jacobian(tmpc[1], tmpc[0]);
  for (long int j = 1; j < state.number_modes/2; j++) {
    tmpc[1][j] = 0.5IL*tmpc[0][j]/j;      // stores phi_k, k = 1, ..., N/2
    tmpc[0][j] =  0.5L*tmpc[0][j]/j;	  // stores psi_k, k = 1, ..., N/2
  }
  memset(tmpc[0]+state.number_modes/2, 0, state.number_modes*sizeof(fftwl_complex)/2);
  memset(tmpc[1]+state.number_modes/2, 0, state.number_modes*sizeof(fftwl_complex)/2);
  compute_zero_mode(tmpc[1], 0, &x0);
  compute_zero_mode(tmpc[0], 0, &y0);
  tmpc[0][0] = x0 + 1.IL*y0;  
  for (long int j = 1; j < state.number_modes/2; j++) {
    tmpc[0][j] = 2.L*tmpc[1][j];
  }
  memset(tmpc[0]+state.number_modes/2, 0, state.number_modes*sizeof(fftwl_complex)/2); // tmpc[0] has FT potential
  for (long int j = state.number_modes/2-1; j > 0; j--) {
    K += cimagl(tmpc[0][j]*conjl(1.IL*j*tmpc[0][j]));
  }
  K = -1.L*PI*K;  
  state.kineticE = K;
  //printf("Kinetic Energy\t\t%.19LE\n", K);
  fftwl_execute(ft0);
  memcpy(out, tmpc[0], state.number_modes*sizeof(fftwl_complex));
}

void convertQtoZ(fftwl_complex *in, fftwl_complex *out) {
  // computes Z from Q by series inversion of equation:	   //
  // 		   z_u Q^2 = 1				   //
  long double overN = 1.L/state.number_modes;
  long double y0 = 0.L, x0 = 0.L;
  long double S0 = 0.L;
  long double P = 0.L;

  for (long int j = 0; j < state.number_modes; j++) {
    tmpc[0][j] = in[j]*in[j]*overN;
  }
  fftwl_execute(ift0);
  inverse(tmpc[0],tmpc[1]);
  fftwl_execute(ft1);
  for (long int j = 0; j < state.number_modes; j++) {
    tmpc[1][j] = (tmpc[1][j] - 1.L)*overN;
  }
  fftwl_execute(ift1);
  div_jacobian(tmpc[1], tmpc[0]);
  for (long int j = state.number_modes/2-1; j > 0; j--) {
    tmpc[1][j] = 0.5IL*tmpc[0][j]/j;			// this stores x_k, k = 1, ..., N/2
    tmpc[0][j] =  0.5L*tmpc[0][j]/j;			// this stores y_k, k = 1, ..., N/2
    S0 += -2.0L*j*creall(tmpc[0][j]*conjl(tmpc[0][j])); 
  }
  compute_zero_mode(tmpc[0], S0, &y0);  // this only takes care of y_0
  compute_zero_mode(tmpc[1],  0, &x0);	// with transformation we also need to
  tmpc[0][0] = x0 + 1.IL*y0;
  for (long int j = 1; j < state.number_modes/2; j++) {
    tmpc[0][j] = 2.L*tmpc[1][j];
  }
  memset(tmpc[0]+state.number_modes/2, 0, state.number_modes/2*sizeof(fftwl_complex));
  tmpc[1][0] = 1.0L*y0;
  for (long int j = 1; j < state.number_modes/2; j++) {
    tmpc[1][j] = -0.5IL*tmpc[0][j];
    tmpc[1][state.number_modes-j] = conjl(tmpc[1][j]);
  }
  div_jacobian(tmpc[1], tmpc[4]); // now S
  fftwl_execute(ft0);
  memcpy(out, tmpc[0], state.number_modes*sizeof(fftwl_complex));
  for (long int j = 0; j < state.number_modes; j++) {
    tmpc[2][j] = cpowl(cimagl(tmpc[0][j]),2)*overN;
  }
  fftwl_execute(ift2);
  for (long int j = state.number_modes/2-1; j > 0; j--) {
    P += 2.L*creall((tmpc[4][j] + j*tmpc[2][j])*conjl(tmpc[1][j]));
  }
  P += creall(tmpc[4][0]*conjl(tmpc[1][0]));
  state.potentialE = P;
  //printf("Potential Energy\t%.19LE\n", P0);
}

void convertZtoQ(fftwl_complex *in, fftwl_complex *out) {
  // computes Q from tilde-Z (z) by series inversion of:   //
  // 		   (z_q q_u + 1) Q^2 = 1		   //
  // 							   //
  // in		-- array with Z-tilde (z)		   //
  // out	-- array with Q				   //
  long double overN = 1.L/state.number_modes;

  memcpy(tmpc[0], in, state.number_modes*sizeof(fftwl_complex));
  fftwl_execute(ift0);
  for (long int j = 0; j < state.number_modes/2; j++) {
    tmpc[0][j] = -1.0IL*j*tmpc[0][j]*overN;
  }
  memset(tmpc[0]+state.number_modes/2, 0, (state.number_modes/2)*sizeof(fftwl_complex));
  fftwl_execute(ft0);
  for (long int j = 0; j < state.number_modes; j++) {
    tmpc[0][j] = (tmpc[0][j]*conf.dq[j] + 1.L)*overN; 
  }
  fftwl_execute(ift0);
  inverse(tmpc[0], tmpc[1]);
  square_ft(tmpc[1], tmpc[0]);
  fftwl_execute(ft0);
  memcpy(out, tmpc[0], state.number_modes*sizeof(fftwl_complex));
}
