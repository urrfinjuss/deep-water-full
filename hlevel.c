#include "header.h"


void convertQtoZ(fftwl_complex *in, fftwl_complex *out) {
  // ----------------------------------------------------- //
  // computes Z from Q by series inversion of equation:	   //
  // 		   z_u Q^2 = 1				   //
  // ----------------------------------------------------- //
  long double overN = 1.L/state.number_modes;

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
  for (long int j = 1; j < state.number_modes; j++) {
    tmpc[0][j] = 1.IL*tmpc[0][j]/j;
  }
  //fftwl_execute(ft0);
  //memcpy(out, tmpc[0], state.number_modes*sizeof(fftwl_complex));
  // tmpc[0] <-- contains the correct y1, y2, ... (but not y0)
  long double y0 = 0.L;
  compute_zero_mode(tmpc[0], &y0);
  //tmpc[0][0] = 1.IL*y0;
  printf("Zero Mode of Z is %.19LE\n", cimagl(tmpc[0][0]));
  fftwl_execute(ft0);
  memcpy(out, tmpc[0], state.number_modes*sizeof(fftwl_complex));
}

void convertZtoQ(fftwl_complex *in, fftwl_complex *out) {
  // ----------------------------------------------------- //
  // computes Q from tilde-Z (z) by series inversion of:   //
  // 		   (z_q q_u + 1) Q^2 = 1		   //
  // 							   //
  // in		-- array with Z-tilde (z)		   //
  // out	-- array with Q				   //
  // ----------------------------------------------------- //
  long double overN = 1.L/state.number_modes;

  //complex_array_out("z-tilde.ph.txt", in);
  memcpy(tmpc[0], in, state.number_modes*sizeof(fftwl_complex));
  fftwl_execute(ift0);
  //complex_array_out("z-tilde.ft.txt", tmpc[0]);  
  for (long int j = 0; j < state.number_modes/2; j++) {
    tmpc[0][j] = -1.0IL*j*tmpc[0][j]*overN;
  }
  //tmpc[0][0] = 0.;
  memset(tmpc[0]+state.number_modes/2, 0, (state.number_modes/2)*sizeof(fftwl_complex));
  fftwl_execute(ft0);
  for (long int j = 0; j < state.number_modes; j++) {
    tmpc[0][j] = (tmpc[0][j]*conf.dq[j] + 1.L)*overN; 
  }
  complex_array_out("zu.ph.txt",tmpc[0]);
  fftwl_execute(ift0);
  complex_array_out("zu.ft.txt",tmpc[0]);
  inverse(tmpc[0], tmpc[1]);
  //square_ft(tmpc[1], tmpc[0]);
  fftwl_execute(ft1);
  complex_array_out("R.ph.txt",tmpc[1]);
  memcpy(out, tmpc[0], state.number_modes*sizeof(fftwl_complex));
}
