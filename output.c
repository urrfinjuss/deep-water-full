#include "header.h"

void real_array_out(char *fname, long double *in) {
  FILE *fh = fopen(fname,"w");
  long double u, q, overN = 1./state.number_modes;

  fprintf(fh, "# 1. q 2. u 3. Array\n\n");
  for (long int j = 0; j < state.number_modes; j++) {
    q = 2.0L*PI*(j*overN - 0.5L) - conf.origin_offset;
    u = conf.image_offset + 2.L*atan2l(conf.scaling*sinl(0.5L*q), cosl(0.5L*q));
    fprintf(fh, "%.19LE\t%.19LE\t%.19LE\n", q+conf.origin_offset, u, in[j]);
  }
  fclose(fh);
}
void complex_array_out(char *fname, fftwl_complex *in) {
  FILE *fh = fopen(fname,"w");
  long double u, q, overN = 1./state.number_modes;

  fprintf(fh, "# 1. q 2. u 3.-4. Array\n\n");
  for (long int j = 0; j < state.number_modes; j++) {
    q = 2.0L*PI*(j*overN - 0.5L) - conf.origin_offset;
    u = conf.image_offset + 2.L*atan2l(conf.scaling*sinl(0.5L*q), cosl(0.5L*q));
    fprintf(fh, "%.19LE\t%.19LE\t%.19LE\t%.19LE\n", q+conf.origin_offset, u, creall(in[j]), cimagl(in[j]));
  }
  fclose(fh);
}

void print_constants() {
  printf("#----------------------------------------------------------------------------#\n");
  printf("#                                                                            #\n");
  printf("#                               Conformal Map                                #\n");
  printf("#                                Parameters:                                 #\n");
  printf("#                                                                            #\n");
  printf("#                         N  = %8ld                                      #\n", state.number_modes);
  printf("#                         L  = %.17LE                       #\n", conf.scaling);
  printf("#                         q* = %.17LE                       #\n", conf.origin_offset);
  printf("#                         u* = %.17LE                       #\n", conf.image_offset);
  printf("#                                                                            #\n");
  printf("#                                                                            #\n");
  printf("#                                 Constants:                                 #\n");
  printf("#                                                                            #\n");
  printf("#       Px = %.17LE\t  Py = %.17LE       #\n", cimagl(state.momentum), creall(state.momentum));
  printf("#        K = %.17LE\t   P = %.17LE       #\n",state.kineticE/PI, state.potentialE/PI);
  printf("#                                                                            #\n");
  printf("#----------------------------------------------------------------------------#\n");
}


