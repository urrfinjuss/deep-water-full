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

