#include "header.h"

void real_array_out(char *fname, long double *in) {
  FILE *fh = fopen(fname,"w");
  long double u, q, overN = 1.L/state.number_modes;

  fprintf(fh, "# 1. q 2. u 3. Array\n\n");
  for (long int j = 0; j < state.number_modes; j++) {
    q = 2.0L*PI*(j*overN - 0.5L) - conf.origin_offset;
    u = conf.image_offset + 2.L*atan2l(conf.scaling*sinl(0.5L*q), cosl(0.5L*q));
    fprintf(fh, "%24.17LE\t%24.17LE\t%24.17LE\n", q+conf.origin_offset, u, in[j]);
  }
  fclose(fh);
}
void complex_array_out(char *fname, fftwl_complex *in) {
  FILE *fh = fopen(fname,"w");
  long double u, q, overN = 1.L/state.number_modes;

  fprintf(fh, "# 1. q 2. u 3.-4. Array\n\n");
  for (long int j = 0; j < state.number_modes; j++) {
    q = 2.0L*PI*(j*overN - 0.5L) - conf.origin_offset;
    u = conf.image_offset + 2.L*atan2l(conf.scaling*sinl(0.5L*q), cosl(0.5L*q));
    fprintf(fh, "%24.17LE\t%24.17LE\t%24.17LE\t%24.17LE\n", q+conf.origin_offset, u, creall(in[j]), cimagl(in[j]));
  }
  fclose(fh);
}

void spec_out(char *fname, fftwl_complex *in1, fftwl_complex *in2) {
  FILE 		*fh = fopen(fname,"w");
  long double 	overN = 1.L/state.number_modes;

  fprintf(fh, "# 1. k 2. |a_k| 3. |b_k|\n\n");
  for (long int j = 0; j < state.number_modes; j++) {
    fprintf(fh, "%ld\t%23.17LE\t%23.17LE\n", j, cabsl(in1[j])*overN, cabsl(in2[j])*overN);
  }
  fclose(fh);
}

void surface_out(char *fname, fftwl_complex *in) {
  FILE *fh = fopen(fname,"w");
  long double u, q, overN = 1.L/state.number_modes;

  fprintf(fh, "# 1. x 2. y\n\n");
  for (long int j = 0; j < state.number_modes; j++) {
    q = 2.0L*PI*(overN*j - 0.5L) - conf.origin_offset;
    u = conf.image_offset + 2.L*atan2l(conf.scaling*sinl(0.5L*q), cosl(0.5L*q));
    fprintf(fh, "%23.17LE\t%23.17LE\n", u + creall(in[j]), cimagl(in[j]));
  }
  fclose(fh);
}

void print_constants() {
  printf("#\t\t\t------------------------------------------------------\t\t\t#\n");
  printf("#\t\t\t\t\t\t\t\t\t\t\t\t#\n");
  printf("#\t\t\tConformal Map\t\t\t\t\t\t\t\t#\n");
  printf("#\t\t\tParameters:\t\t\t\t\t\t\t\t#\n");
  printf("#\t\t\t\t\t\t\t\t\t\t\t\t#\n");
  printf("#\t\t\tFourier modes N:%8ld\t\t\t\t\t\t#\n", state.number_modes);
  printf("#\t\t\tL-scaling:\t%.16LE\t\t\t\t\t#\n", conf.scaling);
  printf("#\t\t\tQ-star:\t\t%18.16LE\t\t\t\t\t#\n", conf.origin_offset);
  printf("#\t\t\tU-star:\t\t%18.16LE\t\t\t\t\t#\n", conf.image_offset);
  printf("#\t\t\t\t\t\t\t\t\t\t\t\t#\n");
  printf("#\t\t\t\t\t\t\t\t\t\t\t\t#\n");
  printf("#\t\t\tConstants:\t\t\t\t\t\t\t\t#\n");
  printf("#\t\t\t\t\t\t\t\t\t\t\t\t#\n");
  printf("#\t\t\tMomentum X:\t\t%23.16LE\t\t\t\t#\n", cimagl(state.momentum));
  printf("#\t\t\tMomentum Y:\t\t%23.16LE\t\t\t\t#\n", creall(state.momentum));
  printf("#\t\t\tKinetic Energy:\t\t%23.16LE\t\t\t\t#\n",  state.kineticE/PI );
  printf("#\t\t\tPotential Energy:\t%23.16LE\t\t\t\t#\n", state.potentialE/PI);
  printf("#\t\t\tTotal Energy:\t\t%23.16LE\t\t\t\t#\n", (state.potentialE+state.kineticE)/PI);
  printf("#\t\t\t\t\t\t\t\t\t\t\t\t#\n");
  printf("#\t\t\tMean Level:\t%23.16LE\t\t\t\t\t#\n", state.mean_level);
  printf("#\t\t\t------------------------------------------------------\t\t\t#\n");
}


