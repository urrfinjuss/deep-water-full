#include "header.h"

static params_ptr 	 input;
static aux_ptr		 extra;

void init_output(params_ptr ginput, aux_ptr gextra) {
  input = ginput; extra = gextra;
}

void debug_msg(char* in, int EXITF) {
  FILE *fh = fopen("output.log","a");
  fprintf(fh, "%s", in);
  fclose(fh);
  if (EXITF) exit(1);
}

void primitive_output(char *fname, fftwl_complex *in) {
  FILE *fh = fopen(fname,"w");
  long double u, q;

  fprintf(fh, "# 1. u 2. Array\n\n");
  for (int j = 0; j < input->N; j++) {
     q = 2.L*PI*(1.L*j/input->N - 0.5L);
     u = input->u + 2.L*atan2l(input->L*sinl(0.5L*(q-extra->q)), cosl(0.5L*(q-extra->q)));
     fprintf(fh, "%.19LE\t%+.19LE\t%+.19LE\n", u, creall(in[j]), cimagl(in[j]) );
   }
  fclose(fh);
}

void basic_output(char *fname, fftwl_complex *in1, fftwl_complex *in2, long double time) {
  FILE *fh = fopen(fname,"w");
  long double u, q;

  fprintf(fh, "# 1. u 2.-3. Q 4.-5. V\n# time = %.12LE\n\n", time);
  for (int j = 0; j < input->N; j++) {
     q = 2.L*PI*(1.L*j/input->N - 0.5L);
     u = input->u + 2.L*atan2l(input->L*sinl(0.5L*(q-extra->q)), cosl(0.5L*(q-extra->q)));
     fprintf(fh, "%.19LE\t%.19LE\t%.19LE\t%.19LE\t%.19LE\n", u, creall(in1[j]), cimagl(in1[j]), creall(in2[j]), cimagl(in2[j]) );
   }
  fclose(fh);
}

void spec_output(char *fname, fftwl_complex *in1, fftwl_complex *in2, long double time){
  FILE *fh = fopen(fname,"w");
  fprintf(fh, "# 1. k 2. |Q_k| 3. |V_k|\n# time = %.12LE\n\n", time);
  for (long int j = 0; j < input->N/2; j++) {
     fprintf(fh, "%ld\t%.19LE\t%.19LE\n", -j, cabsl(in1[j])/input->N, cabsl(in2[j])/input->N);
  }
  fclose(fh);
}

