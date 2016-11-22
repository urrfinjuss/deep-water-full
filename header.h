#include <tgmath.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <fftw3.h>

#define MIN(a,b) ((a) < (b) ? a : b)

#define FMODE FFTW_MEASURE
#define MOVE_MESH 1	// 0 if singularity tracking is off , 1 otherwise
#define PI acosl(-1.0)
// --------  Structures
typedef struct input {
  char restart_name[128];		// restart filename 
  char txt_format[64];			// format of text file (pade, ascii, none)
  char time_march[64];			// time-marching (rk4)
  long double gravity;			// free fall acceleration
  long double surface_tension;		// surface tension
  long double tolerance;		// tolerance for refinement
  long double kineticE;			// kinetic energy
  long double potentialE;		// potential energy
  fftwl_complex momentum;		// momentum P = px + i*py
  unsigned long int refinement_counter;	// refinement counter
  unsigned long int number_poles;	// number of poles
  unsigned long int number_modes;	// number of grid points
} params, *params_ptr;

typedef struct conformal_mapping {
  fftwl_complex *w;			// Fourier multipliers
  long double 	*dq;			// dq/du
  long double 	scaling;		// conformal map scaling factor
  long double 	image_offset;		// accumulation center in u-plane
  long double 	origin_offset;		// accumulation center in q-plane
} map, *map_ptr;

// -------- Global Variables
extern params state;
extern map conf;
extern long double 	**tmpr;
extern fftwl_complex 	**tmpc;
extern fftwl_complex	**data;
extern fftwl_plan 	ft0, ift0;
extern fftwl_plan 	ft1, ift1, ift2;

// --------  Functions
// memory.c
extern void allocate_memory();
extern void deallocate_memory();
extern void remap(map_ptr new_map, unsigned long int N); 
extern void init_memory();
extern void fft_shift(fftwl_complex *in);

// input.c
extern void load_ascii();
extern void set_initial_data();
extern void load_parameters(int argc, char *argv[]);
extern void read_input(char *fname);

// array_func.c
extern void init_arrayf();
extern void div_jacobian(fftwl_complex *in, fftwl_complex *out);
extern void inverse(fftwl_complex *a, fftwl_complex *x);
extern void linear_solve(fftwl_complex *a, fftwl_complex *b, fftwl_complex *x);
extern void square_ft(fftwl_complex *Z, fftwl_complex *x);
extern void compute_zero_mode(fftwl_complex *in, long double S0, long double *out);
extern void compute_zero_mode_complex(fftwl_complex *in, fftwl_complex S0, fftwl_complex *out);

// hlevel.c
extern void project(fftwl_complex *in, fftwl_complex *out);
extern void convertZtoQ(fftwl_complex *in, fftwl_complex *out);
extern void convertQtoZ(fftwl_complex *in, fftwl_complex *out);
extern void restore_potential(fftwl_complex *inQ, fftwl_complex *inV, fftwl_complex *out);

// mapping.c
extern void set_mapping();
extern void map_quality(fftwl_complex *in, unsigned int *QC_pass);

// output.c
extern void real_array_out(char* fname, long double *in);
extern void complex_array_out(char *fname, fftwl_complex *in);
extern void print_constants();





