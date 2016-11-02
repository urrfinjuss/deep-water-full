#include <tgmath.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <fftw3.h>

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

// --------  Functions
// memory.c
extern void allocate_memory();
extern void deallocate_memory();
extern void backup_arrays();
extern void init_memory();

// input.c
extern void load_ascii();
extern void load_parameters(int argc, char *argv[]);
extern void read_input(char *fname);

// mapping.c
extern void set_mapping();

// output.c
extern void real_array_out(char* fname, long double *in);
extern void complex_array_out(char *fname, fftwl_complex *in);





