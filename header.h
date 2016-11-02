#include <tgmath.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <fftw3.h>

#define FMODE FFTW_MEASURE
#define EXIT_TRUE 1
#define EXIT_FALSE 0
#define SPEC_OUTPUT 1   // whether to generate all the spectrum files
#define MOVE_MESH 1	// 0 if singularity tracking is off , 1 otherwise
#define PI acosl(-1.0)
#define RK4 1
#define DIRK4 0

typedef struct work {
  fftwl_complex *Q;
  fftwl_complex *V;
} data, *data_ptr;

typedef struct integrals {
  long double P, K, Mx, My, ml, H0, y0;
  long double H, C, time, Xloc;
} consq, *consq_ptr;

typedef struct more_arrays {
  fftwl_complex *U, *dU;
  fftwl_complex *dQ, *dV;
  fftwl_complex *B;
  fftwl_complex *w;
  fftwl_complex b0;
  long double q, *newdQ, *newdU; //*dq;
} aux, *aux_ptr;

typedef struct rk4_arrays {
  fftwl_complex *k0q, *k1q, *k2q, *k3q, *tmpq;
  fftwl_complex *k0v, *k1v, *k2v, *k3v, *tmpv;
  long double dt, time, tshift;
  long double D;					// hyperviscosity coeff  alpha/kc^6
  long int nsteps, nskip;
} rk4_data, *rk4_ptr;

typedef struct dirk4_arrays {
  fftwl_complex *k1q, *k2q, *tmpq1, *tmpq2;
  fftwl_complex *k1v, *k2v, *tmpv1, *tmpv2;
  fftwl_complex *k1q_prev, *k2q_prev;
  fftwl_complex *k1v_prev, *k2v_prev;
  long double dt, time, tshift;
  long double D, fp_tolerance;					
  long int nsteps, nskip, max_iter;
} dirk4_data, *dirk4_ptr;


typedef struct input {
  char rname[128];		// filename with binary initial data and input structure
  long double g, s;		// gravity and surface tension
  unsigned long int N;		// number of grid points
  unsigned long int ascii;	// flags if reading a text file
  long double L, u;		// transformation parameters
  long double tl;		// tolerance for refinement
  long int refN;		// refinement counter
  long int d;			// number of poles
} params, *params_ptr;

// --------  Variables

extern params state;
extern long double 	**tmpr;
extern fftwl_complex 	**tmpc;

// --------  Functions
// memory.c
extern void allocate_memory();
extern void deallocate_memory();
extern void backup_arrays();
extern void init_memory();
