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
  long double L, u;		// transformation parameters
  long double tl;		// tolerance for refinement
  long int refN;		// refinement counter
  long int d;			// number of poles
} params, *params_ptr;


// parameters.c
extern void init_parameters(params_ptr ginput, aux_ptr gextra, data_ptr garray);
extern void load_parameters(int argc, char* argv[]);
extern void load_data();

// output.c
extern void init_output(params_ptr ginput, aux_ptr gextra);
extern void debug_msg(char* in, int EXITF);
extern void primitive_output(char *fname, fftwl_complex *in);
extern void basic_output(char *fname, fftwl_complex *in1, fftwl_complex *in2, long double time);
extern void spec_output(char *fname, fftwl_complex *in1, fftwl_complex *in2, long double time);

// llevel.c
extern void init_lowlevel(); 
extern void init_timemarching(); 
extern void initialize_auxiliary_arrays();
extern void initialize_data();
extern void inverseQ(fftwl_complex *in, fftwl_complex *out);
extern void hfilter();
extern void get_integrals();
extern void get_surface(char* fname);
extern void get_spectrum(char* fname);
extern fftwl_complex check_resolution();
extern void compute_aux_arrays(fftwl_complex *inQ, fftwl_complex *inV);
extern void move_mesh(long double uin, long double Lin, long int Nin);
extern void set_aux();
extern void resolution_monitor(long int *iter);
extern void read_pade();
extern void simulate();
extern void find_max_height(fftwl_complex *in, long double *x_peak, long double *y_peak, long double *curv);

// exrk4.c
extern void init_exrk4_module(params_ptr ginput, aux_ptr gextra, data_ptr garray, consq_ptr gmotion);
extern void modify_rk4skip(long double factor);
long double init_rk4();
extern void evolve_rk4();
extern void free_rk4();

// dirk4.c
extern void init_dirk4_module(params_ptr ginput, aux_ptr gextra, data_ptr garray, consq_ptr gmotion);
extern void modify_dirk4skip(long double factor);
long double init_dirk4();
extern void dirk4_step();
extern void evolve_dirk4();
extern void free_dirk4();

