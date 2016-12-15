#include "header.h"

//static fftwl_complex *coeffQ, *coeffP;
static long double *M, *Ens;
static fftwl_complex *arrayQ,  *arrayP;
static fftwl_complex **A, **B;
//static fftwl_complex **C, **D;
static fftwl_complex *W, *tmp;
static fftwl_complex **Gees, **Cees;

void init_pade(){
  Gees = fftwl_malloc(4*sizeof(fftwl_complex *));
  Cees = fftwl_malloc(4*sizeof(fftwl_complex *));
  Ens = fftwl_malloc(4*sizeof(long double));
  A = fftwl_malloc(4*sizeof(fftwl_complex *));
  B = fftwl_malloc(4*sizeof(fftwl_complex *));
}
void pade_real_out(char *fname, long double *in) {
  unsigned long 	N = state.number_modes;
  long double		q, overN = 1.L/N;
  FILE *fh = fopen(fname, "w");
  fprintf(fh, "# 1. q 2. Array\n\n");
  for (long int j = 0; j < N-1; j++) {
    q = PI*(2.0L*(j+1)*overN - 1.0L);
    fprintf(fh, "%.18LE\t%.18LE\n", q, in[j]);
  }
  fclose(fh);
}

void pade_complex_out(char *fname, fftwl_complex *in) {
  unsigned long 	N = state.number_modes;
  long double		q, overN = 1.L/N;
  FILE *fh = fopen(fname, "w");
  fprintf(fh, "# 1. q 2. re 3. im\n\n");
  for (long int j = 0; j < N-1; j++) {
    q = PI*(2.0L*(j+1.L)*overN - 1.0L);
    fprintf(fh, "%.18LE\t%.18LE\t%.18LE\n", q, creall(in[j]), cimagl(in[j]));
  }
  fclose(fh);
}

void set_weight() {
  unsigned long 	N = state.number_modes;
  long double		q, overN = 1.L/N;
  for (long int j = 0; j < N-1; j++) {
    q = PI*(2.0L*(j+1.L)*overN - 1.0L);
    M[j] = 0.5L/powl(cosl(0.5L*q), 2)/creall(arrayQ[j]*conjl(arrayQ[j]));
  }
}

void prepare_array(fftwl_complex *in) {
  unsigned long N = state.number_modes;
  for (unsigned long j = 0; j < N-1; j++) {
    W[j] = in[j+1] - in[0];
  }
  set_weight();
}

fftwl_complex dot(fftwl_complex *in1, fftwl_complex *in2) {
   unsigned long N = state.number_modes;
   fftwl_complex  rvalue = 0.L;
   for (unsigned long j = 0; j < N-1; j++) {
      rvalue += in1[j]*conjl(in2[j])*M[j];
   }
   return rvalue;
}

void allocate_pade(unsigned long nD) {
  unsigned long 	N = state.number_modes;
  long double 		q, overN = 1.L/N;
  fftwl_complex		xi;

  //pade_arr = fftwl_malloc(2*sizeof(pade)); 
  M = fftwl_malloc((N-1)*sizeof(long double));
  W = fftwl_malloc((N-1)*sizeof(fftwl_complex));
  tmp = fftwl_malloc((N-1)*sizeof(fftwl_complex));
  arrayQ = fftwl_malloc((N - 1)*sizeof(fftwl_complex));
  arrayP = fftwl_malloc((N - 1)*sizeof(fftwl_complex));
  //C = fftwl_malloc(4*sizeof(fftwl_complex *));
  //D = fftwl_malloc(4*sizeof(fftwl_complex *));
  for (int k = 0; k < 4; k++) {
    Gees[k] = fftwl_malloc((N - 1)*sizeof(fftwl_complex));
    Cees[k] = fftwl_malloc((2*nD + 1)*sizeof(fftwl_complex));
    A[k] = fftwl_malloc((N-1)*sizeof(fftwl_complex));
    B[k] = fftwl_malloc((N-1)*sizeof(fftwl_complex));
    //C[k] = fftwl_malloc((N-1)*sizeof(fftwl_complex));
    //D[k] = fftwl_malloc((N-1)*sizeof(fftwl_complex));
  }
  for (unsigned long j = 0; j < N-1; j++) arrayQ[j] = 1.0L;
  for (int k = 0; k < nD; k++) {
    xi = -1.0L + 1.0L*cexpl(0.5L*PI*(k + 0.5L)/nD); 
    for (unsigned long j = 0; j < N-1; j++) {
      q = PI*(2.0L*(j + 1.L)*overN - 1.0L);
      arrayQ[j] = arrayQ[j]*(tanl(0.5L*q) + 1.IL*xi);  // do renormaliztion
    }
  }
  prepare_array(data[0]);
  //pade_complex_out("Q0.txt", W);
  //pade_real_out("Weight.txt", M);
  //printf("Dot Prod: Re = %.12LE\tIm = %.12LE\n", creall(dot(W,W)), cimagl(dot(W,W)));
}

void deallocate_pade() {
   for (int k = 0; k < 4; k++) {
     //fftwl_free(D[k]);
     //fftwl_free(C[k]);
     fftwl_free(B[3-k]);
     fftwl_free(A[3-k]);
     fftwl_free(Cees[3-k]);
     fftwl_free(Gees[3-k]);
   }
   fftwl_free(arrayP); 
   fftwl_free(arrayQ);
   fftwl_free(tmp);
   fftwl_free(W);
   fftwl_free(M);
   //exit(1);
}

void evaluate_poly_array(unsigned long nD) {
  unsigned long N = state.number_modes;
  long double	s, overN = 1.L/N;
  
  for (unsigned long j = 0; j < N-1; j++) {
    s = tanl(0.5L*PI*(2.L*(j+1.L)*overN - 1.L));
    // --------   step 0
    A[0][j] = 0.L;
    B[0][j] = 1.L;
    // --------   step 1
    A[1][j] = -1.L;
    B[1][j] = -Cees[0][1];
    // --------   step 2
    A[2][j] = s*A[0][j] - Cees[0][2]*A[1][j] - Cees[1][2]*A[0][j];
    B[2][j] = s*B[0][j] - Cees[0][2]*B[1][j] - Cees[1][2]*B[0][j];
    // --------   step 3
    A[3][j] = s*A[1][j] - Cees[0][3]*A[2][j] - Cees[1][3]*A[1][j] - Cees[2][3]*A[0][j];
    B[3][j] = s*B[1][j] - Cees[0][3]*B[2][j] - Cees[1][3]*B[1][j] - Cees[2][3]*B[0][j];
    // --------   step 4+
    for (unsigned long k = 4; k < 2*nD + 1; k++) {
      tmp[j] = s*A[2][j] - Cees[0][k]*A[3][j] - Cees[1][k]*A[2][j] - Cees[2][k]*A[1][j] - Cees[3][k]*A[0][j];	
      A[0][j] = A[1][j];
      A[1][j] = A[2][j];
      A[2][j] = A[3][j];
      A[3][j] = tmp[j]; 

      tmp[j] = s*B[2][j] - Cees[0][k]*B[3][j] - Cees[1][k]*B[2][j] - Cees[2][k]*B[1][j] - Cees[3][k]*B[0][j];
      B[0][j] = B[1][j];
      B[1][j] = B[2][j];
      B[2][j] = B[3][j];
      B[3][j] = tmp[j]; 
    }  
  }
  memcpy(arrayQ, B[3], (N-1)*sizeof(fftwl_complex));
  memcpy(arrayP, A[3], (N-1)*sizeof(fftwl_complex));
}

void fmul_sub(fftwl_complex *in1, fftwl_complex *in2, fftwl_complex cnum, fftwl_complex *out) {
  unsigned long N = state.number_modes;
  for (unsigned long j = 0; j < N-1; j++) {
    out[j] = in1[j] - cnum*in2[j];
  }
}

void sigma_mul(fftwl_complex *in, fftwl_complex *out) {
  unsigned long N = state.number_modes;
  long double	s, overN = 1.L/N;
  for (unsigned long j = 0; j < N-1; j++) {
    s = tanl(0.5L*PI*(2.L*(j+1.L)*overN - 1.L));
    out[j] = s*in[j];
  }
}

void gram_schmidt(unsigned long nD){
  unsigned long N = state.number_modes;
  memset(Ens, 0, 4*sizeof(long double));
  for (unsigned k = 0; k < 4; k++) {
    memset(Cees[k], 0, (2*nD + 1)*sizeof(fftwl_complex));
  }
  // step 0
  memcpy(Gees[0], W, (N-1)*sizeof(fftwl_complex));
  Ens[0] = 1.0L/dot(Gees[0],Gees[0]);
  // step 1
  for (unsigned long j = 0; j < N-1; j++) tmp[j] = 1.L;
  Cees[0][1] = dot(tmp, Gees[0])*Ens[0];
  fmul_sub(tmp, Gees[0], Cees[0][1], Gees[1]);
  Ens[1] = 1.0L/dot(Gees[1],Gees[1]);
  // step 2
  sigma_mul(Gees[0], tmp);
  Cees[0][2] = dot(tmp, Gees[1])*Ens[1];
  fmul_sub(tmp, Gees[1], Cees[0][2], tmp);
  Cees[1][2] = dot(tmp, Gees[0])*Ens[0];
  fmul_sub(tmp, Gees[0], Cees[1][2], Gees[2]);
  Ens[2] = 1.0L/dot(Gees[2],Gees[2]);
  // step 3
  sigma_mul(Gees[1], tmp);
  Cees[0][3] = dot(tmp, Gees[2])*Ens[2];
  fmul_sub(tmp, Gees[2], Cees[0][3], tmp);
  Cees[1][3] = dot(tmp, Gees[1])*Ens[1];
  fmul_sub(tmp, Gees[1], Cees[1][3], tmp);
  Cees[2][3] = dot(tmp, Gees[0])*Ens[0];
  fmul_sub(tmp, Gees[0], Cees[2][3], Gees[3]);
  Ens[3] = 1.0L/dot(Gees[3],Gees[3]);
  // step 4+
  for (long int j = 4; j < 2*nD + 1; j++) {
      sigma_mul(Gees[2], tmp);
      Cees[0][j] = dot(tmp, Gees[3])*Ens[3];
      fmul_sub(tmp, Gees[3], Cees[0][j], tmp);
      Cees[1][j] = dot(tmp, Gees[2])*Ens[2];
      fmul_sub(tmp, Gees[2], Cees[1][j], tmp);
      Cees[2][j] = dot(tmp, Gees[1])*Ens[1];
      fmul_sub(tmp, Gees[1], Cees[2][j], tmp);
      Cees[3][j] = dot(tmp, Gees[0])*Ens[0];
      fmul_sub(tmp, Gees[0], Cees[3][j], tmp);
   
      memcpy(Gees[0], Gees[1], (N-1)*sizeof(fftwl_complex));
      memcpy(Gees[1], Gees[2], (N-1)*sizeof(fftwl_complex));
      memcpy(Gees[2], Gees[3], (N-1)*sizeof(fftwl_complex));
      memcpy(Gees[3], tmp, (N-1)*sizeof(fftwl_complex));

      Ens[0] = Ens[1];
      Ens[1] = Ens[2];
      Ens[2] = Ens[3];	
      Ens[3] = 1.0L/dot(tmp, tmp);
  }
}

void compute_rational(unsigned long nD, unsigned long n_max_iter) {
  pade_data.n_lins = 0;
  pade_data.n_poles = nD;

  allocate_pade(nD);
  gram_schmidt(nD);
  evaluate_poly_array(nD);
  find_l2_error(&pade_data);  
  pade_data.n_lins++;

  //print_pade(&pade_data);

  for (unsigned int j = 1; j < n_max_iter; j++) {
    set_weight();
    gram_schmidt(nD);
    evaluate_poly_array(nD);
    find_l2_error(&pade_data);  
    pade_data.n_lins++;
    //print_pade(&pade_data);
  }
  set_weight();
  //for (unsigned int j = 0; j < N-1; j++){
    //tmp[j] = arrayP[j]/arrayQ[j] - W[j];
  //}
  //pade_real_out("m5.txt", M);
  //pade_complex_out("rational.txt", tmp);
  //deallocate_pade();
  //exit(1);
}

void optimal_pade() {
  FILE *fh = fopen("pc_rate.txt","w");
  unsigned long nd = 1, l_iters = 8;
  pade best_pade;  
  best_pade.l2_rel_err = 1.L;
  
  fprintf(fh, "# 1. nD, number poles 2. Error\n\n");
  //for (unsigned int nd = 1; nd < 32; nd++) {
  while (nd < 32) {
    nd++;
    best_pade.n_poles = nd;
    pade_data.n_lins = 0;
    compute_rational(nd, l_iters);
    deallocate_pade();
    if (pade_data.l2_rel_err < best_pade.l2_rel_err) {
      best_pade.l2_rel_err = pade_data.l2_rel_err;
      best_pade.l2_abs_err = pade_data.l2_abs_err;
      best_pade.l2_nrm = pade_data.l2_nrm;
      best_pade.n_poles = nd;
      best_pade.n_lins = l_iters;
      printf("Relative Error (nd = %3lu) = %11.5LE\n", nd, best_pade.l2_rel_err);
      fprintf(fh, "%.3lu\t%11.5LE\n", nd, best_pade.l2_rel_err);
    } else {
      nd = nd-2;
      best_pade.n_poles = nd;
      pade_data.n_lins = 0;
      compute_rational(nd, l_iters);
      newton_search(nd);
      deallocate_pade();
      best_pade.l2_rel_err = pade_data.l2_rel_err;
      best_pade.l2_abs_err = pade_data.l2_abs_err;
      best_pade.l2_nrm = pade_data.l2_nrm;
      best_pade.n_poles = nd;
      best_pade.n_lins = l_iters;
      printf("Safe Relative Error (nd = %3lu) = %11.5LE\n", nd, best_pade.l2_rel_err);
      break;
    } 
  }
  fclose(fh);
  exit(1);
}

void print_pade(pade_ptr inp) { 
  printf("Iter #%3u:\t", inp->n_lins); 
  printf("%.18LE\n", inp->l2_rel_err); 
}

void find_l2_error(pade_ptr inp){
  unsigned long N = state.number_modes;
  fftwl_complex	tmp;
  long double overN = 2.L*PI/N;
  inp->l2_nrm = 0.L;
  inp->l2_abs_err = 0.L;

  for (unsigned int j = 0; j < N-1; j++) {
    tmp = (arrayP[j]/arrayQ[j] - W[j]);
    inp->l2_nrm += creall(W[j]*conjl(W[j]));
    inp->l2_abs_err += creall(tmp*conjl(tmp));    
  }
  inp->l2_nrm = sqrtl(inp->l2_nrm)*overN;
  inp->l2_abs_err = sqrtl(inp->l2_abs_err)*overN;
  inp->l2_rel_err = (inp->l2_abs_err)/(inp->l2_nrm);
}


void evaluate_poly(fftwl_complex *in, unsigned long nD, fftwl_complex *outQ, fftwl_complex *outQp, fftwl_complex *outP){
  fftwl_complex  	tmp_C, tmp_B, tmp_A;
  fftwl_complex		P[4], Q[4], dQ[4];
  fftwl_complex		s;
   
  s = ctanl(0.5L*(*in));
  P[0] = 0.L;
  dQ[0] = 0.L;
  Q[0] = 1.L;
  P[1] = -1.L;
  dQ[1] = 0.L;
  Q[1] = -Cees[0][1];

  P[2]  = s*P[0] - Cees[0][2]*P[1] - Cees[1][2]*P[0];
  dQ[2] = Q[0] + s*dQ[0] - Cees[0][2]*dQ[1] - Cees[1][2]*dQ[0];
  Q[2]  = s*Q[0] - Cees[0][2]*Q[1] - Cees[1][2]*Q[0];

  P[3]  = s*P[1] - Cees[0][3]*P[2] - Cees[1][3]*P[1] - Cees[2][3]*P[0];
  dQ[3] = Q[1] + s*dQ[1] - Cees[0][3]*dQ[2] - Cees[1][3]*dQ[1] - Cees[2][3]*dQ[0];
  Q[3]  = s*Q[1] - Cees[0][3]*Q[2] - Cees[1][3]*Q[1] - Cees[2][3]*Q[0];
  for (unsigned long k = 4; k < 2*nD + 1; k++) {
    tmp_C = s*P[2] - Cees[0][k]*P[3] - Cees[1][k]*P[2] - Cees[2][k]*P[1] - Cees[3][k]*P[0];
    tmp_B = s*Q[2] - Cees[0][k]*Q[3] - Cees[1][k]*Q[2] - Cees[2][k]*Q[1] - Cees[3][k]*Q[0];
    tmp_A = Q[2] + s*dQ[2] - Cees[0][k]*dQ[3] - Cees[1][k]*dQ[2] - Cees[2][k]*dQ[1] - Cees[3][k]*dQ[0];
    P[0] = P[1];
    P[1] = P[2];
    P[2] = P[3];
    P[3] = tmp_C;
    Q[0] = Q[1];
    Q[1] = Q[2];
    Q[2] = Q[3];
    Q[3] = tmp_B;
    dQ[0] = dQ[1];
    dQ[1] = dQ[2];
    dQ[2] = dQ[3];
    dQ[3] = tmp_A; 
  } 
  *outP  = P[3];  
  *outQ  = Q[3];
  *outQp = dQ[3];
}


void newton_search(unsigned long nD) {
  fftwl_complex *roots = fftwl_malloc(nD*sizeof(fftwl_complex));
  fftwl_complex *residues = fftwl_malloc(nD*sizeof(fftwl_complex));
  fftwl_complex z = 1.0L+0.5IL;
  fftwl_complex p = 0.L;
  fftwl_complex f = 1.L, df;
  fftwl_complex inv_sum;
  long int counter = 0;
  char str[80];
  FILE *fh;

  for (int j = 0; j < nD; j++) {  
    sprintf(str, "./roots/root_%03d.txt", j);
    fh = fopen(str, "w");
    z = 1.0L+0.5IL;
    f = 1.L; 
    inv_sum = 1.L;
    while ( cabsl(f*inv_sum) > 5.0E-29L) {
      evaluate_poly(&z, nD, &f, &df, &p);
      if ( j == 0) {
	inv_sum = 1.L;
        z = z - f/df;
      } else {     
        inv_sum = 0.L;
        for (int k = 0; k < j; k++) {
          inv_sum += 2.L/(z - roots[k]);
        }
        z = z - f/(df - f*inv_sum);
      }
      counter++;
      fprintf(fh, "%4ld\t%26.18LE\t%26.18LE\t%26.18LE\t%26.18LE\t%26.18LE\n", counter, creall(z), cimagl(z), creall(p/df), cimagl(p/df), cabsl(f*inv_sum));
      if (counter == 200) break;
    }
    printf("Root %3d: Iteration %4ld |f| = %14.8LE\n", j, counter, cabsl(f*inv_sum));
    counter = 0;
    fclose(fh);
    residues[j] = p/df;
    roots[j] = z;
  }
  verify_pade(residues, roots, nD);
  fh = fopen("./roots/roots.txt", "w");
  fprintf(fh, "# 1. root # 2.-3. z_k 4.-5. g_k\n\n");
  for (int j = 0; j < nD; j++) {
    fprintf(fh, "%4d\t%26.18LE\t%26.18LE\t%26.18LE\t%26.18LE\n",j, creall(roots[j]), cimagl(roots[j]), creall(residues[j]), cimagl(residues[j]));
  }
  fclose(fh);
  fftwl_free(residues);
  fftwl_free(roots);
}

void verify_pade(fftwl_complex *residues, fftwl_complex *roots, unsigned int nD) {
   unsigned long 	N = state.number_modes;
   long double 		overN = 1.L/N, s;
   long double		tmp = 0.L;
   
   fftwl_complex *pade = fftwl_malloc((N-1)*sizeof(fftwl_complex));
   memset(pade, 0, (N-1)*sizeof(fftwl_complex));
   for (int k = 0; k < N-1; k++ ) {
     s = PI*(2.L*(k+1)*overN - 1.L);
     s = tanl(0.5L*s);
     for (int j = 0; j < nD; j++) {
       pade[k] += residues[j]/(s - roots[j]);
     }
     tmp += powl(cabsl(W[k]-pade[k]),2); 
   }
   tmp = sqrtl(tmp)*overN;
   pade_complex_out("target.txt", W);
   pade_complex_out("pade.txt", pade);
   fftwl_free(pade);
}
