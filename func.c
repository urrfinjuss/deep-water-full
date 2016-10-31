#include "header.h"


//static rk4_data rk4;
static params input, rs_input;
static data array, bs_tmp, HS_bs_tmp;
static aux extra;
static consq motion;
static fftwl_complex *tmp0, *tmp1, *tmp2, *tmp3, *qB, *vB;
static fftwl_complex q, b0, b1;
static fftwl_plan p0f, p1f, p2f;
static fftwl_plan p0b, p1b;
static long double *k, A0, A1;
static long double D0, D1;
static long int ref_nog;

void debug_msg(char* in, int EXITF) {
  FILE *fh = fopen("output.log","a");
  fprintf(fh, "%s", in);
  fclose(fh);
  if (EXITF) exit(1);
}

//-------------------

/*void init_rk4() { 
   //rk4.dt = 0.25L*PI*input.L/input.N;
   //rk4.dt = 45.L*PI*powl(input.L/input.N, 1.5L); // too weak, instability
   //rk4.dt = 32768.L*PI*powl(input.L/input.N,2);  // adjusted to 32768 from 8096  for VZ (needs corrections)
   //rk4.dt = 0.25L*PI*powl(input.L/input.N,1);      // overturning Pavel (no st)
   if ( input.s ) {
     rk4.dt = fminl(0.25L*PI*powl(input.L/input.N,1), 0.125L*powl(input.L/input.N, 1.5L))/sqrtl(input.s);
     //printf("Using Capillary CFL\n");
   } else {
     //rk4.dt = 0.25L*PI*powl(input.L/input.N,1);      // overturning Pavel (no st) run1
     //rk4.dt = 0.125L*fminf(powl(input.L,1),1.0L)*PI*powl(1024.0L/input.N, 1.0L)/1024.L;      // overturning Pavel (no st) run4
     rk4.dt = 1.L*PI*fminf(1./input.L, input.L)/input.N;
   }
   //rk4.dt = 0.25L*PI*powl(input.L/input.N,1.5L);      // overturning Pavel (with st)
   rk4.tshift += rk4.time;
   rk4.time = 0.L;
   rk4.nsteps = 100000;
   rk4.D = 0.25L*powl(95.L*input.N/256.L,-12)*rk4.dt;  // adjusted to 0.25L from 1.0L
   
   if (!(rk4.tmpq = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(rk4.k0q = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(rk4.k1q = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(rk4.k2q = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(rk4.k3q = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);

   if (!(rk4.tmpv = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(rk4.k0v = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(rk4.k1v = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(rk4.k2v = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(rk4.k3v = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
}

void free_rk4() { 
  fftwl_free(rk4.tmpq);	  fftwl_free(rk4.tmpv);
  fftwl_free(rk4.k0q);	  fftwl_free(rk4.k0v);
  fftwl_free(rk4.k1q);	  fftwl_free(rk4.k1v);
  fftwl_free(rk4.k2q);	  fftwl_free(rk4.k2v);
  fftwl_free(rk4.k3q);	  fftwl_free(rk4.k3v);
}*/


void compute_aux_arrays(fftwl_complex *inQ, fftwl_complex *inV) {
  //  --- verified 
  memcpy(tmp0, inQ, (input.N)*sizeof(fftwl_complex));
  memcpy(tmp1, inV, (input.N)*sizeof(fftwl_complex));
  fftwl_execute(p0f);
  fftwl_execute(p1f);
  for (int j = 0; j < input.N; j++) {
    tmp0[j] = 1.IL*k[j]*tmp0[j]/input.N;
    tmp1[j] = 1.IL*k[j]*tmp1[j]/input.N;
  }
  fftwl_execute(p0b);
  fftwl_execute(p1b);
  memcpy(extra.dQ, tmp0,(input.N)*sizeof(fftwl_complex));
  memcpy(extra.dV, tmp1,(input.N)*sizeof(fftwl_complex));    

  for (int j = 0; j < input.N; j++) {
    tmp1[j] = 2.L*creall(inV[j]*conjl(inQ[j]*inQ[j]));
    tmp0[j] = inV[j]*conjl(inV[j]) + 4.L*input.s*cimagl(tmp0[j]*conjl(inQ[j]));
  }
  fftwl_execute(p1b);  
  fftwl_execute(p0b);
  tmp1[input.N/2] = 0.L; extra.b0 = 0.L;
  for (long int j = input.N/2 - 2; j > -1; j--) {
    extra.b0 = extra.b0 + 0.5L*(tmp1[input.N-j-1]*conjl(extra.w[j]) - tmp1[j+1]*extra.w[j] );
  }
  memset(tmp2,0,sizeof(fftwl_complex)*input.N);
  memset(tmp1+input.N/2,0,sizeof(fftwl_complex)*input.N/2);
  memset(tmp0+input.N/2,0,sizeof(fftwl_complex)*input.N/2);
  tmp1[0] = tmp1[0]/2;
  for (int j = 0; j < (input.N)/2; j++) {
    tmp0[j] = -1.IL*k[j]*tmp0[j]/(input.N);
    tmp1[j] = tmp1[j]/(input.N);
    tmp2[j] = -1.IL*k[j]*tmp1[j];
  }
  fftwl_execute(p2f);  // tmp2 contains dU
  fftwl_execute(p1f);  // tmp1 contains  U
  fftwl_execute(p0f);  // tmp0 contains  B
  for (int j = 0; j < input.N; j++) tmp1[j] = tmp1[j] + extra.b0;
  memcpy( extra.B, tmp0, input.N*sizeof(fftwl_complex));
  memcpy( extra.U, tmp1, input.N*sizeof(fftwl_complex));
  memcpy(extra.dU, tmp2, input.N*sizeof(fftwl_complex));
}

/*void rk4_step() {
  memcpy(rk4.tmpq, array.Q, input.N*sizeof(fftwl_complex));
  memcpy(rk4.tmpv, array.V, input.N*sizeof(fftwl_complex));

  compute_aux_arrays(rk4.tmpq, rk4.tmpv);
  for (long int j = 0; j < input.N; j++) {
    rk4.k0q[j] = 0.5IL*(2.0L*extra.dQ[j]*extra.U[j] - extra.dU[j]*rk4.tmpq[j])/extra.dq[j];
    rk4.k0v[j] = 1.0IL*(extra.U[j]*extra.dV[j] - rk4.tmpq[j]*rk4.tmpq[j]*extra.B[j])/extra.dq[j] + input.g*(rk4.tmpq[j]*rk4.tmpq[j] - 1.0L);
    rk4.tmpq[j] = array.Q[j] + 0.5L*rk4.dt*rk4.k0q[j];
    rk4.tmpv[j] = array.V[j] + 0.5L*rk4.dt*rk4.k0v[j];
  }

  compute_aux_arrays(rk4.tmpq, rk4.tmpv);
  for (long int j = 0; j < input.N; j++) {
    rk4.k1q[j] = 0.5IL*(2.0L*extra.dQ[j]*extra.U[j] - extra.dU[j]*rk4.tmpq[j])/extra.dq[j];
    rk4.k1v[j] = 1.0IL*(extra.U[j]*extra.dV[j] - rk4.tmpq[j]*rk4.tmpq[j]*extra.B[j])/extra.dq[j] + input.g*(rk4.tmpq[j]*rk4.tmpq[j] - 1.0L);
    rk4.tmpq[j] = array.Q[j] + 0.5L*rk4.dt*rk4.k1q[j];
    rk4.tmpv[j] = array.V[j] + 0.5L*rk4.dt*rk4.k1v[j];
  }

  compute_aux_arrays(rk4.tmpq, rk4.tmpv);
  for (long int j = 0; j < input.N; j++) {
    rk4.k2q[j] = 0.5IL*(2.0L*extra.dQ[j]*extra.U[j] - extra.dU[j]*rk4.tmpq[j])/extra.dq[j];
    rk4.k2v[j] = 1.0IL*(extra.U[j]*extra.dV[j] - rk4.tmpq[j]*rk4.tmpq[j]*extra.B[j])/extra.dq[j] + input.g*(rk4.tmpq[j]*rk4.tmpq[j] - 1.0L);
    rk4.tmpq[j] = array.Q[j] + 1.0L*rk4.dt*rk4.k2q[j];
    rk4.tmpv[j] = array.V[j] + 1.0L*rk4.dt*rk4.k2v[j];
  }

  compute_aux_arrays(rk4.tmpq, rk4.tmpv);
  for (long int j = 0; j < input.N; j++) {
    rk4.k3q[j] = 0.5IL*(2.0L*extra.dQ[j]*extra.U[j] - extra.dU[j]*rk4.tmpq[j])/extra.dq[j];
    rk4.k3v[j] = 1.0IL*(extra.U[j]*extra.dV[j] - rk4.tmpq[j]*rk4.tmpq[j]*extra.B[j])/extra.dq[j] + input.g*(rk4.tmpq[j]*rk4.tmpq[j] - 1.0L);
    array.Q[j] += 1.0L*rk4.dt*(rk4.k0q[j] + 2.0L*rk4.k1q[j] + 2.0L*rk4.k2q[j] + rk4.k3q[j])/6.0L;
    array.V[j] += 1.0L*rk4.dt*(rk4.k0v[j] + 2.0L*rk4.k1v[j] + 2.0L*rk4.k2v[j] + rk4.k3v[j])/6.0L;
    tmp0[j] = array.Q[j]/input.N;
    tmp1[j] = array.V[j]/input.N;
  }
  fftwl_execute(p0b);
  fftwl_execute(p1b);
  memset(tmp0+input.N/2,0,sizeof(fftwl_complex)*input.N/2);
  memset(tmp1+input.N/2,0,sizeof(fftwl_complex)*input.N/2);
  for (long int j = 0; j < input.N/2; j++) {
    tmp0[j] = tmp0[j]*expl(-rk4.D*powl(j,12));
    tmp1[j] = tmp1[j]*expl(-rk4.D*powl(j,12)); 
  }
  fftwl_execute(p0f);
  fftwl_execute(p1f);
  memcpy(array.Q, tmp0, input.N*sizeof(fftwl_complex));
  memcpy(array.V, tmp1, input.N*sizeof(fftwl_complex));

}*/


//-----  Output 

void basic_output(char *fname){
  FILE *fh = fopen(fname,"w");
  long double u;

  fprintf(fh, "# 1. u 2.-3. Q 4.-5. V\n# time = %.12LE\n\n", rk4.time+rk4.tshift);
  for (int j = 0; j < input.N; j++) {
     q = 2.L*PI*(1.L*j/input.N - 0.5L);
     u = input.u + 2.L*atan2l(input.L*sinl(0.5L*(q-extra.q)), cosl(0.5L*(q-extra.q)));
     fprintf(fh, "%.19LE\t%.19LE\t%.19LE\t%.19LE\t%.19LE\n", u, creall(array.Q[j]), cimagl(array.Q[j]), creall(array.V[j]), cimagl(array.V[j]) );
   }
  fclose(fh);
}

void spec_output(char *fname){
#if SPEC_OUTPUT
  memcpy(tmp0, array.Q, input.N*sizeof(fftwl_complex));
  memcpy(tmp1, array.V, input.N*sizeof(fftwl_complex));
  fftwl_execute(p1b);
  fftwl_execute(p0b);
  FILE *fh = fopen(fname,"w");
  fprintf(fh, "# 1. k 2. |Q_k| 3. |V_k|\n# time = %.12LE\n\n", rk4.time+rk4.tshift);
  for (long int j = 0; j < input.N/2; j++) {
     fprintf(fh, "%ld\t%.19LE\t%.19LE\n", -j, cabsl(tmp0[j])/input.N, cabsl(tmp1[j])/input.N);
  }
  fclose(fh);
#endif
}

void spec_output_nofft(char *fname){
#if SPEC_OUTPUT
  FILE *fh = fopen(fname,"w");
  fprintf(fh, "# 1. k 2. |Q_k| 3. |V_k|\n# time = %.12LE\n\n", rk4.time+rk4.tshift);
  for (long int j = 0; j < input.N/2; j++) {
     fprintf(fh, "%ld\t%.19LE\t%.19LE\n", -j, cabsl(tmp0[j])/input.N, cabsl(tmp1[j])/input.N);
  }
  fclose(fh);
#endif
}

void get_complex_advection(char *fname) {
  long double q, u;
  FILE *fh = fopen(fname,"w");

  for (long int j = 0; j < input.N; j++) {
    tmp0[j] = array.Q[j]*array.Q[j]*conjl(array.V[j]);
    tmp0[j] += conjl(tmp0[j]);
    tmp0[j] = tmp0[j]/(input.N);
  }
  fftwl_execute(p0b);
  memset(tmp0+input.N/2,0,sizeof(fftwl_complex)*input.N/2);
  tmp0[0] = 0.L;
  fftwl_execute(p0f);
  fprintf(fh, "# 1. k 2. re U 3. im U\n# time = %.12LE\n\n", rk4.time+rk4.tshift);
  for (long int j = 0; j < input.N; j++) {
     q = 2.L*PI*(1.L*j/input.N - 0.5L);
     u = input.u + 2.L*atan2l(input.L*sinl(0.5L*(q-extra.q)), cosl(0.5L*(q-extra.q)));
     fprintf(fh, "%.12LE\t%.19LE\t%.19LE\n", u, creall(tmp0[j]), cimagl(tmp0[j]));
   }
  fclose(fh);
  
}

void get_integrals(char *fname) {
  long double u;

  for (long int j = 0; j < input.N; j++) {
    tmp0[j] = array.V[j]*extra.dq[j]/cpowl(array.Q[j], 2.0L)/input.N;
    tmp1[j] = -1.0IL*extra.dq[j]*(1.0L - cpowl(array.Q[j], -2.0L))/input.N;
  }
  fftwl_execute(p0b);
  fftwl_execute(p1b);  
  motion.K = 0.L;
  tmp1[0] = 0.L;
  memset(tmp1+input.N/2,0,sizeof(fftwl_complex)*input.N/2);
  for (long int j = input.N/2-1; j > 0; j--) {
	motion.K += 0.5L*creall(tmp0[j]*conjl(tmp0[j])/k[j]);
        tmp1[j] = tmp1[j]/k[j];
  }
  fftwl_execute(p1f);  
  for (long int j = 0; j < input.N; j++) {
	tmp1[j] = tmp1[j];  			// zt
	tmp0[j] = cimagl(tmp1[j])/input.N;	
  }
  fftwl_execute(p0f);  
  for (long int j = 1; j < input.N; j++) {
	tmp0[j] = fabsl(k[j])*tmp0[j];   // |k|y
  }
  fftwl_execute(p0b);
  motion.ml = 0.L;
  motion.y0 = 0.L;
  motion.P = 0.L;
  for (long int j = 0; j < input.N; j++) motion.y0 += -cimagl(tmp1[j])*(extra.dq[j] + tmp0[j])/input.N;	

  FILE *fh = fopen(fname,"w");
  fprintf(fh, "# 1. u 2. x 3. y\n# N = %ld\tL = %.19LE\tT = %.19LE\n\n", input.N, input.L, rk4.tshift+rk4.time);
  for (long int j = 0; j < input.N; j++) {
     q = 2.L*PI*(1.L*j/input.N - 0.5L);
     u = input.u + 2.L*atan2l(input.L*sinl(0.5L*(q-extra.q)), cosl(0.5L*(q-extra.q)));
     tmp2[j] = u + (tmp1[j]-creall(tmp1[0])) + 1.0IL*motion.y0;
     fprintf(fh, "%.19LE\t%.19LE\t%.19LE\n", u, creall(tmp2[j]), cimagl(tmp2[j]) );
     motion.ml += (extra.dq[j] + tmp0[j])*cimagl(tmp2[j]);
     motion.P += input.g*((extra.dq[j] + tmp0[j])*cimagl(tmp2[j])*cimagl(tmp2[j]));
  }
  fclose(fh);
  motion.ml = motion.ml/input.N;
  motion.P = motion.P/input.N;
  motion.H = cimagl(tmp2[input.N/2]-tmp2[0]);
  motion.C = 2.L*cimagl(tmp2[input.N/2]-tmp2[input.N/2-1])/powl(creall(tmp2[input.N/2] - tmp2[input.N/2-1]),2);

}

long double track_singularity(){
  long double Min = 1.L, tmpl = 1.L, Umin = 0.L, u;

  for (long int j = 0; j < input.N; j++) {
    tmpl = creall(array.Q[j]*conjl(array.Q[j]));
    q = 2.L*PI*(1.L*j/input.N - 0.5L);
    u = input.u + 2.L*atan2l(input.L*sinl(0.5L*(q-extra.q)), cosl(0.5L*(q-extra.q)));
    if (tmpl < Min) {
	Min = tmpl;
        Umin = u;
    }
  }
#if MOVE_MESH
  return Umin;
#else
  return 0.L;
#endif
}

fftwl_complex check_resolution() {
  long double EQ0 = 0.L, EQ1 = 0.L, EQ2 = 0.L, EQ4 = 0.L, EQ8 = 0.L;
  long double EV0 = 0.L, EV1 = 0.L, EV2 = 0.L, EV4 = 0.L, EV8 = 0.L; 
  long double rq0, rq1, rq2, rq3;   
  long double rv0, rv1, rv2, rv3;


  memcpy(tmp0, array.Q, input.N*sizeof(fftwl_complex));
  memcpy(tmp1, array.V, input.N*sizeof(fftwl_complex));
  fftwl_execute(p0b);
  fftwl_execute(p1b);
  //spec_output_nofft("amispectrum.txt");
  for (long int j = input.N/2; j > -1; j--) {
    EQ0 += creall(tmp0[j]*conjl(tmp0[j]));
    EV0 += creall(tmp1[j]*conjl(tmp1[j]));
    if (j >= input.N/8) {
       EQ1 += creall(tmp0[j]*conjl(tmp0[j]));
       EV1 += creall(tmp1[j]*conjl(tmp1[j]));
       if (j >= 3*input.N/16) {
          EQ2 += creall(tmp0[j]*conjl(tmp0[j]));
          EV2 += creall(tmp1[j]*conjl(tmp1[j]));
          if (j >= 7*input.N/32) {
   	    EQ4 += creall(tmp0[j]*conjl(tmp0[j]));
            EV4 += creall(tmp1[j]*conjl(tmp1[j]));
	    if (j >= 15*input.N/64) {
	       EQ8 += creall(tmp0[j]*conjl(tmp0[j]));
               EV8 += creall(tmp1[j]*conjl(tmp1[j]));
	    }
          }
       }
    }
  }
  rq0 = sqrtl(EQ1/EQ0);
  rq1 = sqrtl(EQ2/EQ0);
  rq2 = sqrtl(EQ4/EQ0);
  rq3 = sqrtl(EQ8/EQ0);

  rv0 = sqrtl(EV1/EV0);
  rv1 = sqrtl(EV2/EV0);
  rv2 = sqrtl(EV4/EV0);
  rv3 = sqrtl(EV8/EV0);
  
  return rq3+1.IL*rv3;
}
void backup_arrays() {
  qB = fftwl_malloc(input.N*sizeof(fftwl_complex));
  vB = fftwl_malloc(input.N*sizeof(fftwl_complex));
  // save arrays before L1 ref
  memcpy(qB, array.Q, input.N*sizeof(fftwl_complex)); 
  memcpy(vB, array.V, input.N*sizeof(fftwl_complex));
}

void restore_arrays() {
  memcpy(array.Q, qB, input.N*sizeof(fftwl_complex)); 
  memcpy(array.V, vB, input.N*sizeof(fftwl_complex));
}

/*long double adaptive_timestep() {
   memcpy(tmp0, array.Q, input.N*sizeof(fftwl_complex));
   memcpy(tmp1, array.V, input.N*sizeof(fftwl_complex));
   fftwl_execute(p0b);
   fftwl_execute(p1b);
   A0 = 0.L; A1 = 0.L;
   D0 = 0.L; D1 = 0.L;
   for (long int j = 0; j < 16; j++) {
     A0 += cabsl(tmp0[input.N/2-j]);
     A1 += cabsl(tmp1[input.N/2-j]);
     D0 += cabsl(tmp0[j]);
     D1 += cabsl(tmp1[j]);
   }
   return cabsl(A0/D0 + 1.IL*A1/D1);
}*/

void resolution_monitor_alt(long int *iter) {  
  // this actually works quite well
  // decrease twice by square root 2, then double # of points
  char fname[256];
  fftwl_complex mon = check_resolution();

  if (( creall(mon)>1e-12*sqrt(input.N))||(cimagl(mon)>1e-12*sqrt(input.N))) {
     sprintf(fname, "spec_br%ld.txt", input.refN+1);
     spec_output_nofft(fname);
     //get_integrals(fname);
     //sprintf(fname, "spec_zt%ld.txt", input.refN+1);
     //spec_output2(fname);

     printf("Before Ref: ML = %.12LE\ty0 = %.12LE\n", motion.ml, motion.y0);
     move_mesh(track_singularity(), input.L, input.N);
     mon = check_resolution();
     *iter = 0;
     if (( creall(mon)>2e-13*sqrt(input.N))||(cimagl(mon)>2e-13*sqrt(input.N))) {
        //if (ref_nog < 3) {  // for VZ 2; for Pavel 4
	     printf("L1, ref #%ld\n", ref_nog);
	     backup_arrays();
	     long double posS = track_singularity();
	     if (input.g) {
	       move_mesh(posS, input.L*sqrtl(2.0L), input.N);  // Lin = L/2 in Matlab
     	       mon = check_resolution();
	       printf("L1 with g: revert. ref #%ld\n", ref_nog);
	       if (( creall(mon)>2e-13*sqrt(input.N))||(cimagl(mon)>2e-13*sqrt(input.N))) {
	         printf("L1 with g: revert failed. ref #%ld\n", ref_nog);
      	         move_mesh(posS, input.L/sqrtl(2.0L), input.N);     
		 restore_arrays();	 	 
   	         move_mesh(posS, input.L/sqrtl(2.0L), input.N);  // Lin = L/2 in Matlab
     	         mon = check_resolution();
   	         rk4.nskip = floor(rk4.nskip*sqrtl(2.L));
 	       } else {
	         rk4.nskip = floor(rk4.nskip/sqrtl(2.L));
		 printf("L1 with g: revert success!. ref #%ld\n", ref_nog);
	       }
             } else {
	       printf("L1 with g: forward ref #%ld\n", ref_nog);
	       move_mesh(posS, input.L/sqrtl(2.0L), input.N);  // Lin = L/2 in Matlab
     	       mon = check_resolution();
	       rk4.nskip = floor(rk4.nskip*sqrtl(2.L));
 	     }
	     sprintf(fname, "spec_l1r%ld.txt", input.refN+1);
	     spec_output_nofft(fname);
   	     //get_integrals(fname);
	     //rk4.nskip = floor(rk4.nskip*sqrtl(2.L)); // Pavel no st
	     ref_nog++;	
	     printf("L1: Using Gravity CFL: dt = %.12LE\n", rk4.dt);
        //} else {
	if (( creall(mon)>2e-13*sqrt(input.N))||(cimagl(mon)>2e-13*sqrt(input.N))) {
	     // revert back if L1 failed
	     printf("L1, ref #%ld failed.\n", ref_nog-1);
   	     move_mesh(posS, input.L*sqrtl(2.0L), input.N);
   	     restore_arrays();	 	 	     
	     move_mesh(track_singularity(), input.L, 2*input.N);
     	     mon = check_resolution();
	     sprintf(fname, "spec_l2r%ld.txt", input.refN+1);
             spec_output_nofft(fname);
             printf("L2: Using Gravity CFL: dt = %.12LE\n", rk4.dt);
	     // experimental (dt slow down)
	     rk4.dt = rk4.dt/8.L;	 
	     //rk4.nskip = floor(rk4.nskip*sqrtl(2.L));  // Pavel no st
	     ref_nog = 0;
        } else printf("L1 with g: success!. ref #%ld\n", ref_nog-1);
        if (( creall(mon)>1e-12*sqrt(input.N))||(cimagl(mon)>1e-12*sqrt(input.N))) {
	   printf("Spectrum too wide. CFL?\n"); 
           sprintf(fname, "spec_last.txt");
           spec_output(fname);
	   exit(1);
        }
        sprintf(fname, "more.txt");
        get_integrals(fname);
        //sprintf(fname, "zt.txt");
        //spec_output(fname);
        printf("After Ref: ML = %.12LE\ty0 = %.12LE\n", motion.ml, motion.y0);
	input.refN++;
     }

    
  }
}


//-------------- Runge-Kutta 4

void evolve_rk4(){
  FILE *fh = fopen("time_rk4.log","w");
  char fname[256];
  long int l = 0;
  long int j = 0;
  //long double R1, R0;

  init_rk4_module(&input, &motion, &array, &extra);

  //rk4.time = 0.L;
  //rk4.tshift = 0.L;
  input.refN = 0;
  sprintf(fname, "data/data_%04ld.txt", l);
  get_integrals(fname);
  //basic_output(fname);
  fprintf(fh, "# 1. t 2. K 3. P 4. H 5. C 6. y0 7. ML\n\n");
  fprintf(fh, "% .19LE\t% .19LE\t% .19LE\t% .19LE\t% .19LE\t% .19LE\t% .19LE\n", rk4.time, motion.K, motion.P, motion.H, motion.C, motion.y0, motion.ml);
  fclose(fh);

  //R0 = fmaxl(adaptive_timestep(),1.0E-17); R1 = R0;
  //printf("T = 0.\tR1 = %.12LE\n", R1);
  check_resolution();
  sprintf(fname, "spec_r%ld.txt", input.refN);
  spec_output_nofft(fname);

  //for (long int j = 0; j < rk4.nsteps; j++) {
  while (1) {
    rk4_step();
    rk4.time = (j+1)*rk4.dt; 
    if (j % rk4.nskip == 0) {
      l++;
      /*R1 = adaptive_timestep();
      if (R1/R0 > 20) {
        R0 = R1;
	rk4.dt = rk4.dt/2.L;
        rk4.tshift += rk4.time;
        rk4.time = 0.L;
	printf("Time-Step Halved.\n");
      }*/
      sprintf(fname, "data/data_%04ld.txt", l);
      get_integrals(fname);
      //sprintf(fname, "data/basic_%04ld.txt", l);
      //basic_output(fname);
      sprintf(fname, "data/spec_%04ld.txt", l);
      spec_output(fname);
      //sprintf(fname, "data/cadvect_%04ld.txt", l);
      //get_complex_advection(fname);
      printf("T = % .8LE\tML = % .12LE\tP = % .12LE\tK = % .12LE\tH = % .12LE\n", rk4.tshift+rk4.time, motion.ml, motion.P, motion.K, motion.H);

      fh = fopen("time_rk4.log","a");
      fprintf(fh, "% .19LE\t% .19LE\t% .19LE\t% .19LE\t% .19LE\t% .19LE\t% .19LE\n", rk4.tshift+rk4.time, motion.K, motion.P, motion.H, motion.C, motion.y0, motion.ml);
      fclose(fh);
      //basic_output(fname);
    }
    j++;
    resolution_monitor_alt(&j);  // switched to alternative

  }
  sprintf(fname, "data/data_%04ld.txt", l+1);
  get_integrals(fname);
  printf("T = % .8LE\tML = % .12LE\tP = % .12LE\tK = % .12LE\tH = % .12LE\n", rk4.tshift+rk4.time, motion.ml, motion.P, motion.K, motion.H);
  fh = fopen("time_rk4.log","a");
  fprintf(fh, "% .19LE\t% .19LE\t% .19LE\t% .19LE\t% .19LE\n", rk4.tshift+rk4.time, motion.K, motion.P, motion.H, motion.C);
  fclose(fh);

}


void set_aux() {
  for (int j = 0; j < input.N; j++) {
    k[j] = 1.L*j;
    if (j > (input.N)/2-1) {
      k[j] = -1.L*(input.N - j);
    } 
  }
  extra.q = 2.0L*atan2l(input.L*sinl(0.5L*input.u),cosl(0.5L*input.u));
  for (int j = 0; j < (input.N)/2-1; j++) extra.w[j] = cexpl(-1.0IL*k[j+1]*(extra.q - 2.0IL*atanhl(input.L) ) )/input.N;	   //good
  for (int j = 0; j < input.N; j++) {
    q = 2.L*PI*(1.0L*j/input.N - 0.5L);
    extra.dq[j] = (2.0L*input.L)/(1.0L + input.L*input.L + (1.0L - input.L*input.L)*cosl(q-extra.q));
  }

}
void clean_up_mesh() {
   fftwl_free(array.Q);
   fftwl_free(array.V);
   fftwl_free(tmp2); 
   fftwl_free(tmp3); 
   fftwl_free(k);
   fftwl_free(extra.U);
   fftwl_free(extra.dU);
   fftwl_free(extra.dQ);
   fftwl_free(extra.dV);
   fftwl_free(extra.B);
   fftwl_free(extra.w);
   fftwl_free(extra.dq); 
#if RK4
   free_rk4();
#endif
}

void reinitialize_mesh() {
   if (!(array.Q = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(array.V = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);  
   if (!(tmp2 = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(tmp3 = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(k = fftwl_malloc(input.N*sizeof(long double)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.U = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.dU = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.dQ = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.dV = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.B = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.w = fftwl_malloc((input.N/2-1)*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.dq = fftwl_malloc((input.N)*sizeof(long double)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   //------------ moved init_rk4 from here
   FILE *fh = fopen("output.log","a");
   fprintf(fh, "Reinitialize Data: initialized FFTW plans with mode %d\n", FMODE);
   fclose(fh);
}

void reinitialize_mesh2() {
   free(tmp0); 
   free(tmp1);
   fftwl_destroy_plan(p0f);
   fftwl_destroy_plan(p1f);
   fftwl_destroy_plan(p2f);
   fftwl_destroy_plan(p0b);
   fftwl_destroy_plan(p1b);
   if (!(tmp0 = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(tmp1 = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   p0f = fftwl_plan_dft_1d(input.N, tmp0, tmp0, FFTW_FORWARD, FMODE);
   p1f = fftwl_plan_dft_1d(input.N, tmp1, tmp1, FFTW_FORWARD, FMODE);
   p2f = fftwl_plan_dft_1d(input.N, tmp2, tmp2, FFTW_FORWARD, FMODE);
   p0b = fftwl_plan_dft_1d(input.N, tmp0, tmp0, FFTW_BACKWARD, FMODE);
   p1b = fftwl_plan_dft_1d(input.N, tmp1, tmp1, FFTW_BACKWARD, FMODE);
#if RK4
   init_rk4();  // moved init_rk4 here
#endif
}


void move_mesh(long double uin, long double Lin, long int Nin) { // Lin = L/2 in Matlab
  long double qin = 2.0L*atan2l(Lin*sinl(0.5L*uin), cosl(0.5L*uin));
  long double beta = tanl(0.5L*(uin-input.u));
  long double q, qr;
  long double u0, u1, q0 = 2.0L*atan2l(input.L*sinl(0.5L*input.u), cosl(0.5L*input.u));
  long int N0 = input.N;

  memcpy(tmp0, array.Q, N0*sizeof(fftwl_complex));
  memcpy(tmp1, array.V, N0*sizeof(fftwl_complex));

  fftwl_execute(p0b);
  fftwl_execute(p1b);
  clean_up_mesh();

  input.N = Nin;

  reinitialize_mesh();
  memset(array.Q,0,sizeof(fftwl_complex)*input.N);
  memset(array.V,0,sizeof(fftwl_complex)*input.N);
  FILE *fh = fopen("move_mesh.txt","w");
  fprintf(fh, "# 0. q 1. u1 2. u2 3.-4. V0 5.-6. V1\n\n");
  for (long int j = 0; j < input.N; j++) {
    q  = PI*(2.0L*j/input.N - 1.0L); 
    u0 = input.u + 2.0L*atan2l(input.L*sinl(0.5L*(q-q0)),cosl(0.5L*(q-q0)));
    u1 = uin + 2.0L*atan2l(Lin*sinl(0.5L*(q-qin)),cosl(0.5L*(q-qin)));
    qr = PI+extra.q + 2.0L*atan2l(beta + Lin*tanl(0.5L*(q-qin)), input.L*(1.0L-Lin*beta*tanl(0.5L*(q-qin)) )   );
    for (long int l = 0; l < N0/2; l++) {
      array.Q[j] += tmp0[l]*cexpl(-1.IL*l*qr)/N0;
      array.V[j] += tmp1[l]*cexpl(-1.IL*l*qr)/N0;
    }
    fprintf(fh, "%.19LE\t%.19LE\t%.19LE\t%.19LE\t%.19LE\n", q, u0, u1, creall(array.V[j]), cimagl(array.V[j]));
  }
  fclose(fh);  
  reinitialize_mesh2();
  input.L = Lin; 
  input.u = uin;
  set_aux();

}


void save_data() {
  printf("Writing binary restart to:\n%s\n", input.rname);
  FILE *fh = fopen(input.rname, "wb");
  fwrite(&input, sizeof(params), 1, fh);
  fwrite(array.Q, sizeof(fftwl_complex), input.N, fh);
  fwrite(array.V, sizeof(fftwl_complex), input.N, fh);
  fclose(fh);
}

void load_data() {
  debug_msg("Load Data: initiated\n", EXIT_FALSE);
  FILE *fh = fopen(input.rname, "rb");

  if (fh != NULL) {
    size_t nr = fread(&rs_input, sizeof(params), 1, fh);
    if (rs_input.N == input.N) {
      //printf("%ld\t%ld\n", nr, rs_input.N);
      debug_msg("Load Data: calling initialize data\n", EXIT_FALSE);
      initialize_data();
      initialize_auxiliary_arrays();
      debug_msg("Initialize Data: complete\n", EXIT_FALSE);
      nr = fread(array.Q, sizeof(fftwl_complex), input.N, fh);
      nr = fread(array.V, sizeof(fftwl_complex), input.N, fh);
      fclose(fh);
      //move_mesh(array.Q, array.V, rs_input.u, rs_input.L);
      printf("Default Mesh Set From Restart (L,u) = (%.10LE,%.10Le)\n", rs_input.L, rs_input.u);
      input.L = rs_input.L;  
      input.u = rs_input.u;
      debug_msg("Load Data: restart arrays read successful\n", EXIT_FALSE);
      debug_msg("Load Data: complete\n", EXIT_FALSE);
    } else {
      fclose(fh);
      debug_msg("Load Data: restart arrays mismatch\n", EXIT_FALSE);
      debug_msg("Load Data: arrays sized to match input data\n", EXIT_FALSE);
      debug_msg("Load Data: calling initialize data\n", EXIT_FALSE);
      input.N = rs_input.N;
      input.L = rs_input.L;
      input.u = rs_input.u;
      initialize_data();
      nr = fread(array.Q, sizeof(fftwl_complex), input.N, fh);
      nr = fread(array.V, sizeof(fftwl_complex), input.N, fh);
      debug_msg("Load Data: complete\n", EXIT_FALSE);
      initialize_auxiliary_arrays();
      debug_msg("Initialize Data: complete\n", EXIT_FALSE);
    }
    printf("Name %s\nN = %ld\ng = %LE\ns = %LE\nl = %LE\nq = %LE\ntl = %LE\n", 
        input.rname, input.N, input.g, input.s, input.L, input.u, input.tl);
  } else {
    debug_msg("Load Data: restart file missing\nLoad Data: complete\n", EXIT_FALSE);
    printf("Restart File Missing:\tstarting new simulation.\n");
    initialize_data();
    initialize_auxiliary_arrays();
    debug_msg("Initialize Data: complete\n", EXIT_FALSE);
    debug_msg("Not using restart data in this run\n", EXIT_FALSE);
  }
  set_aux();
}

void dump_input() {
  FILE *fh = fopen("output.log","a");
  fprintf(fh, "Load Parameters: restart name = %s\n", input.rname);
  fprintf(fh, "Load Parameters: number of points = %ld\n", input.N);
  fprintf(fh, "Load Parameters: gravity g = %.15Le\n", input.g);
  fprintf(fh, "Load Parameters: surface tension s = %.15Le\n", input.s);
  fprintf(fh, "Load Parameters: transformation l = %.15Le\n", input.L);
  fprintf(fh, "Load Parameters: transformation u0 = %.15Le\n", input.u);
  fprintf(fh, "Load Parameters: refinement tolerance tol = %.15Le\n", input.tl);
  fprintf(fh, "Load Parameters: complete\n");
  fclose(fh);
}

void read_input(char *fname) {
  FILE *fh = fopen(fname,"r");
  char line[512], name[128], value[128];
  if (fh == NULL) {
    debug_msg("Load Parameters: input file missing\n", EXIT_FALSE);
    debug_msg("Load Parameters: complete\n", EXIT_TRUE);
  } else {
    while (fgets(line, 512, fh)!=NULL) {
      sscanf(line, "%s\t%s", name, value);
      if (strcmp(name,"resname=") == 0) sprintf(input.rname,"%s", value);
      if (strcmp(name,"npoints=") == 0) input.N = strtol(value, NULL, 10);
      if (strcmp(name,"gravity=") == 0) input.g = strtold(value, NULL);
      if (strcmp(name,"surface=") == 0) input.s = strtold(value, NULL);
      if (strcmp(name,"transfl=") == 0) input.L = strtold(value, NULL);
      if (strcmp(name,"transfu=") == 0) input.u = strtold(value, NULL);
      if (strcmp(name,"toleran=") == 0) input.tl = strtold (value, NULL);
    }
    dump_input(); 
  }
}

void load_parameters(int argc, char* argv[]) {
  FILE *fh = fopen("output.log","w");
  fprintf(fh, "Output Log File:\n\n");
  fprintf(fh, "Load Parameters: initiated\n");
  fclose(fh);
  if (argc == 1) {
    fh = fopen("output.log","a");
    fprintf(fh, "Load Parameters: requires name of file with parameters as a single input argument\n");
    fprintf(fh, "Load Parameters: complete\n");
    fclose(fh);
    exit(1);
  } else {
    char str[256];
    sprintf(str,"%s", argv[1]);
    fh = fopen("output.log","a");
    fprintf(fh, "Load Parameters: reading parameters from %s\n", str);
    fclose(fh);
    read_input(str);
  }
}

void initialize_auxiliary_arrays() {
   if (!(tmp0 = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Initialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(tmp1 = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Initialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(tmp2 = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Initialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(tmp3 = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Initialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(k = fftwl_malloc(input.N*sizeof(long double)))) debug_msg("Initialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   debug_msg("Initialize Data: memory allocation of auxiliary arrays successful\n", EXIT_FALSE);
   if (!(extra.U = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Initialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.dU = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Initialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.dQ = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.dV = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Reinitialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.B = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Initialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.w = fftwl_malloc((input.N/2-1)*sizeof(fftwl_complex)))) debug_msg("Initialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   if (!(extra.dq = fftwl_malloc((input.N)*sizeof(long double)))) debug_msg("Initialize Data: memory allocation failed\nInitialize Data: complete\n", EXIT_TRUE);
   //------------ 
   p0f = fftwl_plan_dft_1d(input.N, tmp0, tmp0, FFTW_FORWARD, FMODE);
   p1f = fftwl_plan_dft_1d(input.N, tmp1, tmp1, FFTW_FORWARD, FMODE);
   p2f = fftwl_plan_dft_1d(input.N, tmp2, tmp2, FFTW_FORWARD, FMODE);
   p0b = fftwl_plan_dft_1d(input.N, tmp0, tmp0, FFTW_BACKWARD, FMODE);
   p1b = fftwl_plan_dft_1d(input.N, tmp1, tmp1, FFTW_BACKWARD, FMODE);
   //------------
#if RK4
   init_rk4();
#endif

   FILE *fh = fopen("output.log","a");
   fprintf(fh, "Initialize Data: initialized FFTW plans with mode %d\n", FMODE);
   fclose(fh);
}

void free_arrays() {
   debug_msg("Free Data: initiated\n",EXIT_FALSE);
   free(array.Q);
   free(array.V);
   free(tmp0);    fftwl_destroy_plan(p0f);
   free(tmp1);    fftwl_destroy_plan(p1f);
   free(tmp2);    fftwl_destroy_plan(p0b);
   free(tmp3);    fftwl_destroy_plan(p1b);
   debug_msg("Free Data: complete\n",EXIT_FALSE);
}


void initialize_data() {
   debug_msg("Initialize Data: initiated\n", EXIT_FALSE);
   if (!(array.Q = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Initialize Data: memory allocation failed\nInitialize Data: complete\n", 1);
   if (!(array.V = fftwl_malloc(input.N*sizeof(fftwl_complex)))) debug_msg("Initialize Data: memory allocation failed\nInitialize Data: complete\n", 1);  
   debug_msg("Initialize Data: memory allocation successful\n", EXIT_FALSE);
}

void simulate() {
   long double u, q0; 
   ref_nog = 0;
   for (int j = 0; j < input.N; j++) {
     q = PI*(2.L*j/input.N - 1.0L);
     u = input.u + 2.L*atan2l(input.L*sinl(0.5L*(q-extra.q)), cosl(0.5L*(q-extra.q)));
     array.Q[j] = 1.L; 
     array.V[j] = -0.05IL*(1.L/ctanl(0.5L*(u-0.12IL)) - 1.IL) + 0.00IL*(1.L/ctanl(0.5L*(u-0.24IL)) - 1.IL); // run 8&9
   }  
   free_arrays();
}

void init_simulate() {
   init_rk4_module(&input, &motion, &array, &extra);
}










