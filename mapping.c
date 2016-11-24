#include "header.h"

static fftwl_complex *saveQ, *saveV;

void set_mapping() {
  long double overN = 1.0L/state.number_modes;
  long double a = (1.0L - powl(conf.scaling, 2))/(1.0L + powl(conf.scaling, 2));
  long double b = 0.5L*(1.0L + powl(conf.scaling, 2))/conf.scaling;

  conf.origin_offset = 2.0L*atan2l(conf.scaling*sinl(0.5L*conf.image_offset), cosl(0.5L*conf.image_offset));
  for (long int j = 0; j < state.number_modes/2-1; j++) {
    conf.w[j] = cexpl(-1.0IL*(j+1)*(conf.origin_offset - 2.0IL*atanhl(conf.scaling)))*overN;
  }
  for (int j = 0; j < state.number_modes; j++) {
    long double q = 2.L*PI*(1.0L*j*overN - 0.5L);
    conf.dq[j] = b*(1.0L + a*cosl(q - conf.origin_offset));
  }
}

void remap(map_ptr new_map, unsigned long int N) {
  long double s = sinl(0.5L*new_map->image_offset);
  long double c = cosl(0.5L*new_map->image_offset);
  long double beta = tanl(0.5L*(new_map->image_offset - conf.image_offset));
  long double overN0 = 1.L/state.number_modes;
  long double overN = 1.L/N;
  unsigned long int N0 = state.number_modes;
  unsigned int QC_Q_pass = 0;
  unsigned int QC_V_pass = 0;
  unsigned int QC_pass = 0;

  /*
  map old_map;
  memcpy(&old_map, &conf, sizeof(map));
  printf("Old Map has N = %ld\n", N0);
  printf("Old scaling %23.16LE\n", old_map.scaling);
  */

  new_map->origin_offset = 2.0L*atan2l(new_map->scaling*s,c);
  memcpy(tmpc[0], data[0], N0*sizeof(fftwl_complex));
  memcpy(tmpc[1], data[1], N0*sizeof(fftwl_complex));
  fftwl_execute(ift0);
  fftwl_execute(ift1);
  saveQ = fftwl_malloc(N0*sizeof(fftwl_complex));
  saveV = fftwl_malloc(N0*sizeof(fftwl_complex));
  memcpy(saveQ, tmpc[0], N0*sizeof(fftwl_complex));
  memcpy(saveV, tmpc[1], N0*sizeof(fftwl_complex));
  /*for (long int j = 0; j < state.number_modes/2; j++) {
    tmpc[0][j] = tmpc[0][j]*overN0;
    tmpc[1][j] = tmpc[1][j]*overN0;
  }
  complex_array_out("preref.Q.ft.txt", tmpc[0]);
  complex_array_out("preref.V.ft.txt", tmpc[1]);*/
  deallocate_memory();
  state.number_modes = N;
  allocate_memory(); 
  
  if (N0 < state.number_modes) {
    unsigned long int mem_offset = state.number_modes-N0/2;
    memcpy(tmpc[0], saveQ, (N0/2)*sizeof(fftwl_complex));
    memcpy(tmpc[1], saveV, (N0/2)*sizeof(fftwl_complex));
    memcpy(tmpc[0]+mem_offset, saveQ+N0/2, (N0/2)*sizeof(fftwl_complex));
    memcpy(tmpc[1]+mem_offset, saveV+N0/2, (N0/2)*sizeof(fftwl_complex));
  } else if (N0 > state.number_modes) {
    unsigned long int mem_offset = N0 - state.number_modes/2;
    memcpy(tmpc[0], saveQ, (state.number_modes/2)*sizeof(fftwl_complex));
    memcpy(tmpc[1], saveV, (state.number_modes/2)*sizeof(fftwl_complex));
    memcpy(tmpc[0]+state.number_modes/2, saveQ+mem_offset, (state.number_modes/2)*sizeof(fftwl_complex));
    memcpy(tmpc[1]+state.number_modes/2, saveV+mem_offset, (state.number_modes/2)*sizeof(fftwl_complex));
  } else {
    memcpy(tmpc[0], saveQ, state.number_modes*sizeof(fftwl_complex));
    memcpy(tmpc[1], saveV, state.number_modes*sizeof(fftwl_complex));
  }
  long double 		q, q_r, q0 = PI + conf.origin_offset;
  fftwl_complex 	w = 1.L;
  for (long int j = 0; j < state.number_modes; j++) {
    q   = new_map->scaling*tanl(1.L*PI*(j*overN - 0.5L) - 0.5L*new_map->origin_offset);
    q_r = q0 + 2.0L*atan2l(beta+q, conf.scaling*(1.0L - beta*q));
    w = cexpl(-1.IL*q_r);
    data[0][j] = tmpc[0][state.number_modes/2-1];
    data[1][j] = tmpc[1][state.number_modes/2-1];
    for (long int l = state.number_modes/2-1; l > 0; l--) {
      data[0][j] = data[0][j]*w + tmpc[0][l-1];
      data[1][j] = data[1][j]*w + tmpc[1][l-1];

    } 
    data[0][j] = data[0][j]*overN0;
    data[1][j] = data[1][j]*overN0;
  }
  conf.scaling = new_map->scaling;
  conf.image_offset = new_map->image_offset;
  set_mapping();
  
  // verify that new map passes QC 
  memcpy(tmpc[0], data[0], state.number_modes*sizeof(fftwl_complex));
  memcpy(tmpc[1], data[1], state.number_modes*sizeof(fftwl_complex));
  fftwl_execute(ift0);
  fftwl_execute(ift1);
  for (long int j = 0; j < state.number_modes/2; j++) {
    tmpc[0][j] = tmpc[0][j]*overN;
    tmpc[1][j] = tmpc[1][j]*overN;
  }
  memset(tmpc[0]+state.number_modes/2, 0, state.number_modes/2*sizeof(fftwl_complex));
  memset(tmpc[1]+state.number_modes/2, 0, state.number_modes/2*sizeof(fftwl_complex));
  map_quality(tmpc[0], &QC_Q_pass);
  map_quality(tmpc[1], &QC_V_pass);
  //complex_array_out("postref.Q.ft.txt", tmpc[0]);
  //complex_array_out("postref.V.ft.txt", tmpc[1]);
  QC_pass = MIN(QC_Q_pass, QC_V_pass);
  printf("Conformal map is %u\n", QC_pass);
  if (QC_pass) {
    fftwl_free(saveQ);
    fftwl_free(saveV); 
  } else {
    
    printf("Bad Quality Map.\nPlaceholder\n");
    exit(1);
  }

}

void map_quality(fftwl_complex *in, unsigned int *QC_pass) {
  long double full_sum 		= 0.0L;
  long double partial_sum	= 0.0L;
  long double qc_ratio		= 1.0L;

  for (long int j = state.number_modes/2-1; j > -1; j--) {
    full_sum += cabsl(in[j]);
    if (j == state.number_modes*7/16) partial_sum = full_sum;
  }
  qc_ratio = partial_sum/sqrtl(1.L + powl(full_sum, 2));
  printf("Full Sum    \t%.19LE\n", full_sum);
  printf("Partial Sum \t%.19LE\n", partial_sum);
  printf("Partial/Full\t%.19LE\n", qc_ratio);
  if (qc_ratio < 1.0E-15L) {
	printf("QC Pass\n");
	*QC_pass = 1;
  } else printf("QC Fail\n");
}








