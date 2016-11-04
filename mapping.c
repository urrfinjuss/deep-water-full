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

  new_map->origin_offset = 2.0L*atan2l(new_map->scaling*s,c);
  memcpy(tmpc[0], data[0], N0*sizeof(fftwl_complex));
  memcpy(tmpc[1], data[1], N0*sizeof(fftwl_complex));
  fftwl_execute(ift0);
  fftwl_execute(ift1);
  saveQ = fftwl_malloc(N0*sizeof(fftwl_complex));
  saveV = fftwl_malloc(N0*sizeof(fftwl_complex));
  memcpy(saveQ, tmpc[0], N0*sizeof(fftwl_complex));
  memcpy(saveV, tmpc[1], N0*sizeof(fftwl_complex));
  deallocate_memory();
  state.number_modes = N;
  allocate_memory(); 
  
  if (N0 < state.number_modes) {
    unsigned long int mem_offset = state.number_modes-N0/2;
    memcpy(tmpc[0], saveQ, (N0/2)*sizeof(fftwl_complex));
    memcpy(tmpc[1], saveV, (N0/2)*sizeof(fftwl_complex));
    memcpy(tmpc[0]+mem_offset, saveQ, (N0/2)*sizeof(fftwl_complex));
    memcpy(tmpc[1]+mem_offset, saveV, (N0/2)*sizeof(fftwl_complex));
  } else if (N0 > state.number_modes) {
    unsigned long int mem_offset = N0 - state.number_modes/2;
    memcpy(tmpc[0], saveQ, (state.number_modes/2)*sizeof(fftwl_complex));
    memcpy(tmpc[1], saveV, (state.number_modes/2)*sizeof(fftwl_complex));
    memcpy(tmpc[0]+mem_offset, saveQ, (state.number_modes/2)*sizeof(fftwl_complex));
    memcpy(tmpc[1]+mem_offset, saveV, (state.number_modes/2)*sizeof(fftwl_complex));
  } else {
    memcpy(tmpc[0], saveQ, state.number_modes*sizeof(fftwl_complex));
    memcpy(tmpc[1], saveV, state.number_modes*sizeof(fftwl_complex));
  }
  fftwl_free(saveQ);
  fftwl_free(saveV); 
  long double q, q_r;
  for (long int j = 0; j < state.number_modes; j++) {
    q   = 2.L*PI*(j*overN - 0.5L) - new_map->origin_offset;
    q_r = PI + conf.origin_offset + 2.0L*atan2l(beta+(new_map->scaling)*tanl(0.5L*q), conf.scaling*(1.0L - new_map->scaling*beta*tanl(0.5L*q)));
    for (long int l = N0/2-1; l > -1; l--) {
      data[0][j] += tmpc[0][l]*cexpl(-1.IL*l*q_r)*overN0; 
      data[0][j] += tmpc[1][l]*cexpl(-1.IL*l*q_r)*overN0;
    } 
  }
  conf.scaling = new_map->scaling;
  conf.image_offset = new_map->image_offset;
  set_mapping();
}









