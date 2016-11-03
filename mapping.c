#include "header.h"

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

