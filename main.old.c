#include "header.h"

int main( int argc, char* argv[]) {
  init_lowlevel();
  load_parameters(argc, argv);
  init_timemarching();
  load_data();
  


  simulate();
  
}
