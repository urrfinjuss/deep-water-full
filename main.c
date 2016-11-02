#include "header.h"

params state;

int main( int argc, char* argv[]) {

  state.N = 256;

  //load_parameters(argc, argv);
  init_memory();
  allocate_memory();
  backup_arrays();
  printf("Complete\n");
  
}
