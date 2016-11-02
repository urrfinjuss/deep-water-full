#include "header.h"

params 	state;
map 	conf;

int main( int argc, char* argv[]) {
  load_parameters(argc, argv);		// read configuration from file	 //
  init_memory();			// initialize memory unit 	 //
  allocate_memory();			// allocate all memory		 //
  set_mapping();			// set conformal mapping 	 //

  real_array_out("conf.dq.txt", conf.dq);

  unsigned int format_flag = 0;
  if (strcmp(state.txt_format,"ascii") == 0 )  format_flag = 1;
  if (strcmp(state.txt_format,"binary") == 0 ) format_flag = 2;
  if (strcmp(state.txt_format,"pade") == 0 )   format_flag = 3;
 
  switch (format_flag) {

    case 1:
      printf("Reading ASCII file\n");
      load_ascii();
      break;

    case 2:
      printf("Reading binary save\n");
      printf("Placeholder\n");
      exit(1);
      break;

    case 3:
      printf("Reading Pade data\n");
      printf("Placeholder\n");
      exit(1);
      break;

    default:
      printf("Unknown text format\n");
      exit(1);
  }

  backup_arrays();
  printf("Complete\n");
  
}
