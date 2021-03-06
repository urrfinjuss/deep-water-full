#include "header.h"

params 		state;
pade		pade_data;
map 		conf;
map 		alt_map;

int main( int argc, char* argv[]) {
  load_parameters(argc, argv);		// read configuration from file	 //
  init_memory();			// initialize memory unit 	 //
  init_timemarching();			// initialize timemarching 	 //
  init_pade();				// initialize rational approx	 //
  allocate_memory();			// allocate all memory		 //

  unsigned int format_flag = 0;
  if (strcmp(state.txt_format,"ascii") == 0 )  		format_flag = 1;
  if (strcmp(state.txt_format,"binary") == 0 ) 		format_flag = 2;
  if (strcmp(state.txt_format,"pade") == 0 )   		format_flag = 3;
  if (strcmp(state.txt_format,"none") == 0 )   		format_flag = 4;
  if (strcmp(state.txt_format,"jon_wilkening") == 0 )   format_flag = 5;
 
  switch (format_flag) {

    case 1:
      printf("Reading ASCII file\n");
      load_ascii();
      set_mapping();
      convertQtoZ(data[0], tmpc[5]);  
      restore_potential(data[0], data[1], tmpc[2]);
      break;

    case 2:
      printf("Reading binary save\n");
      printf("Placeholder\n");
      exit(1);
      break;

    case 3:
      printf("Reading Pade data\n");
      set_mapping();
      load_pade();
      break;

    case 4:
      printf("Starting new simulation\n");
      set_mapping();
      set_initial_data();
      convertQtoZ(data[0], tmpc[5]);  
      restore_potential(data[0], data[1], tmpc[2]);
      break;

    case 5: 
      printf("Setting data from Jon Wilkening file format.\n");
      set_mapping();
      set_initial_JW();
      convertQtoZ(data[0], tmpc[5]);  
      restore_potential(data[0], data[1], tmpc[2]);
      break;

    default:
      printf("Unknown text format\n");
      exit(1);
  }
  print_constants();
  evolve_rk6(data[0], data[1]);

}
