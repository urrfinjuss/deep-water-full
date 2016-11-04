#include "header.h"

params 	state;
map 	conf;

int main( int argc, char* argv[]) {
  load_parameters(argc, argv);		// read configuration from file	 //
  init_memory();			// initialize memory unit 	 //
  allocate_memory();			// allocate all memory		 //

  unsigned int format_flag = 0;
  if (strcmp(state.txt_format,"ascii") == 0 )  format_flag = 1;
  if (strcmp(state.txt_format,"binary") == 0 ) format_flag = 2;
  if (strcmp(state.txt_format,"pade") == 0 )   format_flag = 3;
  if (strcmp(state.txt_format,"none") == 0 )   format_flag = 4;
 
  switch (format_flag) {

    case 1:
      printf("Reading ASCII file\n");
      load_ascii();
      set_mapping();
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

    case 4:
      printf("Starting new simulation\n");
      set_mapping();
      set_initial_data();
      break;

    default:
      printf("Unknown text format\n");
      exit(1);
  }


  real_array_out("conf.dq.txt", conf.dq);
  complex_array_out("zread.txt", data[0]);
  complex_array_out("vread.txt", data[1]);
  /*
  fftwl_execute(ift0);
  for (long int j = 0; j < state.number_modes/2; j++) {
    tmpc[0][j] = -1.0IL*j*tmpc[0][j]/state.number_modes;
    tmpc[0][state.number_modes-j-1] = 1.0IL*(j+1)*tmpc[0][state.number_modes-j-1];
  }
  complex_array_out("zq_direct.txt", tmpc[0]);
  */
  
  
  complex_array_out("zt-original.txt", data[0]);
  convertZtoQ(data[0], data[1]);
  
  //  Correct answer is here:
  memcpy(tmpc[0], data[0], state.number_modes*sizeof(fftwl_complex));
  // 
 
  complex_array_out("qread.txt", data[1]);
  convertQtoZ(data[1], data[0]);  
  complex_array_out("zt-recovered.txt", data[0]);

  backup_arrays();
  printf("Complete\n");
  
}
