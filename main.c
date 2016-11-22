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
  
  complex_array_out("zt-original.txt", data[0]);
  convertZtoQ(data[0], data[0]);
  convertQtoZ(data[0], tmpc[5]);  
  //printf("Mapping 1: Potential Energy\t%.19LE\n", state.potentialE);
  complex_array_out("zt-recovered.txt", tmpc[5]);
  
  //complex_array_out("inQ.txt", data[0]);
  //complex_array_out("inV.txt", data[1]);
  
  //fftwl_complex 	kinetic_energy;
  //compute_hamiltonian(data[0], data[1], &kinetic_energy);
  restore_potential(data[0], data[1], tmpc[2]);
  print_constants();
  //printf("Mapping 1: Kinetic Energy\t%.19LE\n", state.kineticE);
  

  map new_map;
  new_map.scaling 	= 0.25L;
  new_map.image_offset 	= 1.00L;

  //complex_array_out("zt-original.txt", data[0]);
  remap(&new_map, 512); 
  //complex_array_out("zt-recovered.txt", data[0]);
  
  restore_potential(data[0], data[1], tmpc[3]);  
  print_constants();
  //convertQtoZ(data[0], tmpc[5]);  
  //printf("Mapping 2: Potential Energy\t%.19LE\n", state.potentialE);
  //printf("Mapping 2: Kinetic Energy\t%.19LE\n", state.kineticE);

  //complex_array_out("Phi.ph.txt", tmpc[3]); 
  for (long int j = 0; j < state.number_modes; j++) {
    tmpc[3][j] = conjl(tmpc[3][j])*tmpc[3][j];
  }
  project(tmpc[3], tmpc[2]);
  //complex_array_out("Proj.absPhi2.ph.txt", tmpc[2]); 
  
}
