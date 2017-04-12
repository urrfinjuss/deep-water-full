#include "header.h"


void load_parameters(int argc, char* argv[]) {
  if (argc != 2) {
    printf("Usage %s filename\n", argv[0]);
    exit(1);
  } else read_input(argv[1]);
}

void read_input(char *fname) {
  FILE *fh = fopen(fname,"r");
  char line[512], name[128], value[128];
  if (fh == NULL) {
    printf("Configuration file missing\n");
    exit(1);
  } else {
    while (fgets(line, 512, fh)!=NULL) {
      sscanf(line, "%s\t%s", name, value);
      if (strcmp(name,"resname=") == 0) sprintf(state.restart_name,"%s", value);
      if (strcmp(name,"txtasci=") == 0) sprintf(state.txt_format,"%s", value);
      if (strcmp(name,"npoints=") == 0) state.number_modes = strtol(value, NULL, 10);
      if (strcmp(name,"gravity=") == 0) state.gravity = strtold(value, NULL);
      if (strcmp(name,"surface=") == 0) state.surface_tension = strtold(value, NULL);
      if (strcmp(name,"transfl=") == 0) conf.scaling = strtold(value, NULL);
      if (strcmp(name,"transfu=") == 0) conf.image_offset = strtold(value, NULL);
      if (strcmp(name,"toleran=") == 0) state.tolerance = strtold (value, NULL);
      if (strcmp(name,"timesim=") == 0) state.final_time = strtold (value, NULL);
      if (strcmp(name,"n_poles=") == 0) state.number_poles = atol(value);
    }
  }
}

void set_initial_data() {
  long double overN = 1.L/state.number_modes;
  long double q, u;
  long double fi = 0.0L*PI;

  long double Q = 0.01L; // sim 5
  //long double Q = 0.20L; // sim 6
  long double C = -0.1L;
  long double a1 = 0.05L;
  long double a2 = 0.10L;
  fftwl_complex D = csqrtl(-0.25L*cpowl(a2-a1,2L)+Q*(a2-a1));
  fftwl_complex w1 = 0.5IL*(a1+a2) + D;
  fftwl_complex w2 = 0.5IL*(a1+a2) - D;

  printf("a1 = %.12LE + I%.12LE\n", creall(a1), cimagl(a1));
  printf("a2 = %.12LE + I%.12LE\n", creall(a2), cimagl(a2));
  printf("D = %.12LE + I%.12LE\n", creall(D), cimagl(D));
  printf("w1 = %.12LE + I%.12LE\n", creall(w1), cimagl(w1));
  printf("w2 = %.12LE + I%.12LE\n", creall(w2), cimagl(w2));

  for (long int j = 0; j < state.number_modes; j++) {
    q = 2.L*PI*(j*overN - 0.5L) - conf.origin_offset;
    u = conf.image_offset - fi + 2.L*atan2l(conf.scaling*sinl(0.5L*q),cosl(0.5L*q));
    //data[0][j] = -sinl(fi)+0.125L*sinl(2.L*fi)-33.0IL/64.L + 1.IL*cexpl(-1.IL*u) + 0.125IL*cexpl(-2.IL*u); // easy Z-tilde
    //data[0][j] = -1.IL*cpowl(0.01L*(1.L/ctanl(0.5L*(u-0.1IL)) - 1.IL),2);
    data[0][j] = 1.L;
    //data[0][j] = 1.L + 0.5L*cexpl(-1.IL*u); // set Q directly
    //data[1][j] = 0.L*cexpl(-1.IL*u);
    //data[1][j] = -0.0005IL*(1.L/ctanl(0.5L*(u-0.03IL)) - 1.IL); 
    data[1][j]  = -0.005IL*(1.L/ctanl(0.5L*(u-0.1IL)) - 1.IL); 
    //data[1][j] += +0.005IL*(1.L/ctanl(0.5L*(u-0.11IL)) - 1.IL); 
    // sims 3 (dipole surf + V)
    data[0][j] = 0.04L/(ctanl(0.5L*u)-0.04IL) - 0.04L/(ctanl(0.5L*u)-0.08IL);
    data[0][j] = cpowl(1.L - 0.5IL*(1.L + cpowl(ctanl(0.5L*u),2))*data[0][j], -0.5L);
    data[1][j] = -0.01IL*(1.L/ctanl(0.5L*(u - 0.8IL))  - 1.IL);
    // sim 5&6
    // formula for Q,V
    data[0][j] = 1.L + 0.5L*(Q*(a2-a1)/D)*(1.L/(2.L*tan(u/2) - w1) - 1.L/(2.L*tan(u/2)-w2)); // R
    data[1][j] = C/Q*(data[0][j] - 1.L);
    data[0][j] = csqrtl(data[0][j]);
    //
    tmpc[0][j] = data[1][j]*overN; 
  }
  fftwl_complex z0;
  fftwl_execute(ift0);
  complex_array_out("V0.txt", data[1]);
  printf("Average V = %.12LE\n", cabsl(tmpc[0][0]));
  compute_zero_mode_complex(tmpc[0], 0.L, &z0);
  tmpc[0][0] = z0;
  printf("Average V = %.12LE\n", cabsl(tmpc[0][0]));
  fftwl_execute(ft0);
  memcpy(data[1], tmpc[0], state.number_modes*sizeof(fftwl_complex));
  complex_array_out("Q0.txt", data[0]);
  complex_array_out("V1.txt", data[1]);
}

void load_ascii() {
  FILE *fh = fopen(state.restart_name, "r");
  if (fh) {
    char line[512], *v[5];
    long int		N;
    long double		T;
    for (int j = 0; j < 5; j++) v[j] = malloc(512);
    if (fgets(line, 512, fh) != NULL);
    if (fgets(line, 512, fh) != NULL) sscanf(line, "# N = %s\tL = %s\tu = %s\tT = %s", v[0], v[1], v[2], v[3]);
    if (fgets(line, 512, fh) != NULL);

    N = strtol(v[0], NULL, 10);
    T = strtold(v[3], NULL);

    conf.scaling = strtold(v[1], NULL);
    conf.image_offset = strtold(v[2], NULL);

    printf("Restart:\ntime = %.19LE\nN modes = %ld\n", T, N); 
    printf("Conformal L = %.19LE\nConformal u = %.19LE\n", conf.scaling, conf.image_offset);
    if (N != state.number_modes) {
      printf("Incompatible Grids\nPlaceholder\n");
      exit(1);
    }
    int counter = 0;
    while (fgets(line, 512, fh) != NULL) {
      if (counter == state.number_modes) {
	printf("Something is wrong with restart file\n"); 
	exit(1);
      }
      sscanf(line, "%s\t%s\t%s\t%s\t%s", v[0], v[1], v[2], v[3], v[4]);
      data[0][counter] = strtold(v[1], NULL)-strtold(v[0], NULL) + 1.0IL*strtold(v[2],NULL);
      data[1][counter] = strtold(v[3], NULL) + 1.0IL*strtold(v[4],NULL);
      counter++;
    }
    fclose(fh);    
    if (counter != state.number_modes) {
	printf("Something is wrong with restart file\n");
        exit(1);
    }
    for (int j = 0; j < 5; j++) free(v[j]);
  } else {
    printf("Missing restart file\n");
    exit(1);
  }
}
