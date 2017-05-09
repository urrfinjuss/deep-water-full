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
  long double q, u, xi;
  long double fi = 0.0L*PI;

  //long double Q = 0.01L; // sim 5
  //long double Q = 0.20L; // sim 6
  //long double Q = 0.001L; sim 10  // sim 8 and 9: Q = 2.5L;
  //long double C = -0.25L; sim 10  // sim 8 and 9: C = -0.02L;
  long double Q = 0.50L; // pirate run 2.5L
  long double C = -0.02L;
  // sim 10
  //fftwl_complex s1 = -1.0IL*ctanl(0.5L*0.050IL);
  //fftwl_complex s2 = -1.0IL*ctanl(0.5L*0.075IL);
  // sim 11
  fftwl_complex s1 = -1.0IL*ctanl(0.5L*0.050IL);
  fftwl_complex s2 = -1.0IL*ctanl(0.5L*0.075IL);
  //
  fftwl_complex Disc = -cpowl(s2-s1,2)*(1+Q*Q) + 2.L*Q*(s2-s1)*(1-s1*s2);
  fftwl_complex r1 = (1.IL*(s1+s2) - csqrtl(Disc))/(2.L-Q*(s2-s1));
  fftwl_complex r2 = (1.IL*(s1+s2) + csqrtl(Disc))/(2.L-Q*(s2-s1));
  //fftwl_complex D = csqrtl(-0.25L*cpowl(a2-a1,2L)+Q*(a2-a1));
  //fftwl_complex w1 = 0.5IL*(a1+a2) + D;
  //fftwl_complex w2 = 0.5IL*(a1+a2) - D;
  long double S = tanhl(0.5L*0.2L);


  printf("s1 = %.12LE + I%.12LE\n", creall(s1), cimagl(s1));
  printf("s2 = %.12LE + I%.12LE\n", creall(s2), cimagl(s2));
  printf("r1 = %.12LE + I%.12LE\n", creall(r1), cimagl(r1));
  printf("r2 = %.12LE + I%.12LE\n", creall(r2), cimagl(r2));
  printf("Discriminant = %.12LE + I%.12LE\n", creall(Disc), cimagl(Disc));
  long double a1 = 0.045L, b1 = 0.050L, c1 = 0.055L;
  long double a2 = 0.070L, b2 = 0.075L, c2 = 0.080L;
  a1 = 0.0050L; a2 = 0.0075L; // sim 11

  long double a3 = 0.2L; // sim 12

  for (long int j = 0; j < state.number_modes; j++) {
    q = 2.L*PI*(j*overN - 0.5L) - conf.origin_offset;
    u = conf.image_offset - fi + 2.L*atan2l(conf.scaling*sinl(0.5L*q),cosl(0.5L*q));
    xi = tanl(0.5L*u);
    // sims (proper periodic)
    data[0][j] = 2.L*(xi-1.IL*s1)*(xi-1.IL*s2)/((2.L-Q*(s2-s1))*(xi-r1)*(xi-r2));
    data[1][j] = C*(data[0][j] - 1.L);
    data[0][j] = csqrtl(data[0][j]);
    // sims (sim 12 lge pole)
    data[0][j] = 2.L/(2.L + Q*S)*(1.L + 1.IL*Q*(1.L - S*S)/((2.L+Q*S)*xi - 1.IL*(2.L*S+Q)));
    data[1][j] = C*(data[0][j] - 1.L);
    data[0][j] = csqrtl(data[0][j]);
    // sims (initial Z-tilde)
    data[0][j] =  clogl(  1.IL*cexpl(-1.IL*(u-1.IL*a3)) - 1.IL ) + 0.5IL*PI;  // sim 12
    //data[0][j] =  clogl(1.IL*csinl(0.5L*(u-1.IL*a1))) -  clogl(1.IL*csinl(0.5L*(u-1.IL*a2)));  // sim 11: pirate run
    data[0][j] = -1.IL*Q*data[0][j]/1.L; 	// necessary (sim 11)
    // pade test for VZi
    tmpc[0][j] = data[1][j]*overN; 
  }
  // if setting Z-tilde then uncomment below
  
  complex_array_out("Z0.txt", data[0]);
  convertZtoQ(data[0], data[0]);  // convert Z-tilde to Q
  
  // end uncomment
  complex_array_out("Q0.txt", data[0]);
  for (long int j = 0; j < state.number_modes; j++) {
    data[1][j] = C*(data[0][j]*data[0][j] - 1.L);  // required when setting Z-tilde
    tmpc[0][j] = data[1][j]*overN;
  }
  complex_array_out("V0.txt", data[1]);
  fftwl_complex z0;
  fftwl_execute(ift0);
  printf("Average V = %.12LE\n", cabsl(tmpc[0][0]));
  if (PADE_TEST) {
    for (long int j = 0; j < state.number_modes; j++) {
      data[0][j] = 1.L/(data[0][j]*data[0][j]);
      data[1][j] = -1.IL*data[1][j]*data[0][j];
    }
    optimal_pade("Zu_initial.pade", data[0]); 
    optimal_pade("Phiu_initial.pade", data[1]); 
    printf("PADE_TEST flag is on!\n");
    exit(1);
  }

  compute_zero_mode_complex(tmpc[0], 0.L, &z0);
  tmpc[0][0] = z0;
  //printf("Average V = %.12LE\n", cabsl(tmpc[0][0]));
  for (long int j = 0; j < state.number_modes; j++) tmpc[1][j] = data[0][j]*overN;
  fftwl_execute(ift1);
  for (long int j = 0; j < state.number_modes; j++) {
    if (cabsl(tmpc[0][j]) < 1.0E-14L) tmpc[0][j] = 0.L;
    if (cabsl(tmpc[1][j]) < 1.0E-14L) tmpc[1][j] = 0.L;
  }
  fftwl_execute(ft0);
  fftwl_execute(ft1);
  memcpy(data[1], tmpc[0], state.number_modes*sizeof(fftwl_complex));
  memcpy(data[0], tmpc[1], state.number_modes*sizeof(fftwl_complex));
  complex_array_out("Q1.txt", data[0]);
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
