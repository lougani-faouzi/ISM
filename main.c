#include <stdio.h>
#include <stdlib.h>
#include "simulation.h"


int main(int argc, char **argv)
{
  
  char *fname = argv[1];
  int N_step = atoi(argv[2]);
  int test = atoi(argv[3]);
  
  particles_t *donne = NULL;
  lennard_jones jon;
  vec_t *tv;
  
  double **r ;

  //read data in filename
  donne = read_data(fname);

  tv = translator_vector_init(N_sym);
  
  //distance
  r = compute_distance(donne);
  //lennard Jones Non periodique
  jon =compute_Lennard_Jones(donne , r);
  printf("Ulj = %lf\n", jon.Ulj);

  lennard_verify(jon, donne);

  print_forces(jon,donne);
  printf("Ulj = %lf\n", jon.Ulj);
  lennard_verify(jon, donne);
  
  print_sum_forces(jon, donne);
  

  return 0;
}

