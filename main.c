#include <stdio.h>
#include <stdlib.h>
#include "velocite.h"
#include "simulation.h"



int main(int argc, char **argv)
{

  //check argments
  if(argc < 2)
  {
    printf("Usage : [bin] [file] [STEP] [Berendsen]\n");
    exit(0);
  }

  char *fname = argv[1];
  ///int N_step = atoi(argv[2]);
  //int test = atoi(argv[3]);
  
  particles_t *p = NULL;
  //lennard_t js;
  vec_t *tv;
  cinet_t Cn;
  
  //double **r ;

  //read data in filename
  p = read_data(fname);

  //print data in teminal
  //print_data(p);
  tv = translator_vector_init(N_sym);
  //update particles
  /*p = update_particle_data(p, tv, N_sym);
  //distance
  r = compute_distance(p);
  //lennard Jones Non periodique
  js =compute_Lennard_Jones(p , r);
  printf("Ulj = %lf\n", js.Ulj);

  lennard_verify(js, p);
  //Part for compute  Lennard Jones periodical
  js = compute_Lennard_Jones_periodic(p, r, N_sym);

  //print forces
  //print_forces(js,p);
  printf("Ulj = %lf\n", js.Ulj);
  lennard_verify(js, p);*/
  
  //print sum forces
  //print_sum_forces(js, p);
  
  printf("*************Compute Dynamics molecule************\n");

  //compute verlet velocity
  Cn = compute_velocity_verlet(p, tv, N_sym);

  //simulation_dynamique_molecule(p,  cn, N_sym, N_step);
  
  free_cinetique(Cn);
  free_particle(p);
  return 0;
}
