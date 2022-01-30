#ifndef _LENNARDJONES_H_
#define _LENNARDJONES_H_
#include <string.h>
#include <stdio.h>
#include <stdlib.h>




//difine different value parameter
#define sigma 3.0
#define eps 0.2
#define L   30
#define R_cut 10.0
#define N_sym 27
#define dt  1
#define CONVERSION_FORCE  0.0001*4.186
#define M_i 18
#define CONSTANT_R  0.00199
#define T0  300
#define gama   0.01


//struct for the differents forces
typedef struct force_s
{
  double fx;
  double fy;
  double fz;
}force_t;

//struct for lennard jones
typedef struct lennard
{
  /* structure pour lennard jones */
  
  double Ulj;
  force_t **frc;
  force_t *som_frc;
  
}lennard_jones;



/*stock the coordinate x y z*/
typedef struct coord_s
{
  double x;
  double y;
  double z;
}coord_t;

//struct for stock the differents value for particles
typedef struct particles_s
{
  int N_particles_total;
  int N_particles_local;
  int type;
  coord_t coord;
}particles_t;


//tanslator vector
typedef struct vec_s
{
  double x;
  double y;
  double z;

}vec_t;

lennard_jones compute_Lennard_Jones(particles_t *p, double **r);
void lennard_verify(lennard_jones jon, particles_t *p);
void print_forces(lennard_jones jon, particles_t *p);
void print_sum_forces(lennard_jones jon, particles_t *p);
lennard_jones compute_Lennard_Jones_periodic(particles_t *p, double **r, int N);
void free_lennard(lennard_jones jon, particles_t *p);





//funtion to return struct particles
particles_t *read_data(char *fname);
//function to print the data of particles
void print_data(particles_t *p);
//function to compute the distance rij
double **compute_distance(particles_t *p);
//funtion transaltor vector
vec_t *translator_vector_init(int N);
//function update distance
particles_t *update_particle_data(particles_t *p, vec_t *tv, int N);
//free memory after usage
void free_vector_translate( vec_t *tv);
void free_distance(double **r, particles_t *p);
void free_particle(particles_t *p);


#endif
