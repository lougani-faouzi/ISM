#ifndef _SIMULATION_H_
#define _SIMULATION_H_

struct particles
{
  int nb_elemt;
  int *type;
  double *x;
  double *y;
  double *z;
  double *fx;
  double *fy;
  double *fz;
  
};
typedef struct particles particles_t;
void affichage(particles_t *);
particles_t *lecture(char *);
double **distance_particules(particles_t *);
void Lennard_Jones(particles_t *, double **);
particles_t *Forces(particles_t *, double **);
void forces_nulles(particles_t *);

#endif 
