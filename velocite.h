#ifndef _VELOCITY_H_
#define _VELOCITY_H_

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "simulation.h"



//struct moment component
typedef struct moment_s
{
  double mx;
  double my;
  double mz;

}moment_t;

//struct result moment and Temperature
typedef struct cinet_s
{
    moment_t *mi;
    int Ndl;
    double Ec;
    double Tc;
    
}cinet_t;

//function velocity verlet;
cinet_t compute_velocity_verlet(particles_t *p, vec_t *tv, int N);
//function init moment
cinet_t init_moment_cinetique(particles_t *p);
//energie cinetique
cinet_t compute_cinetique_energie(cinet_t cn, particles_t *p);
//recalibrated function
cinet_t compute_first_recalibrated(cinet_t cn, particles_t *p);
//second recalibrated
cinet_t compute_second_recalibrated(cinet_t cn, particles_t *p);

//thermostat berendsen
cinet_t compute_mc_berendsen(cinet_t cn , particles_t *p);
//free memory
void free_cinetique(cinet_t cn);
#endif
