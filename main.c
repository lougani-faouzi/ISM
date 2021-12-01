#include "simulation.h"
#include <stdio.h>
#include <math.h>

int main(){
  printf("\n TD1 lougani FAOUZI m2 IHPS -module ISM \n");
 
  //lecture des donnees 
  particles_t *p=lecture("particule.xyz");
  
  // affichage 
  affichage(p);
  
   
  //on calcul le potentiel 
  Lennard_Jones(p,distance_particules(p));
  
  
  //calculer les forces et verifier si le systeme est stable ou pas  
  forces_nulles(Forces(p,distance_particules(p)));
  

  return 0;
}
