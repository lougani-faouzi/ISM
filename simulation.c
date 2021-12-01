#include "simulation.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


particles_t *lecture(char *fname)
{
  // lecture apartir d'un fichier de coordonnée et allocation de la structure 
  int val,val2;
  int type;
  double x,y,z;
  int nb_elemt =0;
  int read;
  char temp[1001];

  FILE *file=fopen(fname, "r");
  
  do
  {
    if(read=='\n') nb_elemt++;
  } while((read=fgetc(file))!= EOF);

  particles_t *p = (particles_t*)malloc(sizeof(particles_t)*nb_elemt);
  rewind(file);
  fscanf(file,"%d   %d",&val,&val2);
  
  p->nb_elemt=nb_elemt-1;
  p->type = aligned_alloc(64,nb_elemt*sizeof(p));
  
  //allication de l espace memoire pour les coordonnées x,y,z des  particules 
  p->x = aligned_alloc(64,nb_elemt*sizeof(p));
  p->y = aligned_alloc(64,nb_elemt*sizeof(p));
  p->z = aligned_alloc(64,nb_elemt*sizeof(p));
  
  p->fx = aligned_alloc(64,nb_elemt*sizeof(p));
  p->fy = aligned_alloc(64,nb_elemt*sizeof(p));
  p->fz = aligned_alloc(64,nb_elemt*sizeof(p));
  
  //initialisation a 0 des forces 
  for (int i = 0; i <nb_elemt ; i++)
  {
    p->fx[i] = 0.0;
    p->fy[i] = 0.0;
    p->fz[i] = 0.0;
  }
  
  for (int i = 0; i < nb_elemt; i++)
  {
    fscanf(file,"%d  %lf   %lf   %lf",&type,&x,&y,&z);
    p->type[i]=type;
    p->x[i]=x;
    p->y[i]=y;
    p->z[i]=z;
  
  }
  
  return p;

}

//affichage des coordonnees
void affichage(particles_t *p)
{
  //check if pointer is NULL
  printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
  printf("Question 1\n");
  printf("PR  \t t  \t x \t\t y \t\t z\n");

  for (int i = 0; i < p->nb_elemt; i++)
  {
   printf("%d \t %d \t %lf \t %lf \t %lf\n",i,p->type[i],p->x[i],p->y[i],p->z[i]);
  }

  printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
}


// calculer les distances entres les particlues et la retourner sous forme d'un tabelau 
double **distance_particules(particles_t *p)
{ 
   double **dist =(double**)malloc(p->nb_elemt*sizeof(double*));
   double **tmp=(double**)malloc(p->nb_elemt*sizeof(double*));

  for (int i = 0; i < p->nb_elemt; i++)
  {
    dist[i] =malloc(p->nb_elemt*sizeof(double));
    tmp[i] =malloc(p->nb_elemt*sizeof(double));
    for (int j = 0; j < p->nb_elemt; j++)
    {
      if(i!=j)
      {
             if(i < j)
             {
        
             tmp[i][j]=    sqrt(((p->x[i] - p->x[j])*(p->x[i] - p->x[j])) +
                                ((p->y[i] - p->y[j])*(p->y[i] - p->y[j])) + 
                                ((p->z[i] - p->z[j])*(p->z[i] - p->z[j])));
             }else if(i>j)
             {
             tmp[i][j] = tmp[j][i];
             }

      dist[i][j]=tmp[i][j];
      }
     
    }
  }
 
  return dist;
}

  

void Lennard_Jones(particles_t *p, double **dist)
{
  // calculer le potentiel lennard jones avec les distances r en paramettres et les particules 
  int i,j;
  double Ulj = 0.0; 
  double val = 0.0;
  double sigma= 0.3;
  double eps=0.2;

  while (i<p->nb_elemt)
  {
    while (j<p->nb_elemt)
    {           // faut fusionnner pour alller plus rapidment 
      if (i<j) // ici il faut supprimer la condition 
      {
        val=val+(pow((sigma/dist[i][j]),12) - 2*(pow((sigma/dist[i][j]),6)));
      }
     j++;  
    }
   i++;
  }
  
  Ulj = 4*eps*val;
  printf("Question 2\n");
  printf("\n le potentiel Lennard_Jones ULJ=%lf  \n ",Ulj);
  printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
}


void forces_nulles(particles_t *p)
{
  double cpt=0.0;
  int i=0;
  while (i<p->nb_elemt)
  {
    cpt= cpt+ p->fx[i] + p->fy[i] +p->fz[i];
    i++;
  }
  if(cpt== 0)
    printf("la sommes des forces agissant sur les particules est nulles\n");
  else
    printf("la sommes des forces agissant sur les particules n est pas nulles = %lf\n", cpt);
    
}

particles_t *Forces(particles_t *p, double **dist)
{
 // calculer les force en prenant en paramettre les distances deja calculees et les particules 
  double sigma=0.3;
  for (int i = 0; i < p->nb_elemt; i++)
  {
    for (int j = 0; j < p->nb_elemt; j++)
    {
      if(i >j)
      {
        p->fx[i] =p->fx[i] +( -48*sigma*((pow((sigma/dist[i][j]),14))-(pow((sigma/dist[i][j]),8)))*(p->x[i] - p->x[j]));
        p->fy[i] =p->fy[i] +( -48*sigma*((pow((sigma/dist[i][j]),14))-(pow((sigma/dist[i][j]),8)))*(p->y[i] - p->y[j]));
        p->fz[i] =p->fz[i] +( -48*sigma*((pow((sigma/dist[i][j]),14))-(pow((sigma/dist[i][j]),8)))*(p->z[i] - p->z[j]));
      }
    }
  }
  
 
  printf("Question 3\n");
 

  return p;
}







