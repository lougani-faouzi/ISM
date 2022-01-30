#include "simulation.h"
#include <math.h>


lennard_jones compute_Lennard_Jones(particles_t *p, double **r)
{
  
  double sommefx = 0.0;
  double sommefy = 0.0;
  double sommefz = 0.0;

  lennard_jones jon;
  double tmp = 0.0;

  jon.frc=(force_t**)malloc(sizeof(force_t*)*p->N_particles_total);
  jon.som_frc=(force_t*)malloc(sizeof(force_t)*p->N_particles_total);
   
  for (int i = 0; i < p->N_particles_total; i++)
  {
      jon.frc[i]=(force_t*)malloc(sizeof(force_t)*p->N_particles_total);
  }
    
  //partie initialisation des forces  
  for (int i = 0; i < p->N_particles_total; i++)
  {
    for (int j = 0; j < p->N_particles_total; j++)
    {
      //initialisation des force pour chaque particule
      jon.frc[i][j].fx = 0.0;
      jon.frc[i][j].fy = 0.0;
      jon.frc[i][j].fz = 0.0;
    }

    //initialisation des sommes des forces 
    jon.som_frc[i].fx = 0.0;
    jon.som_frc[i].fy = 0.0;
    jon.som_frc[i].fz = 0.0;
  }

    for (int i = 0; i < p->N_particles_total; i++)
    {
      for (int j = i+1; j < p->N_particles_total; j++)
       {
        
           tmp +=((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)
           /(r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j]))
           - 2*((sigma*sigma*sigma*sigma*sigma*sigma)/(r[i][j] * r[i][j] *r[i][j]));

           //calculs des differentes forces avec la formule de cours 
           jon.frc[i][j].fx = - 48*eps*(((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)
                      /(r[i][j]*r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j]))
                      -((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)/
                      (r[i][j]*r[i][j]*r[i][j]*r[i][j])))*(p[i].coord.x - p[j].coord.x);
           jon.frc[i][j].fy = - 48*eps*(((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)
                      /(r[i][j]*r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j]))
                      -((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)/
                      (r[i][j]*r[i][j]*r[i][j]*r[i][j])))*(p[i].coord.y - p[j].coord.y);
           jon.frc[i][j].fz = - 48*eps*(((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)
                      /(r[i][j]*r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j]))
                      -((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)/
                      (r[i][j]*r[i][j]*r[i][j]*r[i][j])))*(p[i].coord.z - p[j].coord.z);

          //mise a jour de la force j avec i 
          jon.frc[j][i].fx =-jon.frc[i][j].fx;
          jon.frc[j][i].fy =-jon.frc[i][j].fy;
          jon.frc[j][i].fz =-jon.frc[i][j].fz;
         
         //mise a jour de la force i avec j 
        jon.frc[i][j].fx =jon.frc[i][j].fx;
        jon.frc[i][j].fy =jon.frc[i][j].fy;
        jon.frc[i][j].fz =jon.frc[i][j].fz;

     
          sommefx += jon.frc[i][j].fx;
          sommefy += jon.frc[i][j].fy;
          sommefz += jon.frc[i][j].fz;
       }

        //affecter les somme des forces 
        jon.som_frc[i].fx = sommefx;
        jon.som_frc[i].fy = sommefy;
        jon.som_frc[i].fz = sommefz;
      }
    
    //calulc de lennard jones terme avec la formule du cours 
    jon.Ulj = 4*eps*tmp;

  return jon;
}


particles_t *update_particle_data(particles_t *p, vec_t *tv, int N)
{

  particles_t *tmp = (particles_t*)malloc(sizeof(particles_t)*p->N_particles_total);
  tmp->N_particles_total = p->N_particles_total;

  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < tmp->N_particles_total; j++)
    {
      tmp[j].coord.x = p[j].coord.x + tv[i].x;
      tmp[j].coord.y = p[j].coord.y + tv[i].y;
      tmp[j].coord.z = p[j].coord.z + tv[i].z;
    }
  }
  
  return tmp;
}


lennard_jones compute_Lennard_Jones_periodic(particles_t *p, double **r, int N)
{
  /* calcul periodique de lennard jones */
  double sommefx = 0.0;
  double sommefy = 0.0;
  double sommefz = 0.0;
  double k = 0.5;

  lennard_jones jon;
  double tmp = 0.0;

  jon.frc=(force_t**)malloc(sizeof(force_t*)*p->N_particles_total);
  jon.som_frc=(force_t*)malloc(sizeof(force_t)*p->N_particles_total);

  for (int i = 0; i < p->N_particles_total; i++)
  {
      jon.frc[i]=(force_t*)malloc(sizeof(force_t)*p->N_particles_total);
  }
  
  
  for (int i = 0; i < p->N_particles_total; i++)
  {
    for (int j = 0; j < p->N_particles_total; j++)
    {
      /*initialisation des force pour chaque particule*/
      jon.frc[i][j].fx = 0.0;
      jon.frc[i][j].fy = 0.0;
      jon.frc[i][j].fz = 0.0;
    }

    /*initialisation des sommes des forces */
    jon.som_frc[i].fx = 0.0;
    jon.som_frc[i].fy = 0.0;
    jon.som_frc[i].fz = 0.0;
  }
  
  for (int i_sym = 0; i_sym < N; i_sym++)
  {
    for (int i = 0; i < p->N_particles_total; i++)
      {
        for (int j = i+1; j < p->N_particles_total; j++)
        {
          
          if (r[i][j] > R_cut*R_cut)
              continue;
          
            tmp +=((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)
           /(r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j]))
           - 2*((sigma*sigma*sigma*sigma*sigma*sigma)/(r[i][j] * r[i][j] *r[i][j]));

            /*calculs des differentes forces avec la formule de cours*/
           jon.frc[i][j].fx = - 48*eps*(((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)
                      /(r[i][j]*r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j]))
                      -((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)/
                      (r[i][j]*r[i][j]*r[i][j]*r[i][j])))*(p[i].coord.x - p[j].coord.x);
           jon.frc[i][j].fy = - 48*eps*(((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)
                      /(r[i][j]*r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j]))
                      -((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)/
                      (r[i][j]*r[i][j]*r[i][j]*r[i][j])))*(p[i].coord.y - p[j].coord.y);
           jon.frc[i][j].fz = - 48*eps*(((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)
                      /(r[i][j]*r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j] * r[i][j]))
                      -((sigma*sigma*sigma*sigma*sigma*sigma*sigma*sigma)/
                      (r[i][j]*r[i][j]*r[i][j]*r[i][j])))*(p[i].coord.z - p[j].coord.z);

          /*mise a jour de la force j avec i */
          jon.frc[j][i].fx =-jon.frc[i][j].fx;
          jon.frc[j][i].fy =-jon.frc[i][j].fy;
          jon.frc[j][i].fz =-jon.frc[i][j].fz;
         
          /*mise a jour de la force i avec j*/ 
          jon.frc[i][j].fx =jon.frc[i][j].fx;
          jon.frc[i][j].fy =jon.frc[i][j].fy;
          jon.frc[i][j].fz =jon.frc[i][j].fz;

       
          sommefx += jon.frc[i][j].fx;
          sommefy += jon.frc[i][j].fy;
          sommefz += jon.frc[i][j].fz;
       }

        /*affcter les sommes des forces*/ 
        jon.som_frc[i].fx = sommefx;
        jon.som_frc[i].fy = sommefy;
        jon.som_frc[i].fz = sommefz;
       }
  }
    
     if(N == 1)
        k = 1;

  /*calulc de lennard jones terme avec la formule du cours*/ 
  jon.Ulj = 4*k*eps*tmp;

  return jon;
}

//print forces for each particles
void print_forces(lennard_jones jon, particles_t *p)
{ 
  /* affichage des forces */
  for (int i = 0; i < p->N_particles_total; i++)
      {
        for (int j = 0; j < p->N_particles_total; j++)
        {
           printf("%lf  \t  %lf  \t  %lf", jon.frc[i][j].fx, jon.frc[i][j].fy, jon.frc[j][i].fz) ;
          
        }
        printf("\n");
       }
}


void print_sum_forces(lennard_jones jon, particles_t *p)
{ 
  /*affichage de la somme des forces */
  
  int i=0;
  printf("\a\a\a ## Somme des forces calculées ##\a\a\a \n");
  while(i < p->N_particles_total)
  {
       printf("Fx[%d] = %lf \t Fy[%d] = %lf  \t  Fz[%d] = %lf \n", i, jon.som_frc[i].fx , i, jon.som_frc[i].fy , i, jon.som_frc[i].fz); 
       i++;
  }
  
}

void lennard_verify(lennard_jones jon, particles_t *p)
{ 
  double Fx = 0.0;
  double Fy = 0.0;
  double Fz = 0.0;
  int i=0;
  int j=0;
  
  while (i < p->N_particles_total)
  {
    while(j < p->N_particles_total)
    {
      Fx=Fx+jon.frc[i][j].fx;
      Fy=Fy+jon.frc[i][j].fy;
      Fz=Fz+jon.frc[i][j].fz;
      j++;
    }
    i++;
  }

  printf("Fx = %lf  Fy = %lf  Fz = %lf\n", Fx,Fy,Fz);
}


void free_lennard(lennard_jones jon, particles_t *p)
{

  int i=0;
  while (i<p->N_particles_total)
  {
    	free(jon.frc[i]);
  	i++;
  }
  
  free(jon.frc);
  free(jon.som_frc);

}


/*read data for coord particles
  param file: name of file
  return : table of coordinate particles
*/
particles_t *read_data(char *fname)
{

  //check if file is null
  if (fname == NULL)
  {
    printf("Error read data %s\n",__func__);
    exit(0);
  }
  FILE *file=fopen(fname, "r+");
  //check file
  if(!file)
    exit(1);
  int nb_elmt =0;
  int c;
  
  //check the number of element in file
  do
  {
    if(c=='\n')
      nb_elmt++;
  } while((c=fgetc(file))!=EOF);

  particles_t *p = (particles_t*)malloc(sizeof(particles_t)*nb_elmt);
 

  rewind(file);
  fscanf(file,"%d",&p->N_particles_total);

  //check if nb_particles_total different number of coordinate
  if (p->N_particles_total != nb_elmt-1)
  {
    printf("Number patricle is different of the number coordinates\n" );
    exit(0);
  }

  for (int i = 0; i < p->N_particles_total; i++)
  {
    fscanf(file,"  %d   %lf   %lf   %lf",&p[i].type,&(p[i].coord.x), &(p[i].coord.y),&(p[i].coord.z));
    
  }
  return p;
}

/*Function for print data
input param: particles
*/
void print_data(particles_t *p)
{
  //check if pointer is NULL
  if(p == NULL)
  {
    printf("Error pointer is null\n" );
    exit(0);
  }
  printf("P  \t| Type  \t |X \t\t |Y \t\t |Z\n");

  for (int i = 0; i < p->N_particles_total; i++)
  {
    printf("%d \t %d \t %lf \t %lf \t %lf\n",i,p[i].type,p[i].coord.x, p[i].coord.y,p[i].coord.z);
  }

  printf("=======================================================\n");
}

/*compute distance different particles
 *params: particles
 *return distance for the differents particles
*/
double **compute_distance(particles_t *p)
{
  //check if pointer is NULL
  if(p == NULL)
  {
    printf("Error pointer is null\n" );
    exit(0);
  }

   double **r = NULL;
   r=(double**)malloc(p->N_particles_total*sizeof(double*));
   for (int i = 0; i < p->N_particles_total; i++)
     r[i] =(double*)malloc(p->N_particles_total*sizeof(double));


  for (int i = 0; i < p->N_particles_total; i++)
  {
    for (int j = i+1; j < p->N_particles_total; j++)
    {
        //compute distance
        r[i][j] = ((p[i].coord.x - p[j].coord.x)*(p[i].coord.x - p[j].coord.x))+
                 ((p[i].coord.y - p[j].coord.y)*(p[i].coord.y - p[j].coord.y))+
                 ((p[i].coord.z - p[j].coord.z)*(p[i].coord.z - p[j].coord.z));
        //update  distance
        r[j][i] = r[i][j];
    }
  }

  return r;
}

/*traslator vector initialisation
*/
vec_t *translator_vector_init(int N)
{
  vec_t *tv = (vec_t*)malloc(sizeof(vec_t)*N);

  for (int i = 0; i < N; i++)
  {
    tv[i].x = (double)((i / 9) - 1)*L;
    tv[i].y = (double)((i / 3)% 3 -1)*L;
    tv[i].z = (double)((i%3) - 1) *L;
  }
  
  return tv;
}

//free memory
void free_vector_translate( vec_t *tv)
{
  free(tv);
}

void free_distance(double **r, particles_t *p)
{
  for (int i = 0; i < p->N_particles_total; i++)
              free(r[i]);

  free(r);
}

void free_particle(particles_t *p)
{
  free(p);
}
