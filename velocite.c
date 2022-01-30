#include "velocite.h"
#include<time.h>
#include<math.h>

#define sign_function(x, y) (y < 0.0 ? -x : x)

/*
*velocity verlet algorithm
*
*/
cinet_t compute_velocity_verlet(particles_t *p, vec_t *tv, int N)
{

    double **r;

    lennard_jones js;

    cinet_t cn;

    cn.mi = (moment_t*)malloc(sizeof(moment_t)*p->N_particles_total);
    particles_t *ps = (particles_t*)malloc(sizeof(particles_t)*p->N_particles_total);

    //update particles
    p = update_particle_data(p, tv, N_sym);
    //distance
    r = compute_distance(p);
    printf("r10 = %lf\n",r[1][0]);
    /*Part for compute  Lennard Jones periodical*/
    js = compute_Lennard_Jones_periodic(p, r, N_sym);

    printf("Ulj = %lf\n", js.Ulj);

    printf("js.som_frc[1].fx = %lf\n", js.som_frc[1].fx);
    //initilization moment
    for (int i = 0; i < p->N_particles_total; i++)
    {
        cn.mi[i].mx = dt * CONVERSION_FORCE * js.som_frc[i].fx /2.0;
        cn.mi[i].my = dt * CONVERSION_FORCE * js.som_frc[i].fy / 2.0;
        cn.mi[i].mz = dt * CONVERSION_FORCE * js.som_frc[i].fz / 2.0;
    }

    printf("cn.mi[1].mx = %lf\n", cn.mi[1].mx);
    
    ps->N_particles_total = p->N_particles_total;
    //update position
    for (int i = 0; i < ps->N_particles_total ; i++)
    {
        ps[i].coord.x = dt * cn.mi[i].mx / M_i;
        ps[i].coord.y = dt * cn.mi[i].my / M_i;
        ps[i].coord.z = dt * cn.mi[i].mz / M_i;
    }
    //update distance compute
    r = compute_distance(ps);
    printf("r10 = %lf\n",r[1][0]);
    //uptdate force 
    js = compute_Lennard_Jones_periodic(ps, r, N_sym);

    printf("Ulj = %lf\n", js.Ulj);

    //update moments
    for (int i = 0; i < p->N_particles_total; i++)
    {
        cn.mi[i].mx = dt * CONVERSION_FORCE * js.som_frc[i].fx ;
        cn.mi[i].my = dt * CONVERSION_FORCE * js.som_frc[i].fy / 2.0;
        cn.mi[i].mz = dt * CONVERSION_FORCE * js.som_frc[i].fz / 2.0;
    }

    printf("cn.mi[1].mx = %lf\n", cn.mi[1].mx);
    //init other parameter
    cn.Ndl = 0;
    cn.Ec = 0.0;
    cn.Tc = 0.0;
    //free
    free_distance(r,p);
    free_lennard(js,p);

    return cn;
}

/*initialisation moment cinetique
*params: particles p
*/
cinet_t init_moment_cinetique(particles_t *p)
{
    cinet_t cn;
    double c, s;

    cn.mi = (moment_t*)malloc(sizeof(moment_t) * p->N_particles_total);

    srand(time(NULL));

    for (int i = 0; i < p->N_particles_total; i++)
    {
        //component in x
        c = (double)rand() / (double)RAND_MAX;
        s = (double)rand() / (double)RAND_MAX;
        cn.mi[i].mx = sign_function(1.0, 0.5 - s) * c;

        //component in y
        c = (double)rand() / (double)RAND_MAX;
        s = (double)rand() / (double)RAND_MAX;
        cn.mi[i].my = sign_function(1.0, 0.5 - s) * c;

        //component in z
        c = (double)rand() / (double)RAND_MAX;
        s = (double)rand() / (double)RAND_MAX;
        cn.mi[i].mz = sign_function(1.0, 0.5 - s) * c;
    }
    
    return cn;
}

/*
*compute energie cinetique
*params: in moment
*/
cinet_t compute_cinetique_energie(cinet_t cn, particles_t *p)
{
    double tmp = 0;

    cn.Ndl = 3 * p->N_particles_total -3;

    for (int i = 0; i < p->N_particles_total; i++)
    {
        tmp += (cn.mi[i].mx * cn.mi[i].mx + cn.mi[i].my * cn.mi[i].my + cn.mi[i].mz * cn.mi[i].mz)/M_i;
    }
    
    cn.Ec = tmp / (2* CONVERSION_FORCE);
    
    //Temperature cinetique
    cn.Tc = cn.Ec / (cn.Ndl * CONSTANT_R);

    return cn;
}

/*
*first recalibrated 
*/
cinet_t compute_first_recalibrated(cinet_t cn, particles_t *p)
{
    double Ra;

    Ra = cn.Ndl * CONSTANT_R * T0 / cn.Ec;

    for (int i = 0; i < p->N_particles_total; i++)
    {
        //component in x
        cn.mi[i].mx = Ra * cn.mi[i].mx;

        //component in y
        cn.mi[i].my = Ra * cn.mi[i].my;

        //component in z
        cn.mi[i].mz = Ra * cn.mi[i].mz;
    }
    
    return cn;
}
/*
*second recalibrated
*/
cinet_t compute_second_recalibrated(cinet_t cn, particles_t *p)
{
    double Px = 0, Py = 0, Pz =0;
    //correction
    for (int i = 0; i < p->N_particles_total; i++)
    {
        Px += cn.mi[i].mx;
        Py += cn.mi[i].my;
        Pz += cn.mi[i].mz;
    }
    
    for (int i = 0; i < p->N_particles_total; i++)
    {
        //component in x
        cn.mi[i].mx = cn.mi[i].mx - Px;

        //component in y
        cn.mi[i].my = cn.mi[i].my - Py;

        //component in z
        cn.mi[i].mz = cn.mi[i].mz - Pz;
    }

    printf("Px = %lf, Py = %lf, Pz = %lf \n", Px, Py, Pz);

    return cn;
}

/*
*compute thermostat berendsen
*params: moments particles
*/
cinet_t compute_mc_berendsen(cinet_t cn, particles_t *p)
{
    for (int i = 0; i < p->N_particles_total; i++)
    {
        //component in x
        cn.mi[i].mx += gama*((T0/cn.Tc) - 1)*cn.mi[i].mx;

        //component in y
        cn.mi[i].my = gama*((T0/cn.Tc) - 1)*cn.mi[i].my;

        //component in z
        cn.mi[i].mz = gama*((T0/cn.Tc) - 1)*cn.mi[i].mz;
    }

    return cn;
}
//free memory after usage
void free_cinetique(cinet_t cn)
{
    free(cn.mi);
}
