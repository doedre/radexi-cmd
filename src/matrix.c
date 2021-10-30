/*	matrix.c
 *
 * All the main calculation is done here.
 *
 * ----------------------------------------------------------------------
 *
 * This file is part of the radexi-cmd program - C language implementation of 
 * RADEX software package. It is used to calculate molecular exitation and 
 * radiative transfer in a homogenuous medium. The purpose of this full 
 * refactoring is to speed up the process of calculation and implement some new
 * features like databases. It will make the program more stable and usable, 
 * which is crucial for newcomers.
 *
 * Documentation for radexi-cmd program is posted at 
 * https://github.com/doedre/radexi-cmd
 *
 * ----------------------------------------------------------------------
 *
 * From the RADEX software package:
 *    Documentation for the RADEX program is posted at
 *	  https://personal.sron.nl/~vdtak/radex/index.shtml 
 *
 *	  Although this program has been thoroughly tested, the authors do not 
 *	  claim that it is free of errors and gives correct results in all 
 *	  situations.
 *	
 *	  Publications using this program should make a reference to the paper:
 *	  A&A 468, 627 (2007).
 *
 * ---------------------------------------------------------------------*/

#include "radexi.h"

#include <math.h>
#include <gsl/gsl_linalg.h>

void
first_iteration (struct rxi_data *rxi)
{
  for (unsigned int i = 0; i < rxi->mi.numof_radtr; i++)
    {
      unsigned int u = rxi->rad_transfer[i].ulev;
      unsigned int l = rxi->rad_transfer[i].llev;

      float coef = hP * sol / kB * rxi->rad_transfer[i].xnu / rxi->mc.Tbg;
      if (coef >= 160)
        coef = 0;
      else 
        coef = 1 / (exp (coef) - 1);

      printf ("----preparing [%d][%d]th element\n", u, l);
      double uu = gsl_matrix_get (rxi->res.rates, u, u) +     \
                  rxi->rad_transfer[i].a_einst * (1 + coef);

      double ll = gsl_matrix_get (rxi->res.rates, l, l) +     \
                  rxi->rad_transfer[i].a_einst *              \
                  rxi->energy_level[u].statw * coef /       \
                  rxi->energy_level[l].statw;

      double ul = gsl_matrix_get (rxi->res.rates, u, l) -     \
                  rxi->rad_transfer[i].a_einst *              \
                  rxi->energy_level[u].statw * coef /       \
                  rxi->energy_level[l].statw;

      double lu = gsl_matrix_get (rxi->res.rates, l, u) -     \
                  rxi->rad_transfer[i].a_einst * (1 + coef);

      gsl_matrix_set (rxi->res.rates, u, u, uu);
      gsl_matrix_set (rxi->res.rates, l, l, ll);
      gsl_matrix_set (rxi->res.rates, u, l, ul);
      gsl_matrix_set (rxi->res.rates, l, u, lu);
    }
}

static float
esc_prob (const float tau, enum Geometry g)
{
  float beta = 0;

  if (g == SPHERE)
    {
      if (fabs (tau / 2) < 0.1)
        beta = (1.5 / tau) - (3 / pow (tau, 3)) + \
                    ((2 / tau) + (2 / pow (tau, 2))) * 1.5 * exp (-tau) / tau;
    }
  else if (g == SLAB)
    {
    }
  else if (g == LVG)
    {
    }

  return beta;
}

void
subsequent_iterations (struct rxi_data *rxi)
{
  for (unsigned int i = 0; i < rxi->mi.numof_radtr; i++)
    {
      unsigned int u = rxi->rad_transfer[i].ulev - 1;
      unsigned int l = rxi->rad_transfer[i].llev - 1;

      /* Calculating source function  */
      float Snu = hP * sol / kB * rxi->rad_transfer[i].xnu / rxi->res.Tex[i];
      if (Snu >= 160)
        Snu = 0;
      else 
        Snu = 2 * hP * sol * pow (rxi->rad_transfer[i].xnu, 3) /  \
              (exp (hP * sol / kB * rxi->rad_transfer[i].xnu /    \
              rxi->res.Tex[i]) - 1);

      rxi->res.tau[i] = rxi->mc.coldens / rxi->mc.line_width *                \
                      (rxi->res.pop[l] * rxi->energy_level[u].statw /         \
                      rxi->energy_level[l].statw - rxi->res.pop[u]) /         \
                      (1.0645 * 8 * M_PI * pow (rxi->rad_transfer[i].xnu, 3)/ \
                      rxi->rad_transfer[i].a_einst);

      float beta = esc_prob (rxi->res.tau[i], SPHERE);
      float coef = rxi->bg.intens[i] * beta /                       \
                   (2 * hP * sol) / pow (rxi->rad_transfer[i].xnu, 3);

      double uu = gsl_matrix_get (rxi->res.rates, u, u) +     \
                  rxi->rad_transfer[i].a_einst * (beta + coef);

      double ll = gsl_matrix_get (rxi->res.rates, l, l) +     \
                  rxi->rad_transfer[i].a_einst *              \
                  rxi->energy_level[u-1].statw * coef /       \
                  rxi->energy_level[l-1].statw;

      double ul = gsl_matrix_get (rxi->res.rates, u, l) -     \
                  rxi->rad_transfer[i].a_einst *              \
                  rxi->energy_level[u-1].statw * coef /       \
                  rxi->energy_level[l-1].statw;

      double lu = gsl_matrix_get (rxi->res.rates, l, u) -     \
                  rxi->rad_transfer[i].a_einst * (beta + coef);

      gsl_matrix_set (rxi->res.rates, u, u, uu);
      gsl_matrix_set (rxi->res.rates, l, l, ll);
      gsl_matrix_set (rxi->res.rates, u, l, ul);
      gsl_matrix_set (rxi->res.rates, l, u, lu);
    }
}

void
main_calculations (struct rxi_data *rxi)
{
  printf ("--making first iteration...\n");
  first_iteration(rxi);

  printf ("--correcting rates for collisional rates...\n");
  // Correct rates for collisional rates
  for (unsigned int i = 0; i < rxi->n_el; i++)
    {
      double ii = gsl_matrix_get (rxi->res.rates, i, i) + \
                  gsl_vector_get (rxi->mi.Ctot, i);
      gsl_matrix_set (rxi->res.rates, i, i, ii);
      for (unsigned int j = 0; j < rxi->n_el; j++)
        {
          if (i != j)
            {
              double ij = gsl_matrix_get (rxi->res.rates, i, j) - \
                          gsl_matrix_get (rxi->mi.C, j, i);
              gsl_matrix_set (rxi->res.rates, i, j, ij);
            }
        }
    }

  // Write rates to the new variable to preserve it in the future. Also test 
  // whether the matrix whould be reduced to exclude radiatively coupled 
  // levels.
/*
  float local_rates[rxi->mi.numof_enlev][rxi->mi.numof_enlev];
  float redcrit = 10 * rxi->mc_par.Tkin / hP / sol * kB;
  unsigned int nred = 0;
  for (unsigned int i = 0; i < rxi->mi.numof_enlev; i++)
    {
      for (unsigned int j = 0; j < rxi->mi.numof_enlev; j++)
          local_rates[i][j] = rxi_res->rates[i][j];

      if (rxi->mi.enlev[i].term <= redcrit)
        nred++;
    }
*/
  // Now separates collisionally coupled levels from those that are coupled 
  // by radiative processes, computing an effective cascade matrix for rates of
  // transfer from one low-lyinglevel to another.
/*  for (unsigned int i = 0; i < nred; i++)
    for (unsigned int j = 0; j < nred; j++)
      for (unsigned int k = nred + 1; k < rxi->mi.numof_enlev; k++)
        local_rates[j][i] += fabs (rxi_res->rates[k][i] * \
                              rxi_res->rates[j][k] / rxi_res->rates[k][k]); 
*/
  
  // Solving reduced system of equations explicitly for the low-lying levels
  printf ("--solving equations...\n");
  gsl_vector *b = gsl_vector_alloc (rxi->mi.numof_enlev);
  gsl_vector_set_all (b, 0);
  gsl_vector_set (b, b->size-1, 1);

  for (unsigned int i = 0; i < rxi->n_el; i++)
    {
      for (unsigned int j = 0; j < rxi->n_el; j++)
        printf ("%.2e  ", gsl_matrix_get (rxi->res.rates, i, j));
    printf ("\n");
    }

  gsl_linalg_HH_svx (rxi->res.rates, b);

  printf ("result=\n");
  gsl_vector_fprintf (stdout, b, "%e");

  gsl_vector_free (b);
  
  /*gsl_matrix_view A = gsl_matrix_view_array (lr, rxi->mi.numof_enlev, rxi->mi.numof_enlev);*/
  /*gsl_vector_view b = gsl_vector_view_array (rhs, rxi->mi.numof_enlev);*/
  /*gsl_vector *x = gsl_vector_alloc (rxi->mi.numof_enlev);*/
  /*int s;*/


  /*gsl_permutation *p = gsl_permutation_alloc (rxi->mi.numof_enlev);*/
  /*gsl_linalg_LU_decomp (&A.matrix, p, &s);*/
  /*gsl_linalg_LU_solve (&A.matrix, p, &b.vector, x);*/

  /*printf("x = \n");*/
  /*gsl_vector_fprintf (stdout, x, "%e");*/

  /*gsl_permutation_free (p);*/
  /*gsl_vector_free (x);*/
}
