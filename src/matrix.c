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
first_iteration (struct radexi_data *rxi, struct radexi_results *rxi_res)
{
  for (unsigned int i = 0; i < rxi->mi.numof_radtr; i++)
    {
      unsigned int u = rxi->mi.radtr[i].ulev - 1;
      unsigned int l = rxi->mi.radtr[i].llev - 1;

      printf ("%d -> %d\n", u, l);

      float coef = hP * sol / kB * rxi->mi.radtr[i].xnu / rxi->mc_par.Tbg;
      if (coef >= 160)
        coef = 0;
      else 
        coef = 1 / (exp (coef) - 1);

      rxi_res->rates[u][u] += rxi->mi.radtr[i].a_einst * (1 + coef);
      rxi_res->rates[l][l] += rxi->mi.radtr[i].a_einst * \
                        rxi->mi.enlev[u].statw * coef / rxi->mi.enlev[l].statw;
      rxi_res->rates[u][l] -= rxi->mi.radtr[i].a_einst * \
                        rxi->mi.enlev[u].statw * coef / rxi->mi.enlev[l].statw;
      rxi_res->rates[l][u] -= rxi->mi.radtr[i].a_einst * (1 + coef);
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
subsequent_iterations (struct radexi_data *rxi, struct radexi_results *rxi_res)
{
  for (unsigned int i = 0; i < rxi->mi.numof_radtr; i++)
    {
      unsigned int u = rxi->mi.radtr[i].ulev - 1;
      unsigned int l = rxi->mi.radtr[i].llev - 1;

      float Snu = hP * sol / kB * rxi->mi.radtr[i].xnu / rxi_res->Tex[i];
      if (Snu >= 160)
        Snu = 0;
      else 
        Snu = 2 * hP * sol * pow (rxi->mi.radtr[i].xnu, 3) / \
            (exp (hP * sol / kB * rxi->mi.radtr[i].xnu / rxi_res->Tex[i]) - 1);

      rxi_res->tau[i] = rxi->mc_par.coldens / rxi->mc_par.line_width * \
                                (rxi_res->xpop[l] * rxi->mi.enlev[u].statw / \
                                rxi->mi.enlev[l].statw - rxi_res->xpop[u]) / \
                        (1.0645 * 8 * M_PI * pow (rxi->mi.radtr[i].xnu, 3) / \
                         rxi->mi.radtr[i].a_einst);

      float beta = esc_prob (rxi_res->tau[i], SPHERE);
      float coef = rxi->bg.intens[i] * beta / \
                   (2 * hP * sol) / pow (rxi->mi.radtr[i].xnu, 3);
 
      rxi_res->rates[u][u] += rxi->mi.radtr[i].a_einst * (beta + coef);
      rxi_res->rates[l][l] += rxi->mi.radtr[i].a_einst * \
                        rxi->mi.enlev[u].statw * coef / rxi->mi.enlev[l].statw;
      rxi_res->rates[u][l] -= rxi->mi.radtr[i].a_einst * \
                        rxi->mi.enlev[u].statw * coef / rxi->mi.enlev[l].statw;
      rxi_res->rates[l][u] -= rxi->mi.radtr[i].a_einst * (beta + coef);
    }
}

void
main_calculations (struct radexi_data *rxi, struct radexi_results *rxi_res)
{
  first_iteration(rxi, rxi_res);
  // Correct rates for collisional rates
  for (unsigned int i = 0; i < rxi->mi.numof_enlev; i++)
    {
      rxi_res->rates[i][i] += rxi->totrate[i];
      for (unsigned int j = 0; j < rxi->mi.numof_enlev; j++)
        {
          if (i != j)
            rxi_res->rates[i][j] -= rxi->crate[j][i];
        }
    }

  // Write rates to the new variable to preserve it in the future. Also test 
  // whether the matrix whould be reduced to exclude radiatively coupled 
  // levels.
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

  // Now separates collisionally coupled levels from those that are coupled 
  // by radiative processes, computing an effective cascade matrix for rates of
  // transfer from one low-lyinglevel to another.
  for (unsigned int i = 0; i < nred; i++)
    for (unsigned int j = 0; j < nred; j++)
      for (unsigned int k = nred + 1; k < rxi->mi.numof_enlev; k++)
        local_rates[i][j] += fabs (rxi_res->rates[k][j] * \
                              rxi_res->rates[i][k] / rxi_res->rates[k][k]); 

  // Solving reduced system of equations explicitly for the low-lying levels
  double lr[rxi->mi.numof_enlev * rxi->mi.numof_enlev];
  double rhs[rxi->mi.numof_enlev];
  unsigned int c = 0; 
  for (unsigned int i = 0; i < rxi->mi.numof_enlev; i++)
    {
      rhs[i] = 1e-30 * rxi->mc_par.total_density;

      for (unsigned int j = 0; j < rxi->mi.numof_enlev; j++)
        lr[c++] = local_rates[i][j];
    }
  
  gsl_matrix_view A = gsl_matrix_view_array (lr, rxi->mi.numof_enlev, rxi->mi.numof_enlev);
  gsl_vector_view b = gsl_vector_view_array (rhs, rxi->mi.numof_enlev);
  gsl_vector *x = gsl_vector_alloc (rxi->mi.numof_enlev);
  int s;


  gsl_permutation *p = gsl_permutation_alloc (rxi->mi.numof_enlev);
  gsl_linalg_LU_decomp (&A.matrix, p, &s);
  gsl_linalg_LU_solve (&A.matrix, p, &b.vector, x);

  printf("x = \n");
  gsl_vector_fprintf (stdout, x, "%g");

  gsl_permutation_free (p);
  gsl_vector_free (x);
}
