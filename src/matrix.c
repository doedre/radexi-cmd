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

/* Calculations begin with this function, as it takes starting parameters and 
 * converts them to make first statistical equilibrium equation */
static void
starting_conditions (struct rxi_data *rxi)
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

      double uu = gsl_matrix_get (rxi->res.rates, u, u) +     \
                  rxi->rad_transfer[i].a_einst * (1 + coef);

      double ll = gsl_matrix_get (rxi->res.rates, l, l) +     \
                  rxi->rad_transfer[i].a_einst *              \
                  rxi->energy_level[u].statw * coef /         \
                  rxi->energy_level[l].statw;

      double ul = gsl_matrix_get (rxi->res.rates, u, l) -     \
                  rxi->rad_transfer[i].a_einst *              \
                  rxi->energy_level[u].statw * coef /         \
                  rxi->energy_level[l].statw;

      double lu = gsl_matrix_get (rxi->res.rates, l, u) -     \
                  rxi->rad_transfer[i].a_einst * (1 + coef);

      gsl_matrix_set (rxi->res.rates, u, u, uu);
      gsl_matrix_set (rxi->res.rates, l, l, ll);
      gsl_matrix_set (rxi->res.rates, u, l, ul);
      gsl_matrix_set (rxi->res.rates, l, u, lu);
    }
}

/* The main goal of RADEX is to speed up calculation time for nonLTE case. This
 * is done with the assumtion of homogenous medium and three different mediums 
 * are included in this function as the escape probability assumption.  */
static float
esc_prob (const double tau, enum Geometry g)
{
  double beta = 0;
  double tau_rad = tau / 2;

  if (g == SPHERE)
    {
      if (fabs (tau_rad) < 0.1)
       beta = 1 - 0.75 * tau_rad + 0.4 * pow (tau_rad, 2) -   \
              pow (tau_rad, 3) / 6 + pow (tau_rad, 4) / 17.5;
      else if (fabs (tau_rad) > 50)
        beta = 0.75 / tau_rad;
      else 
        beta = 0.75 / tau_rad * (1 - 1 / (2 * pow (tau_rad, 2)) + \
               (1 / tau_rad + 1 / (2 * pow (tau_rad, 2))) * exp (-2 * tau_rad));
    }
  else if (g == SLAB)
    {
    }
  else if (g == LVG)
    {
    }

  return beta;
}

/* Using results of first iteration to renew starting conditions and make them
 * more precise */
static int
refresh_starting_conditions (struct rxi_data *rxi)
{
  int thick_lines = 0;
  for (unsigned int i = 0; i < rxi->mi.numof_radtr; i++)
    {
      unsigned int u = rxi->rad_transfer[i].ulev;
      unsigned int l = rxi->rad_transfer[i].llev;
      double xt = pow (rxi->rad_transfer[i].xnu, 3);

/*
      [> Calculating source function  <]
      float Snu = fk * rxi->rad_transfer[i].xnu /   \
                  gsl_vector_get (rxi->res.Tex, i);
      if (Snu >= 160)
        Snu = 0;
      else 
        Snu = 2 * hP * sol * pow (rxi->rad_transfer[i].xnu, 3) /  \
              (exp (fk * rxi->rad_transfer[i].xnu /               \
                    gsl_vector_get (rxi->res.Tex, i)) - 1);
*/
      rxi->res.tau[i] = rxi->mc.coldens / rxi->mc.line_width / 1e5 *          \
            (gsl_vector_get (rxi->res.pop, l) * rxi->energy_level[u].statw /  \
            rxi->energy_level[l].statw - gsl_vector_get (rxi->res.pop, u)) /  \
            (1.0645 * 8 * M_PI * xt / rxi->rad_transfer[i].a_einst);

      if (rxi->res.tau[i] > 1e-2) thick_lines++;

      float beta = esc_prob (rxi->res.tau[i], SPHERE);
      float coef = rxi->bg.intens[i] * beta /                 \
                   (2 * hP * sol) / pow (rxi->rad_transfer[i].xnu, 3);

      const double uu = 
                  gsl_matrix_get (rxi->res.rates, u, u) +     \
                  rxi->rad_transfer[i].a_einst * (beta + coef);

      const double ll = 
                  gsl_matrix_get (rxi->res.rates, l, l) +     \
                  rxi->rad_transfer[i].a_einst *              \
                  rxi->energy_level[u].statw * coef /         \
                  rxi->energy_level[l].statw;

      const double ul = 
                  gsl_matrix_get (rxi->res.rates, u, l) -     \
                  rxi->rad_transfer[i].a_einst *              \
                  rxi->energy_level[u].statw * coef /         \
                  rxi->energy_level[l].statw;

      const double lu = 
                  gsl_matrix_get (rxi->res.rates, l, u) -     \
                  rxi->rad_transfer[i].a_einst * (beta + coef);

      gsl_matrix_set (rxi->res.rates, u, u, uu);
      gsl_matrix_set (rxi->res.rates, l, l, ll);
      gsl_matrix_set (rxi->res.rates, u, l, ul);
      gsl_matrix_set (rxi->res.rates, l, u, lu);
    }
  return thick_lines;
}

void
main_calculations (struct rxi_data *rxi)
{
  unsigned int iter = 0;
  int thick_lines = 1;
  double stop_condition = 0;
  // Old level populations
  gsl_vector *prev_pop = gsl_vector_calloc (rxi->n_el);
  gsl_vector *b = gsl_vector_calloc (rxi->mi.numof_enlev);
  do 
    {
      gsl_matrix_set_all (rxi->res.rates, -1e-30 * rxi->mc.total_density);
      printf ("--making %u iteration...\n", iter);
      if (iter == 0)
        starting_conditions (rxi);
      else 
        thick_lines = refresh_starting_conditions (rxi);

      stop_condition = 0;
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

      printf ("--solving equations...\n");
      gsl_vector_set_all (b, 0);
      gsl_vector_set (b, b->size-1, 1);
      for (unsigned int i = 0; i < rxi->n_el; i++)
        gsl_matrix_set (rxi->res.rates, rxi->res.rates->size1-1, i, 1);

      gsl_linalg_HH_svx (rxi->res.rates, b);

      // Calculating total population
      printf ("--calculating total population...\n");
      double total_pop = 0;
      for (unsigned int i = 0; i < rxi->n_el; i++)
        total_pop += gsl_vector_get (b, i);

      printf ("total pop: %.2e\n", total_pop);

      /* Keep old populations for underrelaxation */
      if (iter != 0)
        gsl_vector_memcpy (prev_pop, rxi->res.pop);

      /* Store new level populations  */
      printf ("--store new level populations...\n");
      for (unsigned int i = 0; i < rxi->n_el; i++)
        {
          const double new_pop_i = gsl_vector_get (b, i) / total_pop;

          gsl_vector_set (rxi->res.pop,
                          i,
                          (fabs (new_pop_i) > 1e-20) ? new_pop_i : 1e-20);
        }

      printf ("pop=\n");
      gsl_vector_fprintf (stdout, rxi->res.pop, "%.2e");

      /* Calculating excitation temperatures for lines  */
      printf ("--calculating excitation temperatures...\n");
      for (unsigned int i = 0; i < rxi->n_rt; i++)
        {
          unsigned int u = rxi->rad_transfer[i].ulev;
          unsigned int l = rxi->rad_transfer[i].llev;

          double new_Tex_i = fk * rxi->rad_transfer[i].xnu /            \
                             log (gsl_vector_get (rxi->res.pop, l) *    \
                                  rxi->energy_level[u].statw /          \
                                  gsl_vector_get (rxi->res.pop, u) /    \
                                  rxi->energy_level[l].statw);

          if (iter == 0)
            {
              gsl_vector_set (rxi->res.Tex, i, new_Tex_i);
              stop_condition = 1;
              continue;
            }

          if (rxi->res.tau[i] > 0.01)
            stop_condition += fabs ((gsl_vector_get (rxi->res.Tex, i) - new_Tex_i) / new_Tex_i);

          gsl_vector_set (rxi->res.Tex, i, 0.5 *              \
                (new_Tex_i + gsl_vector_get(rxi->res.Tex, i)));

          rxi->res.tau[i] = rxi->mc.coldens / rxi->mc.line_width *            \
            (gsl_vector_get (rxi->res.pop, l) * rxi->energy_level[u].statw /  \
            rxi->energy_level[l].statw - gsl_vector_get (rxi->res.pop, u)) /  \
                    (1.0645 * 8 * M_PI * pow (rxi->rad_transfer[i].xnu, 3) /  \
                    rxi->rad_transfer[i].a_einst);
        }

      /* Doing underrelaxation  */
      printf ("--doing underrelaxation...\n");
      /*if (iter != 0)
        for (unsigned int i = 0; i < rxi->n_rt; i++)
          {
            double new_pop_i = 0.3 * gsl_vector_get (rxi->res.pop, i) +   \
                               0.7 * gsl_vector_get (prev_pop, i);
            gsl_vector_set (rxi->res.pop, i, new_pop_i);
          }
*/
      printf ("Tex=\n");
      gsl_vector_fprintf (stdout, rxi->res.Tex, "%.2e");
      iter++;

      printf ("\n thick_lines = %d, stop_condition = %.3e\n", thick_lines, stop_condition);
    } while (thick_lines != 0 && stop_condition >= 1e-6);

  gsl_vector_free (prev_pop);
  gsl_vector_free (b);
}
