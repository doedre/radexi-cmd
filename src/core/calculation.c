/**
 * @file core/calculation.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include "core/calculation.h"

#include "rxi_common.h"
#include "core/background.h"
#include "utils/database.h"
#include "utils/debug.h"

static void
set_starting_conditions (struct rxi_calc_data *data, const int n_radtr)
{
  DEBUG ("Set starting conditions; matrix size: %zu", data->einst->size1);
  for (int i = 0; i < n_radtr; ++i)
    {
      const unsigned int u = data->up[i] - 1;
      const unsigned int l = data->low[i] - 1;

      DEBUG ("%d -> %d", u, l);

      double coef = RXI_FK
                    * (gsl_vector_get (data->term, u)
                      - gsl_vector_get (data->term, l))
                    / data->temp_bg[0];

      if (coef >= 160)
        coef = 0;
      else 
        coef = 1 / (exp (coef) - 1);

      const double uu = gsl_matrix_get (data->rates, u, u)
                        + gsl_matrix_get (data->einst, u, l) * (1 + coef);

      const double ll = gsl_matrix_get (data->rates, l, l)
                        + gsl_matrix_get (data->einst, u, l)
                          * gsl_vector_get (data->weight, u) * coef
                          / gsl_vector_get (data->weight, l);

      const double ul = gsl_matrix_get (data->rates, u, l)
                        - gsl_matrix_get (data->einst, u, l)
                          * gsl_vector_get (data->weight, u) * coef
                          / gsl_vector_get (data->weight, l);

      const double lu = gsl_matrix_get (data->rates, l, u)
                        - gsl_matrix_get (data->einst, u, l) * (1 + coef);

      gsl_matrix_set (data->rates, u, u, uu);
      gsl_matrix_set (data->rates, l, l, ll);
      gsl_matrix_set (data->rates, u, l, ul);
      gsl_matrix_set (data->rates, l, u, lu);
    }
}

static int
refresh_starting_conditions (struct rxi_calc_data *data, const int n_radtr)
{
  int thick_lines = 0;
  gsl_matrix_set_all (data->rates, 1e-30);
  for (int i = 0; i < n_radtr; ++i)
    {
      const unsigned int u = data->up[i] - 1;
      const unsigned int l = data->low[i] - 1;

      double energy = gsl_vector_get (data->term, u)
                      - gsl_vector_get (data->term, l);
      const double tau = rxi_calc_optical_depth (data->col_dens[0],
          data->line_width[0], energy,
          gsl_matrix_get (data->einst, u, l),
          gsl_vector_get (data->weight, u), gsl_vector_get (data->weight, l),
          gsl_vector_get (data->pop, u), gsl_vector_get (data->pop, l));

      gsl_matrix_set (data->tau, u, l, tau);
      if (tau > 1e-2)
        ++thick_lines;

      const double beta = rxi_calc_escape_prob (tau, data->geom[0]);

      const double coef =
              (gsl_matrix_get (data->bgfield, u, l) * beta)
          / //---------------------------------------------
               (2 * RXI_HP * RXI_SOL * gsl_pow_3 (energy));

      DEBUG ("%d %d: coef = %.3e | beta %.3e", u, l, coef, beta);

      const double uu = gsl_matrix_get (data->rates, u, u)
                        + gsl_matrix_get (data->einst, u, l) * (beta + coef);

      const double ll = gsl_matrix_get (data->rates, l, l)
                        + gsl_matrix_get (data->einst, u, l) 
                          * gsl_vector_get (data->weight, u)
                          / gsl_vector_get (data->weight, l) * coef;

      const double ul = gsl_matrix_get (data->rates, u, l)
                        - gsl_matrix_get (data->einst, u, l) 
                          * gsl_vector_get (data->weight, u)
                          / gsl_vector_get (data->weight, l) * coef;

      const double lu = gsl_matrix_get (data->rates, l, u)
                        - gsl_matrix_get (data->einst, u, l) * (beta + coef);

      gsl_matrix_set (data->rates, u, u, uu);
      gsl_matrix_set (data->rates, l, l, ll);
      gsl_matrix_set (data->rates, u, l, ul);
      gsl_matrix_set (data->rates, l, u, lu);
    }
  return thick_lines;
}

RXI_STAT
rxi_calc_data_init (struct rxi_calc_data *calc_data,
                    const struct rxi_input_data *inp_data,
                    const struct rxi_db_molecule_info *mol_info)
{
  DEBUG ("Calculation data initialization for %s", inp_data->name);

  calc_data->temp_kin[0] = inp_data->temp_kin;
  calc_data->temp_bg[0] = inp_data->temp_bg;
  calc_data->col_dens[0] = inp_data->col_dens;
  calc_data->line_width[0] = inp_data->line_width;
  calc_data->geom[0] = inp_data->geom;

  RXI_STAT status = RXI_OK;
  struct rxi_db_molecule_enlev *mol_enl;
  status = rxi_db_molecule_enlev_malloc (&mol_enl, mol_info->numof_enlev);
  if (status != RXI_OK)
    goto error;
  status = rxi_db_read_molecule_enlev (inp_data->name, mol_enl);
  if (status != RXI_OK)
    {
      rxi_db_molecule_enlev_free (mol_enl);
      goto error;
    }

  DEBUG ("Molecule enlev parameters were read");

  struct rxi_db_molecule_radtr *mol_rt;
  status = rxi_db_molecule_radtr_malloc (&mol_rt, mol_info->numof_radtr);
  if (status != RXI_OK)
    {
      rxi_db_molecule_enlev_free (mol_enl);
      goto error;
    }
  status = rxi_db_read_molecule_radtr (inp_data->name, mol_rt);
  if (status != RXI_OK)
    {
      rxi_db_molecule_enlev_free (mol_enl);
      rxi_db_molecule_radtr_free (mol_rt);
      goto error;
    }

  DEBUG ("Molecule radtr parameters were read");

  struct rxi_db_molecule_coll_part **mol_cp = malloc (
      inp_data->n_coll_partners * sizeof (**mol_cp));

  for (int8_t i = 0; i < inp_data->n_coll_partners; ++i)
    {
      int8_t cp = cptonum (mol_info, inp_data->coll_part[i]);
      status = rxi_db_molecule_coll_part_malloc (&mol_cp[i],
          mol_info->numof_coll_trans[cp], mol_info->numof_coll_temps[cp]);
      if (status != RXI_OK)
        {
          rxi_db_molecule_enlev_free (mol_enl);
          rxi_db_molecule_radtr_free (mol_rt);
          free (*mol_cp);
          goto error;
        }
      status = rxi_db_read_molecule_coll_part (inp_data->name,
          inp_data->coll_part[i], mol_info->numof_coll_temps[cp], mol_cp[i]);
      if (status != RXI_OK)
        {
          rxi_db_molecule_enlev_free (mol_enl);
          rxi_db_molecule_radtr_free (mol_rt);
          free (*mol_cp);
          goto error;
        }

      DEBUG ("Molecule collision transfer parameters were read");
    }

  status = rxi_calc_data_fill (inp_data, mol_info, mol_enl, mol_rt, mol_cp,
                               calc_data);
  if (status != RXI_OK)
    {
      rxi_db_molecule_enlev_free (mol_enl);
      rxi_db_molecule_radtr_free (mol_rt);
      for (int8_t i = 0; i < inp_data->n_coll_partners; ++i)
        rxi_db_molecule_coll_part_free (mol_cp[i]);
      goto error;
    }

  rxi_calc_bgfield (calc_data, mol_rt, mol_info->numof_radtr);
  set_starting_conditions (calc_data, mol_info->numof_radtr);

  rxi_db_molecule_enlev_free (mol_enl);
  rxi_db_molecule_radtr_free (mol_rt);

  for (int8_t i = 0; i < inp_data->n_coll_partners; ++i)
    rxi_db_molecule_coll_part_free (mol_cp[i]);

  return status;

error:
  return status;
}

static double
interpolate_cp_rate (const double kin_temp, const double *temps,
                     const double *rates, const size_t n_temps)
{
  DEBUG ("Calculating cp rate for given temperature");
  double lcoef = 0;
  double ucoef = 0;
  double ltemp = 0;
  double utemp = 0;
  for (size_t i = 0; i < n_temps; ++i)
    {
      if (kin_temp > temps[i])
        continue;

      lcoef = rates[i];
      ltemp = temps[i];
      if (i == n_temps - 1)
        break;

      ucoef = rates[i + 1];
      utemp = temps[i + 1];
      break;
    }

  DEBUG ("lcoef: %.3e | ltemp: %.3e | ucoef: %.3e | utemp: %.3e", lcoef, ltemp,
         ucoef, utemp);
  if ((ltemp < 1e-30) && (utemp < 1e-30))
    {
      // Case when kinetic temperature is lower than minimum collision rate
      // temperature
      return rates[n_temps - 1];
    }
  else if (ltemp > utemp)
    {
      // Case when kinetic temperature is bigger than maximum collision rate
      // temperature
      return lcoef;
    }
  else
    {
      // Make linear interpolation to needed temperature
      return lcoef +
              (ucoef - lcoef) * (kin_temp - ltemp)
          / //------------------------------------
                        (utemp - ltemp);
    }
}

static double*
get_matrix_row (const gsl_matrix *matrix, const size_t n_row)
{
  DEBUG ("Get matrix row %zu with size %zu", n_row, matrix->size2);
  double *row = malloc (matrix->size2 * sizeof (*row));
  for (size_t i = 0; i < matrix->size2; ++i)
    {
      row[i] = gsl_matrix_get (matrix, n_row, i);
    }
  
  return row;
}

RXI_STAT
rxi_calc_data_fill (const struct rxi_input_data *inp_data,
                    const struct rxi_db_molecule_info *mol_info,
                    const struct rxi_db_molecule_enlev *mol_enlev,
                    const struct rxi_db_molecule_radtr *mol_radtr, 
                    struct rxi_db_molecule_coll_part **mol_cp,
                    struct rxi_calc_data *calc_data)
{
  DEBUG ("Setting terms and molecular weights");

  for (int i = 0; i < mol_info->numof_enlev; ++i)
    {
      gsl_vector_set (calc_data->term, mol_enlev->level[i] - 1,
                      mol_enlev->term[i]);
      gsl_vector_set (calc_data->weight, mol_enlev->level[i] - 1,
                      mol_enlev->weight[i]);
    }

  DEBUG ("Setting Einstein coefs and energies");

  for (int i = 0; i < mol_info->numof_radtr; ++i)
    {
      calc_data->up[i] = mol_radtr->up[i];
      calc_data->low[i] = mol_radtr->low[i];
      gsl_matrix_set (calc_data->einst, mol_radtr->up[i] - 1,
                      mol_radtr->low[i] - 1, mol_radtr->einst[i]);
      gsl_matrix_set (calc_data->energy, mol_radtr->up[i] - 1,
                      mol_radtr->low[i] - 1, mol_radtr->up_en[i]);
    }

  DEBUG ("Setting collision rates");

  for (int p = 0; p < inp_data->n_coll_partners; ++p)
    {
      // Get index number (from .info file) of entered collisional partner
      int8_t cp = cptonum (mol_info, inp_data->coll_part[p]);

      // Collisional coefficients will be written here
      gsl_matrix *coll_coef = gsl_matrix_calloc (
          mol_info->numof_enlev, mol_info->numof_enlev);

      // Set collisional coefficients for current partner
      double *temps_line = get_matrix_row (mol_info->coll_temps, cp);
      for (int i = 0; i < mol_info->numof_coll_trans[cp]; ++i)
        {
          double *rates_line = get_matrix_row (mol_cp[p]->coll_rates, i);
          double coef = interpolate_cp_rate (inp_data->temp_kin,
              temps_line, rates_line, mol_info->numof_coll_temps[cp]);
          gsl_matrix_set (coll_coef, mol_cp[p]->up[i] - 1,
                          mol_cp[p]->low[i] - 1, coef);

          free (rates_line);
        }
      free (temps_line);

      // And here coefficients become collisional rates (but not final ones)
      gsl_matrix_scale (coll_coef, inp_data->coll_part_dens[p]);
      gsl_matrix_add (calc_data->coll_rates, coll_coef);

      gsl_matrix_free (coll_coef);
    }

  // Cannot do this with common gsl matrix operations
  for (int i = 0; i < mol_info->numof_enlev; ++i)
    {
      for (int j = 0; j < mol_info->numof_enlev; ++j)
        {
          gsl_matrix_set (calc_data->rates, i, j, 1e-30);
          double ediff = gsl_vector_get (calc_data->term, i)
                         - gsl_vector_get (calc_data->term, j);
          if (ediff < 0)
            continue;

          double rate = rxi_calc_crate (gsl_vector_get (calc_data->weight, i),
              gsl_vector_get (calc_data->weight, j), ediff, inp_data->temp_kin,
              gsl_matrix_get (calc_data->coll_rates, i, j));

          gsl_matrix_set (calc_data->coll_rates, j, i, rate);
        }
    }

  // Calculate total collisional rates
  for (int i = 0; i < mol_info->numof_enlev; ++i)
    {
      gsl_vector *conv = gsl_vector_alloc (mol_info->numof_enlev);
      gsl_matrix_get_col (conv, calc_data->coll_rates, i);
      gsl_vector_add (calc_data->tot_rates, conv);
      gsl_vector_free (conv);
    }

  return RXI_OK;
}

RXI_STAT
rxi_calc_find_rates (struct rxi_calc_data *data, const int n_enlev,
                     const int n_radtr)
{
  unsigned int iter = 0;
  int thick_lines = 1;
  double stop_condition = 0;
  gsl_vector *prev_pop = gsl_vector_calloc (n_enlev);
  do 
    {
      if (iter != 0)
        thick_lines = refresh_starting_conditions (data, n_radtr);

      stop_condition = 0;
      // Correct rates for collisional rates and prepare new matrix
      // for calculations
      for (int i = 0; i < n_enlev; ++i)
        {
          const double ii = gsl_matrix_get (data->rates, i, i)
                            + gsl_vector_get (data->tot_rates, i);
          gsl_matrix_set (data->rates, i, i, ii);

          for (int j = 0; j < n_enlev; ++j)
            {
              if (i == j)
                continue;

              const double ij = gsl_matrix_get (data->rates, i, j)
                                - gsl_matrix_get (data->coll_rates, j, i);
              gsl_matrix_set (data->rates, i, j, ij);
            }
        }

      // Prepare for calculations
      gsl_vector *b = gsl_vector_alloc (n_enlev);
      gsl_vector_set_all (b, 1);
      gsl_matrix_set_row (data->rates, data->rates->size1 - 1, b);
      gsl_vector_set_all (b, 0);
      gsl_vector_set (b, b->size - 1, 1);
      gsl_vector *x = gsl_vector_alloc (n_enlev);

      for (int i = 0; i < n_enlev; ++i)
        {
          for (int j = 0; j < n_enlev; ++j)
            {
              printf ("%9.2e  ", gsl_matrix_get (data->rates, i, j));
            }
          printf ("\n");
        }

      gsl_permutation *p = gsl_permutation_alloc (n_enlev);
      int s;
      gsl_linalg_LU_decomp (data->rates, p, &s);
      gsl_linalg_LU_solve (data->rates, p, b, x);
      //gsl_linalg_HH_solve (data->rates, b, x);
/*
      for (int i = 0; i < n_enlev; ++i)
        {
          for (int j = 0; j < n_enlev; ++j)
            {
              printf ("%9.2e  ", gsl_matrix_get (calc_rates, i, j));
            }
          printf ("\n");
        }
*/
      gsl_vector_fprintf (stdout, x, "%.3e");
      double total_pop = 0;
      for (int i = 0; i < n_enlev; ++i)
        total_pop += gsl_vector_get (x, i);

      DEBUG ("total population %f", total_pop);

      gsl_vector_memcpy (prev_pop, data->pop);
      for (int i = 0; i < n_enlev; ++i)
        {
          const double new_pop_i = gsl_vector_get (x, i) / total_pop;
          gsl_vector_set (data->pop, i, fabs (new_pop_i));
        }

      if (iter == 0)
        gsl_vector_memcpy (prev_pop, data->pop);

      for (int i = 0; i < n_radtr; ++i)
        {
          const unsigned int u = data->up[i] - 1;
          const unsigned int l = data->low[i] - 1;

          const double new_excit_temp_i = RXI_FK
              * (gsl_vector_get (data->term, u)
                 - gsl_vector_get (data->term, l)) 
              / log (gsl_vector_get (data->pop, l)
                     * gsl_vector_get (data->weight, u)
                     / gsl_vector_get (data->pop, u)
                     / gsl_vector_get (data->weight, l));

          if (iter == 0)
            {
              gsl_matrix_set (data->excit_temp, u, l, new_excit_temp_i);
              stop_condition = 1;
            }
          else
            {
              gsl_matrix_set (data->excit_temp, u, l, 
                  0.5 * (new_excit_temp_i
                    + gsl_matrix_get(data->excit_temp, u, l)));
            }
          const double new_tau = rxi_calc_optical_depth (data->col_dens[0],
              data->line_width[0],
              gsl_vector_get (data->term, u) - gsl_vector_get (data->term, l),
              gsl_matrix_get (data->einst, u, l),
              gsl_vector_get(data->weight, u), gsl_vector_get (data->weight, l),
              gsl_vector_get(data->pop, u), gsl_vector_get (data->pop, l));

          if (new_tau > 0.01)
            stop_condition += fabs ((gsl_matrix_get (data->excit_temp, u, l)
                                     - new_excit_temp_i) / new_excit_temp_i);

          gsl_matrix_set (data->tau, u, l, new_tau);
          DEBUG ("%d -> %d: %.3e", u, l, new_tau);
        }

      for (int i = 0; i < n_enlev; ++i)
        {
          const double new_pop_i = 
                      0.3 * gsl_vector_get (data->pop, i) +
                      0.7 * gsl_vector_get (prev_pop, i);
          gsl_vector_set (data->pop, i, new_pop_i);
        }

      gsl_vector_free (b);
      gsl_vector_free (x);
      gsl_permutation_free (p);
      ++iter;
      DEBUG ("%d: Thick lines: %d | Stopping cond: %.3e", iter, thick_lines,
             stop_condition);
    } while ((thick_lines != 0 && stop_condition / thick_lines >= 1e-6) &&
              iter < 300);

  gsl_vector_free (prev_pop);

  return RXI_OK;
}

double
rxi_calc_crate (const double istat, const double jstat, const double ediff,
                const double kin_temp, const double crate)
{
  return
          istat * exp (- RXI_FK * ediff / kin_temp) * crate
      / //-------------------------------------------------
                          jstat;
}

double
rxi_calc_escape_prob (const double tau, const GEOMETRY geom)
{
  double beta = 0;
  const double tau_rad = tau / 2;

  if (geom == SPHERE)
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
  else if (geom == SLAB)
    {
      if (fabs (3 * tau) < 0.1)
        beta = 1 - 1.5 * (tau + tau * 2);
      else if (fabs (3 * tau) > 50)
        beta = 1 / (3 * tau);
      else 
        beta = (1 - exp (-3 * tau)) / (3 * tau);
    }
  else if (geom == LVG)
    {
      if (fabs (tau_rad) < 0.01)
        beta = 1;
      else if (fabs (tau_rad) < 7)
        beta = 2 * (1 - exp (-2.34 * tau_rad)) / (4.68 * tau_rad);
      else 
        beta = 2 / (tau_rad * 4 * sqrt (log (tau_rad / sqrt (M_PI))));
    }

  return beta;

}

double
rxi_calc_optical_depth (const double coldens, const double line_width,
    const double energy, const double einst, const double ustat,
    const double lstat, const double upop, const double lpop)
{
  /*DEBUG ("1e-5 * (%.3e * %.3e / %.3e - %.3e) * %.3e * %.3e\n\
                                                                 ----------------------------------------------------------------------------\n\
                                                                  %.3e * 1.0645 * 8 * pi * %.3e\n",
  lpop, ustat, lstat, upop, einst, coldens,
  gsl_pow_3 (energy), line_width);
*/
  return
          (1e-5 * (lpop * ustat / lstat - upop) * einst * coldens)
      / //--------------------------------------------------------
          (gsl_pow_3 (energy) * 1.0645 * 8.0 * M_PI * line_width);
}
