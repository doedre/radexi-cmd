/**
 * @file core/calculation.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "core/calculation.h"

#include "rxi_common.h"
#include "utils/database.h"
#include "utils/debug.h"

RXI_STAT
rxi_calc_data_init (const struct rxi_input_data *inp_data)
{
  DEBUG ("Calculation data initialization for %s", inp_data->name);

  RXI_STAT status = RXI_OK;
  struct rxi_db_molecule_info mol_info;
  status = rxi_db_molecule_info_malloc (&mol_info);
  if (status != RXI_OK)
    goto error;
  status = rxi_db_read_molecule_info (inp_data->name, &mol_info);
  if (status != RXI_OK)
    goto error;

  struct rxi_db_molecule_enlev *mol_enl;
  status = rxi_db_molecule_enlev_malloc (&mol_enl, mol_info.numof_enlev);
  if (status != RXI_OK)
    goto error;
  status = rxi_db_read_molecule_enlev (inp_data->name, mol_enl);
  if (status != RXI_OK)
    goto error;

  DEBUG ("Molecule information: Name: %s | Weight: %.3f | EL: %d | RT: %d\
| CP: %d", mol_info.name, mol_info.weight, mol_info.numof_enlev,
      mol_info.numof_radtr, mol_info.numof_coll_part);

  for (int i = 0; i < mol_info.numof_coll_part; ++i)
    {
      DEBUG ("Collision partner %d: Name: %d | Trans: %d | Temps: %d", i,
          mol_info.coll_part[i], mol_info.numof_coll_trans[i],
          mol_info.numof_coll_temps[i]);
    }
/*
  DEBUG ("Molecule enlev parameters were read");
  for (int i = 0; i < mol_info->numof_enlev; ++i)
    {
      DEBUG ("%d : Level: %d | Weight: %.3f | Energy: %.3f | Qnum: %s", i,
             mol_enl->level[i], mol_enl->weight[i], mol_enl->term[i],
             mol_enl->qnum[i]);
    }

  struct rxi_db_molecule_radtr *mol_rt;
  status = rxi_db_molecule_radtr_malloc (&mol_rt, mol_info->numof_radtr);
  if (status != RXI_OK)
    goto error;
  status = rxi_db_read_molecule_radtr (inp_data->name, mol_rt);
  if (status != RXI_OK)
    goto error;

  struct rxi_db_molecule_coll_part *mol_cp;
  status = rxi_db_molecule_coll_part_malloc (&mol_cp,
      mol_info->numof_coll_trans[0], mol_info->numof_coll_temps[0]);
  if (status != RXI_OK)
    goto error;
  status = rxi_db_read_molecule_coll_part (inp_data->name,
      inp_data->coll_part[0], mol_info->numof_coll_temps[0], mol_cp);
  if (status != RXI_OK)
    goto error;

  DEBUG ("Molecule collision transfer parameters were read");
  for (int i = 0; i < mol_info->numof_coll_trans[0]; ++i)
    {
      DEBUG ("%d : Up: %d | Low: %d | cp [3]: %.3e", i, mol_cp->up[i],
          mol_cp->low[i], gsl_matrix_get (mol_cp->coll_rates, i, 3));
    }

  struct rxi_calc_data *calc_data;
  status = rxi_calc_data_malloc (&calc_data, mol_info->numof_enlev,
                                 mol_info->numof_radtr);
  if (status != RXI_OK)
    goto error;
  status = rxi_calc_data_fill (inp_data, mol_info, mol_enl, mol_rt, &mol_cp,
                               calc_data);*/
/*
  rxi_db_molecule_enlev_free (mol_enl);
  rxi_db_molecule_radtr_free (mol_rt);
  rxi_db_molecule_coll_part_free (mol_cp);
*/

  rxi_db_molecule_info_free (mol_info);
  return RXI_OK;

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

      ucoef = rates[i];
      utemp = temps[i];
      break;
    }

  if (ltemp > utemp)
    {
      // Case when kinetic temperature is bigger than maximum collision rate
      // temperature
      return lcoef;
    }
  else if ((ltemp < 1e-30) && (utemp < 1e-30))
    {
      // Case when kinetic temperature is lower than minimum collision rate
      // temperature
      return rates[0];
    }
  else
    {
      // Make linear interpolation to needed temperature
      
      return 
                (lcoef + (ucoef - lcoef) * (kin_temp - ltemp))
        / //------------------------------------------------------
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
      DEBUG ("Matrix element %.3e", row[i]);
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
      DEBUG ("Set %d element", mol_enlev->level[i] - 1);
      gsl_vector_set (calc_data->term, mol_enlev->level[i] - 1,
                      mol_enlev->term[i]);
      gsl_vector_set (calc_data->weight, mol_enlev->level[i] - 1,
                      mol_enlev->weight[i]);
    }

  DEBUG ("Setting Einstein coefs and energies");

  for (int i = 0; i < mol_info->numof_radtr; ++i)
    {
      DEBUG ("Set %dx%d element", mol_radtr->up[i] - 1, mol_radtr->low[i] - 1);
      gsl_matrix_set (calc_data->einst, mol_radtr->up[i] - 1,
                      mol_radtr->low[i] - 1, mol_radtr->einst[i]);
      gsl_matrix_set (calc_data->energy, mol_radtr->up[i] - 1,
                      mol_radtr->low[i] - 1, mol_radtr->up_en[i]);
    }

  DEBUG ("Setting collision rates");

  for (int p = 0; p < inp_data->n_coll_partners; ++p)
    {
      int8_t cp = cptonum (mol_info, inp_data->coll_part[p]);
      for (int i = 0; i < mol_info->numof_coll_trans[cp]; ++i)
        {
          double *rates_line = get_matrix_row (mol_cp[0]->coll_rates, i);
          double *temps_line = get_matrix_row (mol_info->coll_temps, i);

          double coef = interpolate_cp_rate (inp_data->temp_kin,
              temps_line, rates_line, mol_info->numof_coll_temps[cp]);
          DEBUG ("interpolated %.3e coefficient", coef);

          free (rates_line);
          free (temps_line);
        }
    }

  return RXI_OK;
}
