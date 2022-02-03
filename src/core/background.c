/**
 * @file core/background.c
 */

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_const_cgsm.h>

#include "core/background.h"

#include "rxi_common.h"
#include "utils/debug.h"

void
rxi_calc_bgfield (struct rxi_calc_data *data,
    const struct rxi_db_molecule_radtr *mol_radtr, const int n_radtr)
{
  DEBUG ("Calculating background field intensity");

  for (int i = 0; i < n_radtr; ++i)
    {
      const unsigned int u = mol_radtr->up[i] - 1;
      const unsigned int l = mol_radtr->low[i] - 1;

      const double energy = gsl_vector_get (data->term, u)
                            - gsl_vector_get (data->term, l);

      const double intens =
              (2 * RXI_HP * RXI_SOL * gsl_pow_3 (energy))
          / //-------------------------------------------
              (exp (RXI_FK * energy / data->temp_bg[0]) - 1);

      gsl_matrix_set (data->bgfield, u, l, intens);
  }
}

