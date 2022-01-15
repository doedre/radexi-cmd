/**
 * @file core/calculation.c
 */

#include <stdlib.h>
#include <stdio.h>

#include "core/calculation.h"

#include "rxi_common.h"
#include "utils/database.h"
#include "utils/debug.h"

RXI_STAT
rxi_calc_data_init (const struct rxi_input_data *inp_data)
{
  DEBUG ("Calculation data initialization for %s", inp_data->name);
  struct rxi_db_molecule_info *mol_info;
  rxi_db_molecule_info_malloc (&mol_info);
  rxi_db_read_molecule_info (inp_data->name, mol_info);

  struct rxi_db_molecule_enlev *mol_enl;
  rxi_db_molecule_enlev_malloc (&mol_enl, mol_info->numof_enlev);
  rxi_db_read_molecule_enlev (inp_data->name, mol_enl);

  DEBUG ("Molecule enlev parameters were read");
  for (int i = 0; i < mol_info->numof_enlev; ++i)
    {
      DEBUG ("%d : Level: %d | Weight: %.3f | Energy: %.3f | Qnum: %s", i,
             mol_enl->level[i], mol_enl->weight[i],mol_enl->term[i],
             mol_enl->qnum[i]);
    }

  struct rxi_db_molecule_radtr *mol_rt;
  rxi_db_molecule_radtr_malloc (&mol_rt, mol_info->numof_radtr);
  rxi_db_read_molecule_radtr (inp_data->name, mol_rt);

  DEBUG ("Molecule radtr parameters were read");
  for (int i = 0; i < mol_info->numof_radtr; ++i)
    {
      DEBUG ("%d : Up: %d | Low: %d | Einst: %.3e | Freq: %.3f | Eu: %.3f", i,
             mol_rt->up[i], mol_rt->low[i], mol_rt->einst[i],
             mol_rt->freq[i], mol_rt->up_en[i]);
    }

  struct rxi_db_molecule_coll_part *mol_cp;
  rxi_db_molecule_coll_part_malloc (&mol_cp,
      mol_info->coll_partners[0].numof_coll_trans,
      mol_info->coll_partners[0].numof_coll_temps);
  rxi_db_read_molecule_coll_part (inp_data->name,
      inp_data->coll_part[0], mol_info->coll_partners[0].numof_coll_temps,
      mol_cp);

  DEBUG ("Molecule collision transfer parameters were read");
  for (int i = 0; i < mol_info->coll_partners[0].numof_coll_trans; ++i)
    {
      DEBUG ("%d : Up: %d | Low: %d | cp [3]: %.3e", i, mol_cp->up[i],
          mol_cp->low[i], gsl_matrix_get (mol_cp->coll_rates, i, 3));
    }

  rxi_db_molecule_info_free (mol_info);
  rxi_db_molecule_enlev_free (mol_enl);
  rxi_db_molecule_radtr_free (mol_rt);
  rxi_db_molecule_coll_part_free (mol_cp);

  return RXI_OK;
}
