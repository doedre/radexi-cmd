/*	main.c
 *
 * The main file, where the 'main' function located. At first it controls
 * flags that were added to the program call, then it decides what to do.
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "rxi_common.h"
#include "core/dialogue.h"
#include "utils/options.h"
#include "utils/database.h"
#include "utils/debug.h"

RXI_STAT usage_dialogue (const struct rxi_options *opts);

RXI_STAT usage_print_version ();
RXI_STAT usage_print_help ();

int
main (int argc, char **argv)
{
  RXI_STAT return_value = RXI_OK;
  struct rxi_options opts;
  int index = rxi_set_options (&opts, argc, argv);

  switch (opts.usage_mode)
    {
    case UM_DIALOGUE:
      return_value = usage_dialogue (&opts);
      break;

    case UM_MOLECULAR_FILE_ADD:
      return_value = rxi_add_molecule (opts.molecule_name, argv[index]);
      break;

    case UM_MOLECULAR_FILE_DELETE:
      return_value = rxi_delete_molecule (opts.molecule_name);
      break;

    case UM_MOLECULAR_FILE_LIST:
      return_value = rxi_list_molecules ();
      break;

    case UM_VERSION:
      return_value = usage_print_version ();
      break;

    case UM_HELP:
      return_value = usage_print_help ();
      break;

    default:
      break;
    }

  printf ("Status: %u\n", return_value);
  return return_value;
}

RXI_STAT
usage_dialogue (const struct rxi_options *opts)
{
  struct rxi_input_data *inp_data = malloc (sizeof (*inp_data));
  RXI_STAT stat = rxi_dialog_input (inp_data, opts);
  CHECK ((stat == RXI_OK) && "Dialog errors");
  DEBUG ("Name: %s, cp: %u", inp_data->name, inp_data->coll_part[0]);

  struct rxi_db_molecule_info *mol_info;
  rxi_db_molecule_info_malloc (&mol_info);
  rxi_db_read_molecule_info (inp_data->name, mol_info);

  struct rxi_db_molecule_enlev *mol_enl;
  rxi_db_molecule_enlev_malloc (&mol_enl, mol_info->numof_enlev);
  stat = rxi_db_read_molecule_enlev (inp_data->name, mol_enl);

  DEBUG ("Molecule enlev parameters were read");
  for (int i = 0; i < mol_info->numof_enlev; ++i)
    {
      DEBUG ("%d : Level: %d | Weight: %.3f | Energy: %.3f | Qnum: %s", i,
             mol_enl->level[i], gsl_vector_get (mol_enl->weight, i),
             gsl_vector_get (mol_enl->energy, i), mol_enl->qnum[i]);
    }

  struct rxi_db_molecule_radtr *mol_rt;
  stat = rxi_db_molecule_radtr_malloc (&mol_rt, mol_info->numof_radtr);
  stat = rxi_db_read_molecule_radtr (inp_data->name, mol_rt);

  DEBUG ("Molecule radtr parameters were read");
  for (int i = 0; i < mol_info->numof_radtr; ++i)
    {
      DEBUG ("%d : Up: %d | Low: %d | Einst: %.3f | Freq: %.3f | Eu: %.3f", i,
             mol_rt->up[i], mol_rt->low[i], gsl_vector_get (mol_rt->einst, i),
             gsl_vector_get (mol_rt->freq, i),
             gsl_vector_get (mol_rt->up_en, i));
    }

  struct rxi_db_molecule_coll_part *mol_cp;
  stat = rxi_db_molecule_coll_part_malloc (&mol_cp,
      mol_info->coll_partners[0].numof_coll_trans,
      mol_info->coll_partners[0].numof_coll_temps);
  stat = rxi_db_read_molecule_coll_part (inp_data->name,
      inp_data->coll_part[0], mol_info->coll_partners[0].numof_coll_temps,
      mol_cp);

  DEBUG ("Molecule collision transfer parameters were read");
  for (int i = 0; i < mol_info->coll_partners[0].numof_coll_trans; ++i)
    {
      DEBUG ("%d : Up: %d | Low: %d | cp [3]: %.3e", i, mol_cp->up[i],
          mol_cp->low[i], gsl_matrix_get (mol_cp->coll_rates, i, 3));
    }

  free (inp_data);
  rxi_db_molecule_info_free (mol_info);
  rxi_db_molecule_enlev_free (mol_enl);
  rxi_db_molecule_radtr_free (mol_rt);
  rxi_db_molecule_coll_part_free (mol_cp);
  return stat;
}

RXI_STAT
usage_print_help ()
{
  printf ("HELP\n");
  return RXI_OK;
}

RXI_STAT
usage_print_version ()
{
  printf ("VERSION\n");
  return RXI_OK;
}
