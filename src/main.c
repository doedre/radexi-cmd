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
#include "core/calculation.h"
#include "core/output.h"
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

  struct rxi_db_molecule_info *info;
  stat = rxi_db_molecule_info_malloc (&info);
  CHECK ((stat == RXI_OK) && "Info memory allocation error");
  stat = rxi_db_read_molecule_info (inp_data->name, info);
  CHECK ((stat == RXI_OK) && "Info file error");

  struct rxi_calc_data *calc_data;
  stat = rxi_calc_data_malloc (&calc_data, info->numof_enlev,
                               info->numof_radtr);
  CHECK ((stat == RXI_OK) && "Calculation data memory allocation error");
  stat = rxi_calc_data_init (calc_data, inp_data, info);
  CHECK ((stat == RXI_OK) && "Calculation data initialization error");
  stat = rxi_calc_find_rates (calc_data, info->numof_enlev, info->numof_radtr);
  CHECK ((stat == RXI_OK) && "Error in rates calculation");

  stat = rxi_out_result (calc_data, opts);
  CHECK ((stat == RXI_OK) && "Error in result printing");

  free (inp_data);
  rxi_db_molecule_info_free (info);
  rxi_calc_data_free (calc_data);
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
