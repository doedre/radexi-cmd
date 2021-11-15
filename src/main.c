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

#include "radexi.h"

/* Physical constants */
const double sol  = 2.99792458e10;    /* speed of light       [cm s-1]      */
const double hP   = 6.6260963e-27;    /* Planck's constant    [erg s]       */
const double kB   = 1.3806505e-16;    /* Boltzman's constant  [erg K-1]     */

const double fk   = hP * sol / kB;

struct rxi_options rxi_opts;
char radexi_path[PATH_MAX];

int
main (int argc, char **argv)
{
  struct rxi_input inp;
  float sf, ef;

  int pathindex = set_rxi_options (&rxi_opts, argc, argv);
  // If can't write default radexi path
  if (pathindex == -2)
    exit (EXIT_FAILURE);

  if (rxi_opts.usage_mode == UM_DIALOGUE)
    {
      if (rxi_opts.user_defined_out_file_path)
        {
          int check_path = check_result_path (rxi_opts.result_path);
          if (check_path < 0)
            {
              printf (BMAG  "  ## " reset
                      WHT   "Wrong path for result file" reset "\n");
              exit (EXIT_SUCCESS);
            }
        }
      start_dialogue (&sf, &ef, &inp);

      read_info_file (&inp);

      struct rxi_data *rxi = rxi_data_calloc (&inp);

      read_data (rxi);

      calculate_bg_field (rxi);

      main_calculations (rxi);

      calculate_results (sf, ef, rxi);

      write_results (sf, ef, rxi, rxi_opts.result_path);

      rxi_data_free (rxi);
    }
  else if (rxi_opts.usage_mode == UM_FILE)
    {
    }
  else
    {
      operate_molecular_files (argv[pathindex]);
    }

	exit (EXIT_SUCCESS);
}
