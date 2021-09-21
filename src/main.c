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

int
main (int argc, char ** argv)
{
  struct rx_options rx_opts;
  struct radexi_data rxi;
  float sf, ef;
  size_t s;
  int pathindex = set_rx_options (&rx_opts, argc, argv);

  if (rx_opts.usage_mode == UM_DIALOGUE)
    {
      start_dialogue (&sf, &ef, &rxi, &s, &rx_opts);
      reading_data (&rxi);
      printf ("name: %s, weight: %f, numofcolpart: %u", rxi.mi.name, 
                                                        rxi.mi.weight,
                                                        rxi.mi.numof_colpart);
    }
  else if (rx_opts.usage_mode == UM_FILE)
    {
    }
  else
    {
      operate_molecular_files (argv[pathindex], &rx_opts);
    }

	exit (EXIT_SUCCESS);
}
