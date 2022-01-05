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

#include "rxi_common.h"
#include "utils/options.h"

RXI_STAT usage_dialogue ();

int
main (int argc, char **argv)
{
  RXI_STAT return_value = RXI_OK;
  struct rxi_options opts;
  int index = rxi_set_options (&opts, argc, argv);

  switch (opts.usage_mode)
    {
    case UM_DIALOGUE:
      return_value = usage_dialogue ();
      break;


    default:
      break;
    }

	exit (return_value);
}

RXI_STAT usage_dialogue ()
{

}
