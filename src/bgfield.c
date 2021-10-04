/*	bgfield.c
 *
 * Calculates background field in three possible variants.
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

int
calculate_bg_field (struct radexi_data *rxi)
{
  for (unsigned int i = 0; i < rxi->mi.numof_radtr; i++)
    {
      rxi->bg.intens[i] = 2 * hP * sol * pow (rxi->mi.radtr[i].xnu, 3)    / \
            (exp (hP * sol / kB * rxi->mi.radtr[i].xnu / rxi->mc_par.Tbg) - 1);
      printf ("bgField: %.3e\n", rxi->bg.intens[i]);
    }

  return 0;
}

