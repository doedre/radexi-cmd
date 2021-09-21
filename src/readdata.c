/*	readdata.c
 *
 * Used to control reading of database for the specified molecule. 
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
#include <string.h>

int 
reading_data (struct radexi_data *rxi)
{
  char *mi_path;
  mi_path = (char *) malloc (2*RXI_MOLECULE_MAX_SIZE + 12);
  sprintf (mi_path, "data/%s/%s.info", rxi->mi.name, rxi->mi.name);
  FILE *mi = fopen (mi_path, "r");

  char *line;
  size_t n = 200;
  line = (char *) malloc (n);
  enum ColPartner col_partner = NO_MOLECULE;
  while (fgets (line, n, mi))
    {
      char *p_ddots = strrchr (line, ':');
      char *parameter = p_ddots ? p_ddots + 1 : line;
      if (strstr (line, "Molecular weight:"))
        rxi->mi.weight = strtof (parameter, NULL);
      else if (strstr (line, "Number of energy levels:"))
        rxi->mi.numof_enlev = strtof (parameter, NULL);
      else if (strstr (line, "Number of radiative transitions:"))
        rxi->mi.numof_radtr = strtof (parameter, NULL);
      else if (strstr (line, "Number of collision partners:"))
        rxi->mi.numof_colpart = strtof (parameter, NULL);
      else if (strstr (line, "Collision partner:"))
        col_partner = strtol (parameter, NULL, 10);
      else if (strstr (line, "Number of collisional transitions:"))
        {
          for (int i = 0; i < RXI_MAX_COLL_PARTNERS; i++)
            {
              if (rxi->mc_par.cps[i].name != col_partner)
                continue;
              rxi->mc_par.cps[i].numof_coltr = strtof (parameter, NULL);
            }
        }
      else if (strstr (line, "Collisional temperatures:"))
        {
          for (int i = 0; i < RXI_MAX_COLL_PARTNERS; i++)
            {
              if (rxi->mc_par.cps[i].name != col_partner)
                continue;
              char *end;
              int j = 0;
              for (float f = strtof (parameter, &end); parameter != end;
                                             f = strtof (parameter, &end))
                rxi->mc_par.cps[i].temps[j++] = f;
            }
        }
      else continue;
    }
  free (line);

  fclose (mi);
  return 0;
}
