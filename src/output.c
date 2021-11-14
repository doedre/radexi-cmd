/*	output.c
 *
 * Writes an output for the user.
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
#include "linenoise.h"

#include <sys/stat.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* Checks if specified result file already exists.
 * Returns:
 * - 0 if specified file doesn't exists and can be written
 * - 1 if specified file exists and will be rewritted
 * - 2 if directory was given instead
 * - -1 if specified file exists and wont be rewritten
 * - -2 if specified file doesn't exists and can't be written
 */
int
check_result_path (const char *result_path)
{
  printf ("--stat check for %s...\n", result_path);
  struct stat sb;
  stat(result_path, &sb);

  if (S_ISREG (sb.st_mode) && !rxi_opts.force_fs)
    {
      printf (BMAG  "  ## " reset
              WHT   "Specified file \"%s\" already exists. Rewrite it?" reset
              BWHT  " (Y/n)" reset "\n"
              , result_path);

      char *choice;
      choice = (char *) malloc (1);
      while ((choice = linenoise ("  >> ")) != NULL)
        {
          if (!strcmp (choice, "y") || !strcmp (choice, "Y"))
            {
              free (choice);
              return 1;
            }
          else 
            {
              free (choice);
              return -1;
            }
        }
    }
  else if (S_ISDIR (sb.st_mode))
    {
      return 2;
    }

  printf ("--opening new file...\n");
  FILE *check;
  check = fopen (result_path, "w");
  if (check)
    {
      fclose (check);
      return 0;
    }
  else 
    {
      return -2;
    }
}

/* Writes calculation results into file.
 * Returns:
 * - 0 on success
 * - -1 on file creation error
 * - -2 on rxi_data problems
 */
int
write_results(const struct rxi_data *rxi, const char *result_path)
{
  struct stat sb;
  stat(result_path, &sb);
  char *path;
  path = malloc (100);
  strcpy (path, result_path);

  if (S_ISDIR (sb.st_mode))
    {
      strcat (path, "/");      
      strcat (path, "radexi_output.txt");
    }

  FILE *result_file;
  result_file = fopen (path, "w");
  if (!result_file)
    {
      return -1;
    }
  
  fprintf (result_file,     "* Radexi version             : "           PROGRAM_VERSION "\n");
  fprintf (result_file,     "* Geometry                   : %d\n",      rxi->mc.geom);
  fprintf (result_file,     "* Molecule                   : %s\n",      rxi->mi.name);
  fprintf (result_file,     "* Kinetic temperature    [K] : %.3f\n",    rxi->mc.Tkin);
  fprintf (result_file,     "* Background temperature [K] : %.3f\n",    rxi->mc.Tbg);
  fprintf (result_file,     "* Column density      [cm-2] : %.3e\n",    rxi->mc.coldens);
  fprintf (result_file,     "* Line width          [km/s] : %.3f\n",    rxi->mc.line_width);
  for (unsigned int i = 0; i < rxi->mc.numof_colpart; i++)
    {
      char cp_name[10];
      conv_int_to_name (rxi->mc.cp[i].name, cp_name);
      fprintf (result_file, "* Density of %9s[cm-3] : %.3e\n", cp_name, rxi->mc.cp[i].dens);
    }

  fprintf (result_file,     "*    LINE       E_UP          FREQ         WAVEL        T_EX        TAU        \
      T_R        POP        POP       FLUX       FLUX\n");
  fprintf (result_file,     "*                [K]         [GHz]          [nm]         [K]                   \
      [K]         UP        LOW   [K km/s]  [erg cm-2 s-1]\n");
  char line_format[120];
  strcpy (line_format, "%3d -> %3d  %8.1f  %12.4f  %12.4f  %10.3f  %.3e  %.3e  %.3e  %.3e  %.3e  %.3e\n");
  for (unsigned int i = 0; i < rxi->n_rt; i++)
    {
      double xt = pow (rxi->rad_transfer[i].xnu, 3);
      fprintf (result_file, line_format, rxi->rad_transfer[i].ulev, rxi->rad_transfer[i].llev, 
                                         rxi->rad_transfer[i].enup, 
                                         rxi->rad_transfer[i].spfreq, 
                                         sol / rxi->rad_transfer[i].spfreq / 1e5, 
                                         gsl_vector_get (rxi->res.Tex, i), 
                                         rxi->res.tau[i], 
                                         rxi->res.Tant[i],
                                         gsl_vector_get (rxi->res.pop, rxi->rad_transfer[i].ulev),
                                         gsl_vector_get (rxi->res.pop, rxi->rad_transfer[i].llev), 
                                         1.0645 * rxi->mc.line_width * rxi->res.Tant[i],
                                         1.0645 * 8 * M_PI * kB * rxi->mc.line_width * rxi->res.Tant[i] * xt);
    }

  fclose (result_file);
  return 0;
}
