/*	dialogue.c
 *
 * Contains realisations for functions, needed in dialogue w/ the user. Primary
 * goal is to collect all the needed info for calculations w/o troubles.
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

#include <string.h>
#include <math.h>

#include "radexi.h"
#include "linenoise.h"


/* States of the dialogue w/ the user. All possible errors are collected here 
 * and used in 'dialogue_usage()' function. */
enum 
{
  PARAMETERS_NORMAL = 1,
  MOLECULE_FAIL,
  FREQ_RANGE_FAIL,
  TKIN_FAIL,
  TBG_FAIL,
  COLDENS_FAIL,
  FWHM_FAIL,
  COLLISION_FAIL,
  GEOMETRY_FAIL,
  EXIT_REQUEST,
};

/* Those two are used to lower the amount of input functions used by this 
 * program. Need to rethink their usage in future.  */
typedef bool (allowed_value)(const char *str, 
                             double *par, 
                             const struct rxi_options *rx_opts);

/* This function takes the state argument, which is one from enum above, and 
 * then it either shows an error or exits (tries to do it good). */
static void
dialogue_usage (int state)
{
  if (state == MOLECULE_FAIL)
    { 
      printf ("\x1B[0;37;40m\x1B[1;31;40m  ## \x1B[0;37;40m");
      printf ("No such molecule was added\x1B[3;37;40m\n");
    }
  else if (state == FREQ_RANGE_FAIL)
    {
      printf ("\x1B[0;37;40m\x1B[1;31;40m  ## \x1B[0;37;40m");
      printf ("Wrong frequency range. Enter two numbers\x1B[3;37;40m\n");
    }
  else if (state == TKIN_FAIL)
    {
      printf ("\x1B[0;37;40m\x1B[1;31;40m  ## \x1B[0;37;40m");
      printf ("Wrong input. Check for temperature range\x1B[3;37;40m\n");
    }
  else if (state == TBG_FAIL)
    {
      printf ("\x1B[0;37;40m\x1B[1;31;40m  ## \x1B[0;37;40m");
      printf ("Wrong input. Check for bg temperature range\x1B[3;37;40m\n");
    }
  else if (state == COLDENS_FAIL)
    {
      printf ("\x1B[0;37;40m\x1B[1;31;40m  ## \x1B[0;37;40m");
      printf ("Wrong input. Check for coldens range\x1B[3;37;40m\n");
    }
  else if (state == FWHM_FAIL)
    {
      printf ("\x1B[0;37;40m\x1B[1;31;40m  ## \x1B[0;37;40m");
      printf ("Wrong input. Check for FWHM range\x1B[3;37;40m\n");
    }
  else if (state == COLLISION_FAIL)
    {
      printf ("\x1B[0;37;40m\x1B[1;31;40m  ## \x1B[0;37;40m");
      printf ("Wrong input. Collision\x1B[3;37;40m\n");
    }
  else if (state == COLLISION_FAIL)
    {
      printf ("\x1B[0;37;40m\x1B[1;31;40m  ## \x1B[0;37;40m");
      printf ("Wrong input. Geometry\x1B[3;37;40m\n");
    }
  else if (state == EXIT_REQUEST)
    {
      printf ("\x1B[0;37;40m");
      exit (EXIT_SUCCESS);
    }
}

/* Add all names and exceptions if no molecule was added  */
static void
molecules_completion (const char *buf, linenoiseCompletions *lc)
{
  if (buf[0] == 'c')
    {
      linenoiseAddCompletion (lc, "ch3oh_a");
      linenoiseAddCompletion (lc, "ch3oh_e");
    }
  else if (buf[0] == 'l')
    {
      linenoiseAddCompletion (lc, "list");
    }
}

/* Add all names and exceptions if no molecule was added  */
static void
geometry_completion (const char *buf, linenoiseCompletions *lc)
{
  if (buf[0] == 'l')
    {
      linenoiseAddCompletion (lc, "lvg");
    }
  else if (buf[0] == 's')
    {
      linenoiseAddCompletion (lc, "sphere");
      linenoiseAddCompletion (lc, "slab");
    }
}

/* Helps w/ kinetic temperature range.  */
static char*
hints_temp (const char *buf, int *color, int *bold)
{
  float f = atof (buf);
  if (f > 10000)
    {
      *color = 31;
      *bold = 3;
      return " > 10000";
    }
  else if (f < 0.1)
    { 
      *color = 31;
      *bold = 3;
      return " < 0.1";
    }

  return NULL;
}

/* Searches if the molecule, entered by the user, is added to the
 * radexi database  */
static bool
isAllowedMolecule (const char * str)
{
  bool res = false;
  if (!strcmp (str, "ch3oh_e") || !strcmp (str, "hco"))
    res = true;
  return res;
}

/* Used to fill molecule name in MC_parameters structure. */
static void
enter_molecule (char *name)
{
  char *line;
  line = (char *) malloc (RXI_MOLECULE_MAX_SIZE);
  linenoiseSetCompletionCallback (molecules_completion);
  linenoiseHistoryLoad (".molecules");

  while ((line = linenoise ("  >> ")) != NULL)
    {
      if (isAllowedMolecule (line))
        {
          strncpy (name, line, RXI_MOLECULE_MAX_SIZE);
          linenoiseHistoryAdd (line);
          linenoiseHistorySave (".molecules");
          break;
        }
      else if (!strcmp (line, "list"))
        {
          printf (" TODO: add list option ");
        }
      else if (!strcmp (line, "quit"))
        {
          dialogue_usage (EXIT_REQUEST);
        }
      else 
        {
          dialogue_usage (MOLECULE_FAIL);
        } 
          
      free (line);
    }
}

/* Parses user's input and fills starting & ending frequencies. Returns 'true'
 * if there was 2 numbers that were parsed and casted to floats w/o an errors.
 * Numbers are written to 'sf' and 'ef'.  */
static bool
isAllowedFrequencies (const char *str, float *sf, float *ef)
{
  bool res = false;
  char *storage;
  storage = (char *) malloc (15);
  strcpy (storage, str);

  /* Parse line and then cast its containments  */
  int numof_tokens = 0;
  *sf = 0;
  *ef = 0;
  for (char *token = strtok (storage, " "); token; token = strtok (NULL, " ")) 
    {
      if (numof_tokens == 0)
        {
          *sf = atof (token);
        }
      else if (numof_tokens == 1)
        {
          *ef = atof (token);
          /* If user prints '+10', then it's added to the starting frequency */
          if ((!strncmp (token, "+", 1)) && (*ef != 0))
            *ef += *sf;
        }
      else break;
      numof_tokens++;
    }

  if (*sf > 0 && *ef > 0)
      res = true;

  free (storage);
  return res;
}

/* Used to fill starting and ending frequencies. */
static void
enter_frequencies (float *sf, float *ef)
{
  size_t size = 15;
  char *line;
  line = (char *) malloc (size);
  linenoiseHistoryLoad (".frequencies");

  while ((line = linenoise ("  >> ")) != NULL)
    {
      if (isAllowedFrequencies (line, sf, ef))
        {
          linenoiseHistoryAdd (line);
          linenoiseHistorySave (".frequencies");
          break;
        }
      else if (!strcmp (line, "quit"))
        {
          dialogue_usage (EXIT_REQUEST);
        }
      else 
        {
          dialogue_usage (FREQ_RANGE_FAIL);
        } 
      free (line);
    }
}
   

/* Casts str to float and checks if it's of the correct range. */
static bool 
isAllowedTemp (const char *str, 
               double *val, 
               const struct rxi_options *rx_opts)
{
  /* Supress unused parameter warning */
  (void) rx_opts;

  bool res = false;
  double f = strtod (str, NULL);

  if ((f > 0.1) && (f < 10000))
    {
      *val = f;
      res = true;
    }

  return res;
}


/* Casts str to float and checks if it's of the correct column density range. 
 * Also works w/ -L flag. */
static bool
isAllowedColdens (const char *str, 
                  double *val, 
                  const struct rxi_options *rx_opts)
{
  bool res = false;
  double f = strtod (str, NULL);

  /* Checks for -l option */
  if (rx_opts->dens_log_scale)
    {
      if ((f > 5) && (f < 25))
        {
          *val = pow (10, f);
          res = true;
        }
    }
  else
    {
      if ((f > 1e5) && (f < 1e25))
        {
          *val = f;
          res = true;
        }
    }

  return res;
}

/* Casts str to float and checks if it's of the correct line widths range. 
 * Also works w/ -H flag. */
static bool
isAllowedWidth (const char *str, 
                double *val, 
                const struct rxi_options *rx_opts)
{
  bool res = false;
  double f = strtod (str, NULL);

  if (rx_opts->hz_width)
    {
      if ((f > 1e-4) && (f < 1e4))
        {
          *val = f;
          res = true;
        }
    }
  else
    {
      if ((f > 1e-3) && (f < 1e3))
        {
          *val = f;
          res = true;
        }
    }

  return res;
}

/* Used for typical float parameter entering. */
static void
enter_double_parameter (double *par, 
                        int fail_state,
                        const struct rxi_options *rx_opts,
                        const char *history_file, 
                        linenoiseCompletionCallback *completion,
                        linenoiseHintsCallback *hints,
                        allowed_value isAllowed)
{
  size_t size = 20;
  char *line;
  line = (char *) malloc (size);
  linenoiseHistoryLoad (history_file);
  linenoiseSetCompletionCallback (completion);
  linenoiseSetHintsCallback (hints);

  while ((line = linenoise ("  >> ")) != NULL)
    {
      if (isAllowed (line, par, rx_opts))
        {
          linenoiseHistoryAdd (line);
          linenoiseHistorySave (history_file);
          break;
        }
      else if (!strcmp (line, "quit"))
        {
          dialogue_usage (EXIT_REQUEST);
        }
      else 
        {
          dialogue_usage (fail_state);
        } 
      free (line);
    }
}

/* Converts collision partner's name into a number (or ColPartner enum).  */
static enum ColPart
conv_name_to_int (const char *str)
{
  enum ColPart res = He;
  if (!strncmp (str, "H2", 2) || !strncmp (str, "h2", 2))
    res = H2;
  else if (!strncmp (str, "PARA_H2", 7) || !strncmp (str, "para_h2", 7))
    res = PARA_H2;
  else if (!strncmp (str, "ORTHO_H2", 8) || !strncmp (str, "ortho_h2", 8))
    res = ORTHO_H2;
  else if (!strncmp (str, "ELECTRONS", 9) || !strncmp (str, "electrons", 9))
    res = ELECTRONS;
  else if (!strncmp (str, "HI", 2) || !strncmp (str, "hi", 2))
    res = HI;
  else if (!strncmp (str, "He", 2) || !strncmp (str, "he", 2))
    res = He;
  else if (!strncmp (str, "HII", 3) || !strncmp (str, "hii", 3))
    res = HII;
  else 
    res = NO_MOLECULE;

  return res;
}

/* Parses a line to get collision partner's name and it's volume density. */
static bool
isAllowedCps (const char *str, 
              struct rxi_input *inp, 
              const struct rxi_options *rx_opts)
{
  bool res = true;
  char *storage;
  storage = (char *) malloc (40);
  strcpy (storage, str);

  /* Parse line and then cast its containments by different molecules  */
  int numof_cps = 0;
  for (char *token = storage; token; token = strchr (token, ','))
    {
      /* Parse molecule by its name and density */
      if (numof_cps != 0)
        token++;
      int numof_t = 0;
      char *storage2;
      storage2 = (char *) malloc (40);
      strcpy (storage2, token);
      for (char *t = strtok (storage2, " "); t; t = strtok (NULL, " "))
        {
          if (numof_t == 0)
              inp->cps[numof_cps].name = conv_name_to_int (t);
          else if (numof_t == 1)
            {
              float d = atof (t);
              if (rx_opts->dens_log_scale)
                {
                  if (d >= 0 && d <= 14)
                    inp->cps[numof_cps].dens = powf (10, d);
                  else
                    res = false;
                }
              else
                {
                  if (d >= 1 && d <= 1e14)
                    inp->cps[numof_cps].dens = d;
                  else 
                    res = false;
                }
            }
          else 
              break;
          numof_t++;
        }
      numof_cps++;
      free (storage2);
    }

  inp->numof_colpart = numof_cps;
  free (storage);
  return res;
}

/* Equal to the enter_float(), but used for collision partner's retrieval.  */
static void
enter_collision_partners (struct rxi_input *inp,
                          int fail_state,
                          const struct rxi_options *rx_opts)
{
  size_t size = 30;
  char *line;
  line = (char *) malloc (size);
  linenoiseHistoryLoad (".colpartners");

  while ((line = linenoise ("  >> ")) != NULL)
    {
      if (isAllowedCps (line, inp, rx_opts))
        {
          linenoiseHistoryAdd (line);
          linenoiseHistorySave (".colpartners");
          break;
        }
      else if (!strcmp (line, "quit"))
        {
          dialogue_usage (EXIT_REQUEST);
        }
      else 
        {
          dialogue_usage (fail_state);
        } 
      free (line);
    }
}

bool isAllowedGeometry (const char *line, struct rxi_input *inp)
{
  bool result = false;
  if (!strncmp (line, "lvg", 3))
    {
      inp->g = LVG;
      result = true;
    }
  else if (!strncmp (line, "slab", 4))
    {
      inp->g = SLAB;
      result = true;
    }
  else if (!strncmp (line, "sphere", 6))
    {
      inp->g = SPHERE;
      result = true;
    }
  
  return result;
}

void 
enter_geometry (struct rxi_input *inp, 
                int fail_state)
{
  size_t size = 10;
  char *line;
  line = (char *) malloc (size);
  linenoiseSetCompletionCallback (geometry_completion);
  linenoiseHistoryLoad (".geometry");

  while ((line = linenoise ("  >> ")) != NULL)
    {
      if (isAllowedGeometry (line, inp))
        {
          linenoiseHistoryAdd (line);
          linenoiseHistorySave (".colpartners");
          break;
        }
      else if (!strcmp (line, "quit"))
        {
          dialogue_usage (EXIT_REQUEST);
        }
      else 
        {
          dialogue_usage (fail_state);
        } 
      free (line);
    }
}

/* This function is called by the main() for parameters entering. */
void
start_dialogue (float *sfreq, 
                float *efreq, 
                struct rxi_input *inp,
                const struct rxi_options *opts) 
{
  if (!opts->quite_start)
    {
      printf ("starting message here\n");
    }

  linenoiseHistorySetMaxLen (10);

  printf ("\x1B[0;37;40m\x1B[1;35;40m  ## \x1B[0;37;40m");
  printf ("Enter molecule's name ");
  printf ("(\x1B[2;37;40m'\x1B[4;37;40mlist\x1B[0;37;40m\x1B[2;37;40m'");
  printf (" to see the variants\x1B[0;37;40m)\x1B[3;37;40m\n");
 
  enter_molecule (inp->name);
      
  printf ("\x1B[0;37;40m\x1B[1;35;40m  ## \x1B[0;37;40m");
  printf ("Enter starting & ending frequencies ");
  printf ("\x1B[1;37;40m[GHz]\x1B[0;37;40m\x1B[3;37;40m\n");

  enter_frequencies (sfreq, efreq);

  printf ("\x1B[0;37;40m\x1B[1;35;40m  ## \x1B[0;37;40m");
  printf ("Enter kinetic temperature ");
  printf ("\x1B[1;37;40m[K]\x1B[0;37;40m\x1B[3;37;40m\n");

  enter_double_parameter (&inp->Tkin, 
                          TKIN_FAIL, 
                          opts, 
                          ".Tkin", 
                          NULL, 
                          hints_temp, 
                          isAllowedTemp);

  printf ("\x1B[0;37;40m\x1B[1;35;40m  ## \x1B[0;37;40m");
  printf ("Enter background temperature ");
  printf ("\x1B[1;37;40m[K]\x1B[0;37;40m\x1B[3;37;40m\n");

  enter_double_parameter (&inp->Tbg, 
                          TBG_FAIL, 
                          opts, 
                          ".Tbg", 
                          NULL, 
                          hints_temp, 
                          isAllowedTemp);

  printf ("\x1B[0;37;40m\x1B[1;35;40m  ## \x1B[0;37;40m");
  printf ("Enter column density for the molecule ");
  printf ("\x1B[1;37;40m[cm-2]\x1B[0;37;40m\n");

  enter_double_parameter (&inp->coldens, 
                          COLDENS_FAIL, 
                          opts, 
                          ".coldens",
                          NULL, 
                          NULL, 
                          isAllowedColdens);

  printf ("\x1B[0;37;40m\x1B[1;35;40m  ## \x1B[0;37;40m");
  printf ("Enter FWHM lines width ");
  if (opts->hz_width) 
    printf ("\x1B[1;37;40m[Hz]\x1B[0;37;40m\n");
  else 
    printf ("\x1B[1;37;40m[km s-1]\x1B[0;37;40m\n");

  enter_double_parameter (&inp->fwhm, 
                          FWHM_FAIL, 
                          opts, 
                          ".fwhm",
                          NULL, 
                          NULL, 
                          isAllowedWidth);

  printf ("\x1B[0;37;40m\x1B[1;35;40m  ## \x1B[0;37;40m");
  printf ("Enter collision partners and their densities ");
  printf ("\x1B[1;37;40m[cm-3]\x1B[0;37;40m\n");

  enter_collision_partners (inp, 
                            COLLISION_FAIL, 
                            opts);

  printf ("\x1B[0;37;40m\x1B[1;35;40m  ## \x1B[0;37;40m");
  printf ("Enter type of the cloud's geometry ");
  printf ("\x1B[1;37;40m[cm-3]\x1B[0;37;40m\n");

  enter_geometry (inp, GEOMETRY_FAIL);

}
