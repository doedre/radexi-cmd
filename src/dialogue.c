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
  EXIT_REQUEST,
};

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
  if (!strcmp (str, "ch3oh_a") || !strcmp (str, "ch3oh_e"))
    res = true;
  return res;
}


/* Used to fill molecule name in MC_parameters structure. */
static void
enter_molecule (struct MC_parameters *mc_pars)
{
  size_t size = 10;
  char *line;
  line = (char *) malloc (size);
  linenoiseSetCompletionCallback (molecules_completion);
  linenoiseHistoryLoad (".molecules");

  while ((line = linenoise ("  >> ")) != NULL)
    {
      if (isAllowedMolecule (line))
        {
          strncpy (mc_pars->molecule, line, 10);
          linenoiseHistoryAdd (line);
          linenoiseHistorySave (".molecules");
          break;
        }
      else if (!strcmp (line, "list"))
        {
          printf (" ch3oh_a and ch3oh_e ");
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
isAllowedTemp (const char *str, float *val, const struct rx_options *rx_opts)
{
  /* Supress unused parameter warning */
  (void) rx_opts;

  bool res = false;
  float f = atof (str);

  if ((f > 0.1) && (f < 10000))
    {
      *val = f;
      res = true;
    }

  return res;
}


/* Casts str to float and checks if it's of the correct range. */
static bool
isAllowedColdens (const char *str, float *val, 
                                              const struct rx_options *rx_opts)
{
  bool res = false;
  float f = atof (str);

  /* Checks for -l option */
  if (rx_opts->dens_log_scale)
    {
      if ((f > 5) && (f < 25))
        {
          *val = powf (10, f);
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

static bool
isAllowedWidth (const char *str, float *val, const struct rx_options *rx_opts)
{
  bool res = false;
  float f = atof (str);

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




static void
enter_float_parameter  (float *par, 
                        int fail_state,
                        const struct rx_options *rx_opts,
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


static enum ColPartner
conv_name_to_int (const char *str)
{
  enum ColPartner res = He;
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

static bool
isAllowedCollisionPartners (const char *str, struct col_partner cps[7], 
                                  size_t *s, const struct rx_options *rx_opts)
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
              cps[numof_cps].name = conv_name_to_int (t);
          else if (numof_t == 1)
            {
              float d = atof (t);
              cps[numof_cps].dens = atof (t);
              if (rx_opts->dens_log_scale)
                {
                  if (d >= 0 && d <= 14)
                    cps[numof_cps].dens = powf (10, d);
                  else
                    res = false;
                }
              else
                {
                  if (d >= 1 && d <= 1e14)
                    cps[numof_cps].dens = d;
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

  *s = numof_cps;
  free (storage);
  return res;
}


static void
enter_collision_partners (struct col_partner *cps[7], 
                          size_t *s,
                          int fail_state,
                          const struct rx_options *rx_opts,
                          const char *history_file, 
                          linenoiseCompletionCallback *completion,
                          linenoiseHintsCallback *hints,
                          allowed_value_cps isAllowed)
{
  size_t size = 30;
  char *line;
  line = (char *) malloc (size);
  linenoiseHistoryLoad (history_file);
  linenoiseSetCompletionCallback (completion);
  linenoiseSetHintsCallback (hints);

  while ((line = linenoise ("  >> ")) != NULL)
    {
      if (isAllowed (line, cps, s, rx_opts))
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

void
start_dialogue (float *sfreq, float *efreq, struct MC_parameters *mc_pars, 
                                               size_t *s,
                                               const struct rx_options *opts) 
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
 
  enter_molecule (mc_pars);
      
  printf ("\x1B[0;37;40m\x1B[1;35;40m  ## \x1B[0;37;40m");
  printf ("Enter starting & ending frequencies ");
  printf ("\x1B[1;37;40m[GHz]\x1B[0;37;40m\x1B[3;37;40m\n");

  enter_frequencies (sfreq, efreq);

  printf ("\x1B[0;37;40m\x1B[1;35;40m  ## \x1B[0;37;40m");
  printf ("Enter kinetic temperature ");
  printf ("\x1B[1;37;40m[K]\x1B[0;37;40m\x1B[3;37;40m\n");

  enter_float_parameter (&mc_pars->Tkin, TKIN_FAIL, opts, ".Tkin", NULL, 
                                                  hints_temp, isAllowedTemp);

  printf ("\x1B[0;37;40m\x1B[1;35;40m  ## \x1B[0;37;40m");
  printf ("Enter background temperature ");
  printf ("\x1B[1;37;40m[K]\x1B[0;37;40m\x1B[3;37;40m\n");

  enter_float_parameter (&mc_pars->Tbg, TBG_FAIL, opts, ".Tbg", NULL, 
                                                  hints_temp, isAllowedTemp);

  printf ("\x1B[0;37;40m\x1B[1;35;40m  ## \x1B[0;37;40m");
  printf ("Enter column density for the molecule ");
  printf ("\x1B[1;37;40m[cm-2]\x1B[0;37;40m\n");

  enter_float_parameter (&mc_pars->coldens, COLDENS_FAIL, opts, ".coldens",
                                    NULL, NULL, isAllowedColdens);

  printf ("\x1B[0;37;40m\x1B[1;35;40m  ## \x1B[0;37;40m");
  printf ("Enter FWHM lines width ");
  if (opts->hz_width) 
    printf ("\x1B[1;37;40m[Hz]\x1B[0;37;40m\n");
  else 
    printf ("\x1B[1;37;40m[km s-1]\x1B[0;37;40m\n");

  enter_float_parameter (&mc_pars->line_width, FWHM_FAIL, opts, ".fwhm",
                                                NULL, NULL, isAllowedWidth);

  printf ("\x1B[0;37;40m\x1B[1;35;40m  ## \x1B[0;37;40m");
  printf ("Enter collision partners and their densities ");
  printf ("\x1B[1;37;40m[cm-3]\x1B[0;37;40m\n");

  enter_collision_partners (&mc_pars->cps, s, COLLISION_FAIL, opts,
      ".collisions", NULL, NULL, isAllowedCollisionPartners); 
}
