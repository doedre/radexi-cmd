/*	options.c
 *
 * Uses 'getopt' to set options for running this program. Long and short 
 * options are produced.
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

#include <getopt.h>   // getopt_long() defined here
#include <limits.h>   // CHAR_MAX defined here
#include <stddef.h>   // NULL defined here

#include "radexi.h"

/* Defines short-like option for long options w/o short equivalent. Starts w/
 * CHAR_MAX to assert that no short option will be used instead.  */
enum 
{
  ADD_MOLECULE_OPTION = CHAR_MAX + 1,
  LIST_MOLECULES_OPTION,
  DELETE_MOLECULE_OPTION,
  VERSION_OPTION
};

/* Defines long options for usage by getopt_long(). */
static struct option const long_options[] = {
//  {"add-molecule", required_argument, NULL, ADD_MOLECULE_OPTION},
//  {"list-molecules", no_argument, NULL, LIST_MOLECULES_OPTION},
//  {"delete-molecule", required_argument, NULL, DELETE_MOLECULE_OPTION},
  {"result", required_argument, NULL, 'r'},
  {"help", no_argument, NULL, 'h'},
  {"version", no_argument, NULL, VERSION_OPTION}
};

int
set_rx_options (struct rx_options *opts, int argc, char ** argv)
{
  set_default_rx_options (opts);
  int option_index = 0;
  int opt;
  while ((opt = getopt_long (argc, argv, 
                              "hor:", 
                              long_options, &option_index)) != -1)
    {
      switch (opt)
        {
        case 'h':
          cmd_usage (EXIT_SUCCESS);
          break;

        case 'o':
          break;

        case 'r':
          opts->user_defined_out_file_path = true;
          break;

        case ADD_MOLECULE_OPTION:
          break;

        case LIST_MOLECULES_OPTION:
          break;

        case DELETE_MOLECULE_OPTION:
          break;

        case VERSION_OPTION:
          break;

        default:
          cmd_usage (RADEXI_FAILURE);

        }
    }

  return optind;
}

void 
set_default_rx_options (struct rx_options *opts)
{
}

void 
cmd_usage (int status)
{
  if (status != EXIT_SUCCESS) 
    {
      fprintf (stderr, ("Try '%s --help' for more information.\n"), 
                PROGRAM_NAME);
    }
  else
    {
      printf ("--help option is in a maintanence...\n");
    }

  exit (status);
}
