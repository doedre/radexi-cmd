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
#include <stddef.h>   // NULL defined here
#include <string.h>
#include <unistd.h>

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
  {"add-molecule",  required_argument,  NULL, ADD_MOLECULE_OPTION},
//  {"list-molecules", no_argument, NULL, LIST_MOLECULES_OPTION},
//  {"delete-molecule", required_argument, NULL, DELETE_MOLECULE_OPTION},
  {"log-density",   no_argument,        NULL, 'L'},
  {"hz-width",      no_argument,        NULL, 'H'},
  {"help",          no_argument,        NULL, 'h'},
  {"result",        required_argument,  NULL, 'r'},
  {"version",       no_argument,        NULL, VERSION_OPTION}
};

int
set_rxi_options (struct rxi_options *opts, int argc, char ** argv)
{
  // Making default path to the radexi folder
  strcat (radexi_path, "/home/");
  char *username = malloc (32 * sizeof (char));
  if (!getlogin_r (username, 32))
    strcat (radexi_path, username);
  else 
    return -2;
  strcat (radexi_path, "/radexi");

  set_default_rxi_options (opts);
  int option_index = 0;
  int opt;
  while ((opt = getopt_long (argc, argv, 
                              "flqHhor:", 
                              long_options, &option_index)) != -1)
    {
      switch (opt)
        {
        case 'L':
          opts->dens_log_scale = true;
          break;

        case 'H':
          opts->hz_width = true;
          break;

        case 'h':
          cmd_usage (EXIT_SUCCESS);
          break;

        case 'o':
          break;

        case 'q':
          opts->quite_start = true;
          break;

        case 'r':
          opts->user_defined_out_file_path = true;
          strcpy (opts->result_path, optarg);
          break;

        case 'f':
          opts->force_fs = true;
          break;

        case ADD_MOLECULE_OPTION:
          opts->usage_mode = UM_MOLECULAR_FILE_ADD;
          strcpy (opts->molecule_name, optarg);
          break;

        case LIST_MOLECULES_OPTION:
          opts->usage_mode = UM_MOLECULAR_FILE_LIST;
          break;

        case DELETE_MOLECULE_OPTION:
          opts->usage_mode = UM_MOLECULAR_FILE_DELETE;
          break;

        case VERSION_OPTION:
          break;

        default:
          cmd_usage (RADEXI_FAILURE);
        }
    }
  
  if (optind == -1)
    opts->usage_mode = UM_FILE;

  return optind;
}

void 
set_default_rxi_options (struct rxi_options *opts)
{
  opts->usage_mode = UM_DIALOGUE;
  opts->force_fs = false;
  opts->quite_start = false;
  opts->cmd_output = false;
  opts->dens_log_scale = false;
  opts->hz_width = false;
  opts->user_defined_out_file_path = false;
  strcat (opts->result_path, ".");
}

void 
cmd_usage (int status)
{
  if (status != EXIT_SUCCESS) 
    {
      fprintf (stderr, ("Try '%s --help' for more information.\n"), PROGRAM_NAME);
    }
  else
    {
      printf ("--help option is in a maintanence...\n");
    }

  exit (status);
}
