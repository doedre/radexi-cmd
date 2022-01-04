/**
 * @file options.c
 */

#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <getopt.h>

#include "options.h"

#include "rxi_common.h"
#include "utils/debug.h"

/// @brief Defines all possible command line options.
///
/// Look for `getopt.h` `man` pages for additional information.
struct option const long_options[] = {
  {"add-molecule",  required_argument,  NULL, ADD_MOLECULE_OPTION},
//  {"list-molecules", no_argument, NULL, LIST_MOLECULES_OPTION},
//  {"delete-molecule", required_argument, NULL, DELETE_MOLECULE_OPTION},
  {"log-density",   no_argument,        NULL, 'L'},
  {"hz-width",      no_argument,        NULL, 'H'},
  {"help",          no_argument,        NULL, 'h'},
  {"result",        required_argument,  NULL, 'r'},
  {"version",       no_argument,        NULL, VERSION_OPTION}
};


void 
rxi_set_default_options (struct rxi_options *opts)
{
  DEBUG ("Set options to default values");
  opts->usage_mode = UM_NONE;
  opts->status = RXI_OK;
  opts->force_fs = false;
  opts->quite_start = false;
  opts->cmd_output = false;
  opts->dens_log_scale = false;
  opts->hz_width = false;
  opts->user_defined_out_file_path = false;
  strcpy (opts->result_path, ".");
}

int
rxi_set_options (struct rxi_options *opts, int argc, char **argv)
{
  DEBUG ("Set options from command line");
  rxi_set_default_options (opts);
  int option_index = 0;
  int opt;
  while ((opt = getopt_long (argc, argv, 
                              ":fLqHhor:", 
                              long_options, &option_index)) != -1)
    {
      switch (opt)
        {
        case 'L':
          DEBUG ("Set -L (--log-density) option");
          opts->dens_log_scale = true;
          break;

        case 'H':
          DEBUG ("Set -H (--hz-width) option");
          opts->hz_width = true;
          break;

        case 'h':
          DEBUG ("Set -h (--help) option");
          opts->usage_mode = UM_HELP;
          break;

        case 'o':
          DEBUG ("Set -o option");
          break;

        case 'q':
          DEBUG ("Set -q option");
          opts->quite_start = true;
          break;

        case 'r':
          DEBUG ("Set -r (--result) option");
          opts->user_defined_out_file_path = true;
          strcpy (opts->result_path, optarg);
          break;

        case 'f':
          DEBUG ("Set -f option");
          opts->force_fs = true;
          break;

        case ADD_MOLECULE_OPTION:
          DEBUG ("Set --add-molecule option");
          if (opts->usage_mode != UM_NONE)
            break;

          opts->usage_mode = UM_MOLECULAR_FILE_ADD;
          strcpy (opts->molecule_name, optarg);
          break;

        case LIST_MOLECULES_OPTION:
          DEBUG ("Set --list-molecule option");
          if (opts->usage_mode != UM_NONE)
            break;

          opts->usage_mode = UM_MOLECULAR_FILE_LIST;
          break;

        case DELETE_MOLECULE_OPTION:
          DEBUG ("Set --delete-molecule option");
          if (opts->usage_mode != UM_NONE)
            break;

          opts->usage_mode = UM_MOLECULAR_FILE_DELETE;
          break;

        case VERSION_OPTION:
          DEBUG ("Set --version option");
          opts->usage_mode = UM_VERSION;
          break;

        case '?':
          DEBUG ("Unknown option");
          fprintf (stderr, "Unknown option was used\n");
          opts->usage_mode = UM_HELP;
          opts->status = RXI_ERR_OPTS;
          break;

        case ':':
          DEBUG ("No argument for option");
          fprintf (stderr, "No argument for `%s' option was specified\n",
                   long_options[option_index].name);
          opts->usage_mode = UM_HELP;
          opts->status = RXI_ERR_OPTS;
          break;

        default:
          DEBUG ("Unspecified error in option handling");
          fprintf (stderr, "Unknown option was used\n");
          opts->usage_mode = UM_HELP;
          opts->status = RXI_ERR_OPTS;
          break;
        }
    }
  
  if ((optind == -1) && (opts->usage_mode == UM_NONE))
    opts->usage_mode = UM_FILE;
  else if (opts->usage_mode == UM_NONE)
    opts->usage_mode = UM_DIALOGUE;

  DEBUG ("Usage mode: %u", opts->usage_mode);
  ASSERT ((opts->usage_mode != UM_NONE) && "Usage mode undefined");

  return optind;
}
