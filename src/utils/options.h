/**
 * @file options.h
 * @brief Wrap for `getopt_long()` (defined in `getopt.h`) to extract options. 
 */

#ifndef RXI_OPTIONS_H
#define RXI_OPTIONS_H

#include <getopt.h>
#include <stddef.h>
#include <limits.h>

#include "rxi_common.h"

/// @brief Defines short-like option for long options w/o short equivalent.
/// 
/// Starts with `CHAR_MAX` (defined in `limits.h`) to assert that no short
/// option will be used instead. Look for @ref `long_options[]` array for
/// additional information (or for `getopt.h` `man` pages).
enum 
{
  ADD_MOLECULE_OPTION = CHAR_MAX + 1,
  LIST_MOLECULES_OPTION,
  DELETE_MOLECULE_OPTION,
  VERSION_OPTION
};

/// @brief Sets all options to their default values.
///
/// For default values see `struct rxi_options`.
/// @param *opts -- pointer to used `struct rxi_options`
void rxi_set_default_options (struct rxi_options *opts);

/// @brief Sets options according to command line arguments.
///
/// Using `getopt_long()` to parse command line arguments.
/// @param *opts -- pointer to used `struct rxi_options`;
/// @param argc -- number of command line options;
/// @param **argv -- command line arguments.
/// @return `optind` argument from command line, which defines free argument.
int rxi_set_options (struct rxi_options *opts, int argc, char **argv);

#endif  // RXI_OPTIONS_H
