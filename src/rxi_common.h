/**
 * @file rxi_common.h
 * @brief Provides common data structures for calculations.
 */

#ifndef RXI_COMMON_H
#define RXI_COMMON_H

#include <stdint.h>
#include <stdbool.h>

#define RXI_PATH_MAX 1024
#define RXI_MOLECULE_MAX 15

/// @brief Status codes for functions that may fail.
///
/// Status codes should be returned by all of the functions, that may fail
/// during memory allocation, filesystem errors or just local mistakes in
/// molecular databases. They are divided by two types: warnings and errors.
/// Errors may lead to unexpected behavior. Warnings are expected not to
/// terminate the program with segmentation faults.
typedef enum RXI_STAT
{
  RXI_OK = 0,             //!< Everything is ok.
  RXI_ERR_ALLOC,          //!< Error on memory allocation.
  RXI_ERR_OPTS,
  RXI_WARN_LIMITS = 10    //!< 
}
RXI_STAT;

/// @brief Used to specify how the program will be used.
///
/// Only used in @ref `struct rxi_options` to set the program in specified
/// state.
enum USAGE_MODE
{
  UM_NONE = 0,                //!< Undefined usage mode.
  UM_DIALOGUE,                //!< Standart dialog with the user.
  UM_FILE,                    //!< Read input parameters from the file.
  UM_MOLECULAR_FILE_ADD,      //!< Add molecular data file from LAMDA.
  UM_MOLECULAR_FILE_DELETE,   //!< Delete local molecular data file.
  UM_MOLECULAR_FILE_LIST,     //!< List local molecular data files.
  UM_HELP,                    //!< Print help information.
  UM_VERSION                  //!< Print version information.
};


/* This structure defines options for all of the possible processes in the 
 * program. It is being filled after reading options by 'getopt_long()' in 
 * 'set_rx_options' functions. */
struct rxi_options
{
  /* How molecular cloud's parameters will be collected: via dialogue or file.
   */
  enum USAGE_MODE usage_mode;

  enum RXI_STAT status;

  /* Molecule name defined by the user when adding one in case of 
   * UM_MOLECULAR_FILE usage_mode.  */
  char molecule_name[RXI_MOLECULE_MAX];

  /* Forces file rewrite if specified files for results of calculation 
   * already exist.
   * -f flag: disabled by default */
  bool force_fs;

  /* Whether resulting table should be printed in the terminal or not.
   * -o flag: disabled by default */
  bool cmd_output;

  /* Print message w/ information on the start or not.
   * -q flag: disabled by default */
  bool quite_start;

  /* If 'true', then instead of writing densities in scientific format user 
   * should only write power of 10.
   * -L of --log-density flag: disabled by default */
  bool dens_log_scale;

  /* if 'true', then use [GHz] line width instead of [km -s].
   * -H og --hz-width flag: disabled by default  */
  bool hz_width;

  /* If user defined a path for the file w/ resilts. 
   * -r or --result flag */
  bool user_defined_out_file_path;

  /* Specified path for user's output file  */
  char result_path[RXI_PATH_MAX];
};

/// @brief Get `$(HOME)/.local/share/radexi/` path.
///
/// @return On success returns `$(HOME)/.local/share/radexi/` path string. On
/// error returns `NULL`.
const char *rxi_database_path ();

/// @brief Get `$(HOME)/.config/radexi/` path.
///
/// @return On success returns `$(HOME)/.config/radexi/` path string. On error
/// returns `NULL`.
const char *rxi_config_path ();

#endif  // RXI_COMMON_H
