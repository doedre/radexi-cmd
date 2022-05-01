/**
 * @file core/dialogue.h
 * @brief Common dialogue entries.
 */

#include "rxi_common.h"

/// @brief Typical dialog with the user to get input information.
///
/// Initiates dialogue with the user to collect starting parameters. Mostly
/// use tools from `utils/cli_tools.h`.
/// @param *inp_data -- structure to fill;
/// @param *opts -- options from command line.
/// @return `RXI_OK` on success.
RXI_STAT rxi_dialog_input (struct rxi_input_data *inp_data);

/// @brief Dialog to get information about the spectral lines.
///
/// Used in finding the best fit using escape probability method.
RXI_STAT
rxi_dialog_best_fit (struct rxi_input_data *inp_data,
                     struct rxi_db_molecule_info **minfo,
                     struct rxi_db_molecule_enlev **enlev,
                     struct rxi_db_molecule_radtr **radtr);
