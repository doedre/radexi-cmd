/**
 * @file output.h
 * @brief Create convinient output.
 */

#include "rxi_common.h"

/// @brief Writes output depending on the user options.
RXI_STAT rxi_out_result (struct rxi_calc_data *data[RXI_MOLECULE_MAX],
                         const struct rxi_options *opts);

/// @brief Checks if specified file can be written.
///
/// @param *path -- specified path.
/// @return `0` if specified file can be written; `1` if it contains a regular
/// file; `2` if it is a directory; `-1` if file already exists.
int rxi_check_path (const char *path);
