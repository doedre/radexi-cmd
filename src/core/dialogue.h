/**
 * @file core/dialogue.h
 * @brief Common dialogue entries.
 */

#include "rxi_common.h"

/// @brief Typical dialog with the user to get input information.
RXI_STAT rxi_dialog_input (struct rxi_input_data *inp_data,
                           const struct rxi_options *opts);
