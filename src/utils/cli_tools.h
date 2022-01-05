/**
 * @file utils/cli_tools.h
 * @brief Controls colored output and readline functionality.
 */

#include <stdbool.h>

#include "rxi_common.h"

/// @brief Wrap for `linenoise()`.
char *rxi_readline (const char *prompt);

/// @brief Accept for `[Y/n]` case.
bool rxi_readline_accept ();

/// @brief Save current history to specified file.
RXI_STAT rxi_history_save (const char *line, const char *filename);

/// @brief Reallocate buffer from linenoise to load only specified history.
RXI_STAT rxi_history_load (const char *filename);
