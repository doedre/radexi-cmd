/**
 * @file utils/database.h
 * @brief Manipulating local database files.
 */

#include "rxi_common.h"

/// @brief Add LAMDA's molecular database file to the local database.
RXI_STAT rxi_add_molecule (const char *name, const char *path);

/// @brief Remove molecule from the local database.
RXI_STAT rxi_delete_molecule (char *name);
