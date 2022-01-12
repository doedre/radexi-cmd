/**
 * @file utils/database.h
 * @brief Manipulating local database files.
 */

#include <dirent.h>

#include "rxi_common.h"

/// @brief Add LAMDA's molecular database file to the local database.
RXI_STAT rxi_add_molecule (const char *name, const char *path);

/// @brief Remove molecule from the local database.
RXI_STAT rxi_delete_molecule (const char *name);

/// @brief Prints molecule names from the local database.
RXI_STAT rxi_list_molecules ();

/// @brief Converts collisional partner number to string.
char *numtoname (COLL_PART cp);

/// @brief Converts string to collisional partner number.
COLL_PART nametonum (const char *name);

/// @brief Helps to iterate through local database molecule names.
int rxi_db_molecule_iter (DIR *dir, char *name);

/// @brief Reads molecular info file.
RXI_STAT rxi_db_read_molecule_info (const char *name,
                                    struct rxi_db_molecule_info *mol_info);

/// @brief Reads energy levels file.
RXI_STAT rxi_db_read_molecule_enlev (const char *name,
                                     struct rxi_db_molecule_enlev *mol_enlev);
