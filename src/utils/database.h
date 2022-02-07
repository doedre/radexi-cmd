/**
 * @file utils/database.h
 * @brief Manipulating local database files.
 */

#ifndef RXI_DATABASE_H
#define RXI_DATABASE_H

#include <dirent.h>

#include "rxi_common.h"

/// @brief Add LAMDA's molecular database file to the local database.
///
/// Function for `--add-molecule` option. Write .csv files to
/// `/home/$(USER)/.local/share/radexi/`.
/// @param *name -- molecule name defined by the user;
/// @param *path -- path to the LAMDA's database file.
/// @return `RXI_OK` on success, `RXI_ERR_ALLLOC` on allocation error,
/// `RXI_WARN_LAMDA` on possible errors in LAMDA database file.
RXI_STAT rxi_add_molecule (const char *name, const char *path);

/// @brief Remove molecule from the local database.
///
/// Function for `--delete-molecule` option. Recursively delete all files
/// in specified molecule's folder in local database
/// (`/home/$(USER)/.local/share/radexi/`).
/// @param *name -- molecule name for deletion.
/// @return `RXI_OK` on success, `RXI_ERR_ALLOC` on allocation error,
/// `RXI_WARN_NOFILE` if no such molecule exist.
RXI_STAT rxi_delete_molecule (const char *name);

/// @brief Prints molecule names from the local database.
///
/// TODO
RXI_STAT rxi_list_molecules ();

/// @brief Converts collisional partner number to string.
///
/// Allocates memory for the returned string, so it should be freed after
/// usage. Convert `COLL_PART` to predefined string names.
/// @param cp -- collision partner of `enum COLL_PART` type.
/// @return Allocated string with collision partner's name.
char *numtoname (COLL_PART cp);

/// @brief Converts string to collisional partner number.
///
/// Converts collision partner's name to `enum COLL_PART`.
/// @param *name -- string with name.
/// @return One of the collisional partners or `NO_MOLECULE` on error.
COLL_PART nametonum (const char *name);

/// @brief Get number for specified collision partner from database info file.
///
/// Searches database (through `struct rxi_db_molecule_info`) for specified
/// collision partner number in database (not `enum COLL_PART`).
/// @param *mol_info -- structure with information about the molecule from
/// local database;
/// @param cp -- collision partner which number needs to be defined.
/// @return Count number for specified collision partner.
int8_t cptonum (const struct rxi_db_molecule_info *mol_info, COLL_PART cp);

/// @brief Iterate through local database molecule names.
///
/// Given database path iterates through it and fills second parameter with
/// new name. Return value depends on the type of file with gotten name.
/// @param *dir -- directory stream opened before usage;
/// @param *name -- string with allocated memory to write filenames into.
/// @return `1` on folder, `-1` on hidden files, `-2` on regular files,
/// `0` on end of file stream.
int rxi_db_molecule_iter (DIR *dir, char *name);

/// @brief Reads molecular info file.
///
/// Searches through local database to fill `struct rxi_db_molecule_info`.
/// @param *name -- molecule name;
/// @param *mol_info -- allocated structure to hold information.
/// @return `RXI_OK` on success, `RXI_ERR_ALLOC` on memory allocation error,
/// `RXI_ERR_FILE` on file errors.
RXI_STAT rxi_db_read_molecule_info (const char *name,
                                    struct rxi_db_molecule_info *mol_info);

/// @brief Reads energy level file.
///
/// Searches through local database to fill `struct rxi_db_molecule_enlev`.
/// @param *name -- molecule name;
/// @param *mol_enlev -- allocated structure to hold information.
/// @return `RXI_OK` on success, `RXI_ERR_ALLOC` on memory allocation error,
/// `RXI_ERR_FILE` on file errors.
RXI_STAT rxi_db_read_molecule_enlev (const char *name,
                                     struct rxi_db_molecule_enlev *mol_enlev);

/// @brief Reads radiation transfer file.
///
/// Searches through local database to fill `struct rxi_db_molecule_radtr`.
/// @param *name -- molecule name;
/// @param *mol_enlev -- allocated structure to hold information.
/// @return `RXI_OK` on success, `RXI_ERR_ALLOC` on memory allocation error,
/// `RXI_ERR_FILE` on file errors.
RXI_STAT rxi_db_read_molecule_radtr (const char *name,
                                     struct rxi_db_molecule_radtr *mol_radtr);

/// @brief Reads collision partner file.
///
/// Searches through local database to fill `struct rxi_db_molecule_coll_part`.
/// @param *mol_name -- molecule name;
/// @param cp -- collision partner's name;
/// @param n_temps -- number of collisional temperatures for specified
/// collisional partner.
/// @param *mol_cp -- allocated structure to hold information.
/// @return `RXI_OK` on success, `RXI_ERR_ALLOC` on memory allocation error,
/// `RXI_ERR_FILE` on file errors.
RXI_STAT rxi_db_read_molecule_coll_part (const char *mol_name,
      const COLL_PART cp, const size_t n_temps,
      struct rxi_db_molecule_coll_part *mol_cp);

#endif  // RXI_DATABASE_H
