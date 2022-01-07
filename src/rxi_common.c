/**
 * @file rxi_common.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "rxi_common.h"

#include <utils/debug.h>

const char*
rxi_database_path ()
{
  DEBUG ("Form database path");

  char *db_path = malloc (RXI_PATH_MAX * sizeof (*db_path));
  CHECK (db_path);
  if (!db_path)
    return NULL;

  const char *home_path = getenv ("HOME");
  CHECK (home_path);
  if (!home_path)
    return NULL;

  strcpy (db_path, home_path);
  strcat (db_path, "/.local/share/radexi/");

  ASSERT (db_path);
  return db_path;
}

const char*
rxi_config_path ()
{
  DEBUG ("Form config path");

  char *config_path = malloc (RXI_PATH_MAX * sizeof (*config_path));
  CHECK (config_path);
  if (!config_path)
    return NULL;

  const char *home_path = getenv ("HOME");
  CHECK (home_path);
  if (!home_path)
    return NULL;

  strcpy (config_path, home_path);
  strcat (config_path, "/.config/radexi/");

  ASSERT (config_path);
  return config_path;
}

RXI_STAT
rxi_db_molecule_info_malloc (struct rxi_db_molecule_info **mol_info)
{
  struct rxi_db_molecule_info *mi = malloc (sizeof (*mi));
  CHECK (mi);
  if (!mi)
    goto malloc_error;

  char *name = malloc (RXI_STRING_MAX * sizeof (*name));
  CHECK (name);
  if (!name)
    {
      free (mi);
      goto malloc_error;
    }

  mi->name = name;
  *mol_info = mi;

  return RXI_OK;

malloc_error:
  *mol_info = NULL;
  return RXI_ERR_ALLOC;
}

void
rxi_db_molecule_info_free (struct rxi_db_molecule_info *mol_info)
{
  free (mol_info->name);
  free (mol_info);
}
