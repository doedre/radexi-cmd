/**
 * @file utils/database.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

#include "utils/database.h"

#include "rxi_common.h"
#include "utils/cli_tools.h"
#include "utils/debug.h"

#include "minIni/minIni.h"

static RXI_STAT
rxi_add_molecule_info (const char *name, const char *db_path, FILE *molfile)
{
  char *info_filename = malloc (RXI_PATH_MAX * sizeof (*info_filename));
  CHECK (info_filename && "Allocation error");
  if (!info_filename)
    return RXI_ERR_ALLOC;

  char *line = malloc (RXI_STRING_MAX * sizeof (*line));
  CHECK (line && "Allocation error");
  if (!line)
    {
      free (info_filename);
      return RXI_ERR_ALLOC;
    }

  strcpy (info_filename, db_path);
  strcat (info_filename, name);
  strcat (info_filename, "/");
  strcat (info_filename, name);
  strcat (info_filename, ".info");

  DEBUG ("Write information to `%s' file", info_filename);

  int8_t ini_stat = 1;
  for (int8_t i = 0; i < 7; i++)
    {
      if (!fgets (line, RXI_STRING_MAX, molfile))
        goto file_error;

      if (i == 1)
        ini_stat = ini_puts ("Information", "name", line, info_filename);
      else if (i == 3)
        ini_stat = ini_puts ("Information", "weight", line, info_filename);
      else if (i == 5)
        ini_stat = ini_puts ("Information", "energy_levels",
                             line, info_filename);

      if (ini_stat == 0)
        goto file_error;
    }
  free (line);
  free (info_filename);
  return RXI_OK;

file_error:
  free (line);
  free (info_filename);
  return RXI_ERR_FILE;
}

RXI_STAT rxi_add_molecule (const char *name, const char *path)
{
  DEBUG ("Start add molecule '%s' from `%s' to local database", name, path);

  RXI_STAT status = RXI_OK;

  FILE *molfile = fopen (path, "r");
  CHECK (molfile && "Error opening LAMDA database file");
  if (!molfile)
    return RXI_ERR_FILE;

  const char *db_path = rxi_database_path ();
  CHECK (db_path && "Allocation error");
  if (!db_path)
    {
      fclose (molfile);
      return RXI_ERR_ALLOC;
    }

  char *db_folder = malloc (RXI_PATH_MAX * sizeof (*db_folder));
  CHECK (db_folder && "Allocation error");
  if (!db_folder)
    {
      fclose (molfile);
      free ((void*)db_path);
      return RXI_ERR_ALLOC;
    }
  strcpy (db_folder, db_path);
  strcat (db_folder, name);

  DEBUG ("Local database folder `%s'", db_folder);

  // Check if local database folder with the same name exist
  struct stat sb;
  if (stat (db_folder, &sb) == -1)
    {
      DEBUG ("Creating new folder");
      if (mkdir (db_folder, 0700) != 0)
        goto file_error;
    }
  else 
    {
      DEBUG ("Ask to rewrite this folder");
      printf ("  ## "
              "Specified molecule name already exists, rewrite it?\n");
      if (!rxi_readline_accept ())
        goto file_error;
    }

  status = rxi_add_molecule_info (name, db_path, molfile);
  CHECK ((status == RXI_OK));

  free (db_folder);
  free ((void*)db_path);
  fclose (molfile);

  return status;

file_error:
  free (db_folder);
  free ((void*)db_path);
  fclose (molfile);
  return RXI_ERR_FILE;
}
