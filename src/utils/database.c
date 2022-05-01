/**
 * @file utils/database.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <gsl/gsl_matrix.h>

#include "utils/database.h"

#include "rxi_common.h"
#include "utils/cli_tools.h"
#include "utils/csv.h"
#include "utils/debug.h"

#include "minIni/minIni.h"

static void
normalise_lamda_comment (char *comment)
{
  CHECK (comment && "NULL comment passed");
  if (!comment)
    return;

  if (comment[0] != '!')
    return;

  DEBUG ("Comment before normalising: %s", comment);

  for (int i = 0; comment[i] != '\0'; ++i)
    {
      if (comment[i] == '\n')
        comment[i] = ' ';

      if (isupper (comment[i]))
        comment[i] = tolower (comment[i]);

      if (comment[i] == ' ')
        {
          memmove (comment + i, comment + i + 1, RXI_STRING_MAX - i - 1);
          --i;
        }
    }

  DEBUG ("Comment after normalising: %s", comment);
}

static RXI_STAT
rxi_save_molecule_info (FILE *molfile,
                        struct rxi_db_molecule_info *molecule_info,
                        int8_t n_partner)
{
  char *line = malloc (RXI_STRING_MAX * sizeof (*line));
  CHECK (line && "Allocation error");
  if (!line)
    return RXI_ERR_ALLOC;

  RXI_STAT status = RXI_OK;
  while (fgets (line, RXI_STRING_MAX, molfile))
    {
      normalise_lamda_comment (line);
      if (!strncmp (line, "!molecule", 9))
        {
          if (!fgets (line, RXI_STRING_MAX, molfile))
            goto file_error;

          memcpy (molecule_info->name, line, RXI_STRING_MAX);
        }
      else if (!strncmp (line, "!molecularweight", 16)
               || !strncmp (line, "!mass", 5))
        {
          if (!fgets (line, RXI_STRING_MAX, molfile))
            goto file_error;

          molecule_info->weight = strtof (line, NULL);

          if (molecule_info->weight == 0)
            goto conversion_error;
        }
      else if (!strncmp (line, "!numberofenergylevels", 21))
        {
          if (!fgets (line, RXI_STRING_MAX, molfile))
            goto file_error;

          molecule_info->numof_enlev = strtol (line, NULL, 10);

          if (molecule_info->numof_enlev == 0)
            goto conversion_error;
        }
      else if (!strncmp (line, "!numberofradiativetransitions", 29))
        {
          if (!fgets (line, RXI_STRING_MAX, molfile))
            goto file_error;

          molecule_info->numof_radtr = strtol (line, NULL, 10);

          if (molecule_info->numof_radtr == 0)
            goto conversion_error;
        }
      else if (!strncmp (line, "!numberofcollpartners", 21)
               || !strncmp (line, "!numberofcollisionpartners", 26))
        {
          if (!fgets (line, RXI_STRING_MAX, molfile))
            goto file_error;

          molecule_info->numof_coll_part = strtol (line, NULL, 10);

          if (molecule_info->numof_coll_part == 0)
            goto conversion_error;
        }
      else if (!strncmp (line, "!partner", 7)
               || !strncmp (line, "!collisionsbetween", 18)
               || !strncmp (line, "!collisionpartner", 17))
        {
          if (!fgets (line, RXI_STRING_MAX, molfile))
            goto file_error;

          molecule_info->coll_part[n_partner] = strtol (&line[0], NULL, 10);

          if (molecule_info->coll_part[n_partner] == 0)
            goto conversion_error;
        }
      else if (!strncmp (line, "!numberofcolltrans", 18)
               || !strncmp (line, "!numberofcollisionaltransitions", 31))
        {
          if (!fgets (line, RXI_STRING_MAX, molfile))
            goto file_error;

          molecule_info->numof_coll_trans[n_partner] = strtol (line, NULL, 10);

          if (molecule_info->numof_coll_trans[n_partner] == 0)
            goto conversion_error;
        }
      else if (!strncmp (line, "!numberofcolltemps", 18)
               || !strncmp (line, "!numberofcollisionaltemperatures", 32)
               || !strncmp (line, "!numberofcollisiontemperatures", 30))
        {
          if (!fgets (line, RXI_STRING_MAX, molfile))
            goto file_error;

          molecule_info->numof_coll_temps[n_partner] = strtol (line, NULL, 10);

          if (molecule_info->numof_coll_temps[n_partner] == 0)
            goto conversion_error;
        }
      else if (!strncmp (line, "!colltemps", 10)
               || !strncmp (line, "!collisionaltemperatures", 24)
               || !strncmp (line, "!collisiontemperatures", 22))
        {
          if (!fgets (line, RXI_STRING_MAX, molfile))
            goto file_error;

          char *start = malloc (RXI_STRING_MAX * sizeof (*start));
          memcpy (start, line, RXI_STRING_MAX);
          char *end;
          int8_t i = 0;
          for (float f = strtof (start, &end);
               start != end;
               f = strtof (start, &end))
            {
              if (f == 0)
                goto conversion_error;

              start = end;
              DEBUG ("Number for writing: %f", f);
              gsl_matrix_set (molecule_info->coll_temps, n_partner, i, f);
              ++i;
            }

          if (molecule_info->numof_coll_temps[n_partner] != i)
            status = RXI_WARN_LAMDA;
        }
      else
        {
          break;
        }
    }

  free (line);
  return status;

file_error:
  free (line);
  return RXI_ERR_FILE;

conversion_error:
  free (line);
  return RXI_ERR_CONV;
}

static RXI_STAT
check_db_molecule_info (const struct rxi_db_molecule_info *mol_info)
{
  if (!mol_info)
    return RXI_ERR_ALLOC;

  if (mol_info->weight < 1 || mol_info->numof_enlev == 0
      || mol_info->numof_radtr == 0 || mol_info->numof_coll_part == 0)
    return RXI_WARN_LAMDA;

  if (mol_info->numof_coll_temps[0] == 0
      || mol_info->numof_coll_trans[0] == 0)
    return RXI_WARN_LAMDA;

  return RXI_OK;
}

static RXI_STAT
rxi_add_molecule_info (const char *db_folder, const char *name,
                       const struct rxi_db_molecule_info *mol_info)
{
  if (check_db_molecule_info (mol_info) != RXI_OK)
    return RXI_WARN_LAMDA;

  char *info_filename = malloc (RXI_PATH_MAX * sizeof (*info_filename));
  CHECK (info_filename && "Allocation error");
  if (!info_filename)
    return RXI_ERR_ALLOC;

  strcpy (info_filename, db_folder);
  strcat (info_filename, "/");
  strcat (info_filename, name);
  strcat (info_filename, ".info");

  DEBUG ("Write information to `%s' file", info_filename);

  int8_t ini_stat = 0;
  ini_stat += !ini_puts ("Information", "name", mol_info->name, info_filename);
  ini_stat += !ini_putf ("Information", "weight", mol_info->weight,
                        info_filename);
  ini_stat += !ini_putl ("Information", "energy_levels", mol_info->numof_enlev,
                         info_filename);
  ini_stat += !ini_putl ("Information", "radiative_transitions",
                         mol_info->numof_radtr, info_filename);
  ini_stat += !ini_putl ("Information", "collision_partners",
                         mol_info->numof_coll_part, info_filename);

  for (int8_t i = 0; i < mol_info->numof_coll_part; ++i)
    {
      char *cp_name = numtoname (mol_info->coll_part[i]);
      if (!cp_name)
        {
          free (info_filename);
          return RXI_ERR_CONV;
        }

      char *section_name = malloc (RXI_STRING_MAX * sizeof (*section_name));
      sprintf (section_name, "Partner %d", i + 1);
      ini_stat += !ini_puts (section_name, "partner", cp_name, info_filename);
      ini_stat += !ini_putl (section_name, "collisional_transitions",
          mol_info->numof_coll_trans[i], info_filename);
      ini_stat += !ini_putl (section_name, "collisional_temperatures",
          mol_info->numof_coll_temps[i], info_filename);

      char *temps = malloc (RXI_STRING_MAX * sizeof (*temps));
      if (!temps)
        {
          free (cp_name);
          free (info_filename);
          free (section_name);
          return RXI_ERR_ALLOC;
        }

      for (int8_t j = 0; j < mol_info->numof_coll_temps[i]; ++j)
        {
          char *t = malloc (9 * sizeof (*t));
          snprintf (t, 9, "%8.2f ",
                    gsl_matrix_get (mol_info->coll_temps, i, j));
          memcpy (temps + 8 * j, t, 9);
          free (t);
        }
      ini_stat += !ini_puts (section_name, "temperatures", temps,
                             info_filename);

      free (cp_name);
      free (temps);
      free (section_name);
    }

  free (info_filename);

  if (ini_stat != 0)
    return RXI_ERR_FILE;
  else
    return RXI_OK;
}

static RXI_STAT
rxi_add_molecule_csv (FILE *molfile, const char *db_folder,
                      const char *filename, size_t nlines)
{
  char *csv_filename = malloc (RXI_PATH_MAX * sizeof (*csv_filename));
  CHECK (csv_filename && "Allocation error");
  if (!csv_filename)
    return RXI_ERR_ALLOC;

  char *line = malloc (RXI_STRING_MAX * sizeof (*line));
  CHECK (line && "Allocation error");
  if (!line)
    {
      free (csv_filename);
      return RXI_ERR_ALLOC;
    }

  strcpy (csv_filename, db_folder);
  strcat (csv_filename, "/");
  strcat (csv_filename, filename);

  FILE *csv = fopen (csv_filename, "w");
  CHECK (csv && "File open error");
  if (!csv)
    {
      free (csv_filename);
      free (line);
      return RXI_ERR_FILE;
    }

  DEBUG ("Write information to `%s'", csv_filename);

  RXI_STAT stat = RXI_OK;
  for (size_t i = 0; i < nlines; ++i)
    {
      if (!fgets (line, RXI_STRING_MAX, molfile))
        {
          stat = RXI_ERR_FILE;
          break;
        }

      stat = rxi_csv_write_line (csv, line);
      if (stat != RXI_OK)
        break;
    }

  CHECK ((stat == RXI_OK) && "CSV error");

  free (csv_filename);
  free (line);
  fclose (csv);
  return stat;
}

RXI_STAT
rxi_add_molecule (const char *name, const char *path)
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

  // LAMDA database file parsing starts here
  struct rxi_db_molecule_info *mol_info;
  status = rxi_db_molecule_info_malloc (&mol_info);
  CHECK ((status == RXI_OK) && "Allocation error");
  if (status != RXI_OK)
    {
      free (db_folder);
      free ((void*)db_path);
      fclose (molfile);
      return status;
    }

  do
    {
      // Saving header
      status = rxi_save_molecule_info (molfile, mol_info, 0);
      CHECK ((status == RXI_OK) && "Information file error");
      if (status != RXI_OK)
        break;

      // Writing energy levels information to enlev.csv
      status = rxi_add_molecule_csv (molfile, db_folder, "enlev.csv",
                                     mol_info->numof_enlev);
      CHECK ((status == RXI_OK) && "Energy levels file error");
      if (status != RXI_OK)
        break;

      // Saving number of radiation transitions
      status = rxi_save_molecule_info (molfile, mol_info, 0);
      CHECK ((status == RXI_OK) && "Information file error");
      if (status != RXI_OK)
        break;

      // Writing radiation transitions info to radtr.csv
      status = rxi_add_molecule_csv (molfile, db_folder, "radtr.csv",
                                     mol_info->numof_radtr);
      CHECK ((status == RXI_OK) && "Radiative transition file error");
      if (status != RXI_OK)
        break;

      // Saving first collisional partner information
      status = rxi_save_molecule_info (molfile, mol_info, 0);
      CHECK ((status == RXI_OK) && "Information file error");
      if (status != RXI_OK)
        break;

      // Cycle to loop through collisional transitions
      for (int8_t i = 0; i < mol_info->numof_coll_part; ++i)
        {
          char *cp_name = numtoname (mol_info->coll_part[i]);
          if (!cp_name)
            {
              status = RXI_WARN_LAMDA;
              break;
            }
          strcat (cp_name, ".csv");

          // Writing collisional partners information to specified file
          status = rxi_add_molecule_csv (molfile, db_folder, cp_name,
              mol_info->numof_coll_trans[i]);
          CHECK ((status == RXI_OK) && "Colision partner file error");
          if (status != RXI_OK)
            break;

          status = rxi_save_molecule_info (molfile, mol_info, i + 1);
          CHECK ((status == RXI_OK) && "Information file error");
          if (status != RXI_OK)
            break;

          free (cp_name);
        }

      status = rxi_add_molecule_info (db_folder, name, mol_info);
      CHECK ((status == RXI_OK) && "Information file error");
    }
  while (false);

  rxi_db_molecule_info_free (mol_info);
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

static RXI_STAT
rxi_delete_molecule_folder (const char *db_folder, const struct stat *sb)
{
  RXI_STAT status = RXI_OK;
  DIR *dir;
  struct dirent *entry;
  struct stat sb_entry;

  // Check if it is a directory
  if (S_ISDIR (sb->st_mode) == 0)
    return RXI_ERR_FILE;

  // Check if it is openable
  if ((dir = opendir (db_folder)) == NULL)
    return RXI_ERR_FILE;

  while ((entry = readdir (dir)) != NULL)
    {
      if (!strcmp (entry->d_name, ".") || !strcmp (entry->d_name, ".."))
        continue;

      char *entry_pathname = malloc (RXI_STRING_MAX * sizeof (*entry_pathname));
      strcpy (entry_pathname, db_folder);
      strcat (entry_pathname, "/");
      strcat (entry_pathname, entry->d_name);

      stat (entry_pathname, &sb_entry);

      DEBUG ("Removing `%s'", entry_pathname);
      // Recursively delete directories
      if (S_ISDIR (sb_entry.st_mode) != 0)
        {
          status = rxi_delete_molecule_folder (entry_pathname, &sb_entry);
          if (status != RXI_OK)
            {
              free (entry_pathname);
              break;
            }

          free (entry_pathname);
          continue;
        }

      if (unlink (entry_pathname) != 0)
        status = RXI_ERR_FILE;

      free (entry_pathname);
    }

  closedir (dir);

  if (rmdir (db_folder) != 0)
    status = RXI_ERR_FILE;

  return status;
}

RXI_STAT
rxi_delete_molecule (const char *name)
{
  DEBUG ("Start delete molecule '%s'", name);

  RXI_STAT status = RXI_OK;

  const char *db_path = rxi_database_path ();
  CHECK (db_path && "Allocation error");
  if (!db_path)
    return RXI_ERR_ALLOC;

  char *db_folder = malloc (RXI_PATH_MAX * sizeof (*db_folder));
  CHECK (db_folder && "Allocation error");
  if (!db_folder)
    {
      free ((void*)db_path);
      return RXI_ERR_ALLOC;
    }
  strcpy (db_folder, db_path);
  strcat (db_folder, name);

  DEBUG ("Local database folder `%s'", db_folder);

  // Check if local database folder with this name exist
  struct stat sb;
  if (stat (db_folder, &sb) != -1)
    {
      DEBUG ("Folder exists, trying to delete");
      status = rxi_delete_molecule_folder (db_folder, &sb);
    }
  else
    {
      status = RXI_WARN_NOFILE;
    }

  free ((void*)db_path);
  free (db_folder);
  return status;
}

RXI_STAT
rxi_list_molecules ()
{
  return RXI_OK;
}

int
rxi_db_molecule_iter (DIR *dir, char *name)
{
  struct dirent *entry;
  struct stat sb;
  if ((entry = readdir (dir)) != NULL)
    {
      strncpy (name, entry->d_name, 15);
      if (strncmp (entry->d_name, ".", 1) == 0)
        return -1;

      stat (entry->d_name, &sb);
      if (S_ISDIR (sb.st_mode))
        return -2;
    }
  else
    return 0;

  return 1;
}

RXI_STAT
rxi_db_read_molecule_info (const char *name,
                           struct rxi_db_molecule_info *mol_info)
{
  const char *db_path = rxi_database_path ();
  CHECK (db_path && "Allocation error");
  if (!db_path)
    return RXI_ERR_FILE;

  char *filename = malloc (RXI_PATH_MAX * sizeof (*filename));
  CHECK (filename && "Allocation error");
  if (!filename)
    {
      free ((void*)db_path);
      return RXI_ERR_ALLOC;
    }
  strcpy (filename, db_path);
  strcat (filename, name);
  strcat (filename, "/");
  strcat (filename, name);
  strcat (filename, ".info");
  DEBUG ("%s", filename);

  struct stat sb;
  stat (filename, &sb);
  if (S_ISDIR (sb.st_mode))
    {
      free ((void*)db_path);
      free (filename);
      return RXI_ERR_FILE;
    }

  ini_gets ("Information", "name", "no_name", mol_info->name, RXI_STRING_MAX,
            filename);
  DEBUG ("Molecule name: %s", mol_info->name);

  mol_info->weight = ini_getf ("Information", "weight", 0., filename);
  DEBUG ("Molecule weight: %f", mol_info->weight);

  mol_info->numof_enlev = ini_getl ("Information", "energy_levels", 0,
                                    filename);
  DEBUG ("Number of energy levels: %d", mol_info->numof_enlev);

  mol_info->numof_radtr = ini_getl ("Information", "radiative_transitions", 0,
                                    filename);
  DEBUG ("Number of radiative transitions: %d", mol_info->numof_radtr);

  mol_info->numof_coll_part = ini_getl ("Information", "collision_partners",
                                        0, filename);
  DEBUG ("Number of collision partners: %d", mol_info->numof_coll_part);

  for (int8_t i = 0; i < mol_info->numof_coll_part; ++i)
    {
      char *section_name = malloc (RXI_STRING_MAX * sizeof (*section_name));
      CHECK (section_name && "Allocation error");
      if (!section_name)
        {
          free ((void*)db_path);
          free (filename);
          return RXI_ERR_ALLOC;
        }

      char *cp_name = malloc (RXI_STRING_MAX * sizeof (*cp_name));
      CHECK (cp_name && "Allocation error");
      if (!cp_name)
        {
          free ((void*)db_path);
          free (filename);
          free (section_name);
          return RXI_ERR_ALLOC;
        }

      sprintf (section_name, "Partner %d", i + 1);
      ini_gets (section_name, "partner", "no_name", cp_name, RXI_STRING_MAX,
                filename);
      mol_info->coll_part[i] = nametonum (cp_name);
      DEBUG ("%d collision partner name: %s -> %d", i, cp_name,
             mol_info->coll_part[i]);
      free (cp_name);

      mol_info->numof_coll_trans[i] = ini_getl (section_name,
          "collisional_transitions", 0, filename);
      DEBUG ("%d number of collisional transitions: %d", i,
             mol_info->numof_coll_trans[i]);

      mol_info->numof_coll_temps[i] = ini_getl (section_name,
          "collisional_temperatures", 0, filename);
      DEBUG ("%d number of collisional temperatures: %d", i,
             mol_info->numof_coll_temps[i]);

      char *start = malloc (RXI_STRING_MAX * sizeof (*start));
      ini_gets (section_name, "temperatures", "no_temps", start,
                RXI_STRING_MAX, filename);
      char *end;
      int8_t j = 0;
      for (float f = strtof (start, &end);
           start != end;
           f = strtof (start, &end))
        {
          start = end;
          gsl_matrix_set (mol_info->coll_temps, i, j, f);
          DEBUG ("%d temperature: %f", i,
                 gsl_matrix_get (mol_info->coll_temps, i, j));
          ++j;
        }

      free (section_name);
    }

  free ((void*)db_path);
  free (filename);

  return RXI_OK;
}

RXI_STAT
rxi_db_read_molecule_enlev (const char *name,
                            struct rxi_db_molecule_enlev *mol_enl)
{
  char *filename = malloc (RXI_PATH_MAX * sizeof (*filename));
  CHECK (filename && "Allocation error");
  if (!filename)
    return RXI_ERR_ALLOC;

  const char *db_path = rxi_database_path ();
  CHECK (db_path && "Allocation error");
  if (!db_path)
    {
      free (filename);
      return RXI_ERR_ALLOC;
    }

  strcpy (filename, db_path);
  strcat (filename, name);
  strcat (filename, "/enlev.csv");
  free ((void*)db_path);

  DEBUG ("Reading enlev from `%s'", filename);

  FILE *enlev_csv = fopen (filename, "r");
  CHECK (enlev_csv && "Error open file");
  if (!enlev_csv)
    {
      free (filename);
      return RXI_ERR_FILE;
    }

  char *buff[RXI_ELEMENTS_MAX];
  for (size_t i = 0; i < RXI_ELEMENTS_MAX; ++i)
    {
      buff[i] = malloc (RXI_QNUM_MAX * sizeof (buff[i]));
      CHECK (buff[i] && "Allocation error");
      if (!buff[i])
        {
          free (filename);
          fclose (enlev_csv);
          for (size_t j = 0; j < i; ++j)
            free (buff[j]);

          return RXI_ERR_ALLOC;
        }
    }

  int n = 0;
  RXI_STAT stat = RXI_OK;
  for (stat = rxi_csv_read_line (enlev_csv, (void**)buff);
       stat == RXI_OK;
       stat = rxi_csv_read_line (enlev_csv, (void**)buff))
    {
      mol_enl->level[n] = strtol (buff[0], NULL, 10);
      mol_enl->term[n] = strtod (buff[1], NULL);
      mol_enl->weight[n] = strtod (buff[2], NULL);
      ++n;
    }

  fclose (enlev_csv);
  free (filename);
  free (buff[0]);

  DEBUG ("Will return %d", stat);
  if (stat != RXI_FILE_END)
    return stat;

  return RXI_OK;
}

RXI_STAT
rxi_db_read_molecule_radtr (const char *name,
                            struct rxi_db_molecule_radtr *mol_radtr)
{
  char *filename = malloc (RXI_PATH_MAX * sizeof (*filename));
  CHECK (filename && "Allocation error");
  if (!filename)
    return RXI_ERR_ALLOC;

  const char *db_path = rxi_database_path ();
  CHECK (db_path && "Allocation error");
  if (!db_path)
    {
      free (filename);
      return RXI_ERR_ALLOC;
    }

  strcpy (filename, db_path);
  strcat (filename, name);
  strcat (filename, "/radtr.csv");
  free ((void*)db_path);

  DEBUG ("Reading radtr from `%s'", filename);

  FILE *radtr_csv = fopen (filename, "r");
  CHECK (radtr_csv && "Error open file");
  if (!radtr_csv)
    {
      free (filename);
      return RXI_ERR_FILE;
    }

  char *buff[RXI_ELEMENTS_MAX];
  for (size_t i = 0; i < RXI_ELEMENTS_MAX; ++i)
    {
      buff[i] = malloc (RXI_QNUM_MAX * sizeof (buff[i]));
      CHECK (buff[i] && "Allocation error");
      if (!buff[i])
        {
          free (filename);
          fclose (radtr_csv);
          for (size_t j = 0; j < i; ++j)
            free (buff[j]);

          return RXI_ERR_ALLOC;
        }
    }

  int n = 0;
  RXI_STAT stat;
  for (stat = rxi_csv_read_line (radtr_csv, (void**)buff);
       stat != RXI_FILE_END;
       stat = rxi_csv_read_line (radtr_csv, (void**)buff))
    {
      if (stat != RXI_OK)
        break;

      mol_radtr->up[n] = strtol (buff[1], NULL, 10);
      mol_radtr->low[n] = strtol (buff[2], NULL, 10);
      mol_radtr->einst[n] = strtod (buff[3], NULL);
      mol_radtr->freq[n] = strtod (buff[4], NULL);
      mol_radtr->up_en[n] = strtod (buff[5], NULL);
      ++n;
    }

  fclose (radtr_csv);
  free (filename);
  free (buff[0]);

  if (stat != RXI_FILE_END)
    return stat;

  return RXI_OK;
}

RXI_STAT
rxi_db_read_molecule_coll_part (const char *mol_name, const COLL_PART cp,
    const size_t n_temps, struct rxi_db_molecule_coll_part *mol_cp)
{
  char *filename = malloc (RXI_PATH_MAX * sizeof (*filename));
  CHECK (filename && "Allocation error");
  if (!filename)
    return RXI_ERR_ALLOC;

  const char *db_path = rxi_database_path ();
  CHECK (db_path && "Allocation error");
  if (!db_path)
    {
      free (filename);
      return RXI_ERR_ALLOC;
    }

  char *cp_name = numtoname (cp);
  CHECK (cp_name && "Allocation error");
  if (!cp_name)
    {
      free (filename);
      free ((void*)db_path);
      return RXI_ERR_ALLOC;
    }

  strcpy (filename, db_path);
  strcat (filename, mol_name);
  strcat (filename, "/");
  strcat (filename, cp_name);
  strcat (filename, ".csv");
  free ((void*)db_path);
  free (cp_name);

  DEBUG ("Reading collision partner from `%s'", filename);

  FILE *radtr_csv = fopen (filename, "r");
  CHECK (radtr_csv && "Error open file");
  if (!radtr_csv)
    {
      free (filename);
      return RXI_ERR_FILE;
    }

  char *buff[RXI_ELEMENTS_MAX];
  for (size_t i = 0; i < RXI_ELEMENTS_MAX; ++i)
    {
      buff[i] = malloc (RXI_QNUM_MAX * sizeof (buff[i]));
      CHECK (buff[i] && "Allocation error");
      if (!buff[i])
        {
          free (filename);
          fclose (radtr_csv);
          for (size_t j = 0; j < i; ++j)
            free (buff[j]);

          return RXI_ERR_ALLOC;
        }
    }

  int n = 0;
  RXI_STAT stat;
  for (stat = rxi_csv_read_line (radtr_csv, (void**)buff);
       stat != RXI_FILE_END;
       stat = rxi_csv_read_line (radtr_csv, (void**)buff))
    {
      if (stat != RXI_OK)
        break;

      mol_cp->up[n] = strtol (buff[1], NULL, 10);
      mol_cp->low[n] = strtol (buff[2], NULL, 10);
      for (size_t i = 0; i < n_temps; ++i)
        gsl_matrix_set (mol_cp->coll_rates, n, i, strtod (buff[i + 3], NULL));
      ++n;
    }

  fclose (radtr_csv);
  free (filename);
  free (buff[0]);

  if (stat != RXI_FILE_END)
    return stat;

  return RXI_OK;
}
