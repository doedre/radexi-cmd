/**
 * @file core/dialogue.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <dirent.h>

#include "core/dialogue.h"

#include "utils/cli_tools.h"
#include "utils/database.h"
#include "utils/debug.h"

static void
print_dialog_error ()
{
  printf ("ERROR\n");
}

static RXI_STAT
get_molecule_name (char *names, char name[10][15], int8_t *numof_molecules)
{
  DEBUG ("Get molecule name");

  char *line = malloc (RXI_MOLECULE_MAX * sizeof (*line));
  CHECK (line && "Allocation error");
  if (!line)
    return RXI_ERR_ALLOC;

  char *mname = malloc (RXI_MOLECULE_MAX * sizeof (*mname));
  CHECK (mname && "Allocation error");
  if (!mname)
    {
      free (line);
      return RXI_ERR_ALLOC;
    }

  RXI_STAT status = RXI_OK;
  status = rxi_history_load ("mname.history");
  CHECK ((status == RXI_OK) && "Can't load history file");
  if (status != RXI_OK)
    {
      free (line);
      free (mname);
      return status;
    }

  *numof_molecules = 0;
  bool is_written = true;
  while ((line = rxi_readline ("  >> ")) != NULL)
    {
      const char *db_path = rxi_database_path ();

      int n = sscanf (line, "%s %s %s %s %s %s %s %s %s %s",
                      name[0], name[1], name[2], name[3], name[4], name[5],
                      name[6], name[7], name[8], name[9]);

      is_written = true;
      for (int i = 0; (i < n) && is_written; ++i)
        {
          DIR *dir;
          if ((dir = opendir (db_path)) == NULL)
            {
              free (line);
              free (mname);
              free ((void*)db_path);
              return RXI_ERR_FILE;
            }

          is_written = false;
          for (int p = rxi_db_molecule_iter (dir, mname); p;
               p = rxi_db_molecule_iter (dir, mname))
            {
              if ((p > 0) && !strncmp (name[i], mname, 15))
                {
                  is_written = true;
                  break;
                }
            }
          free (dir);

          if (is_written == false)
            break;
        }
      free ((void*)db_path);

      if (!is_written)
        {
          print_dialog_error ();
          continue;
        }

      rxi_history_save (line, "mname.history");
      *numof_molecules = n;
      strcpy(names, line);
      break;
    }
  free (mname);
  free (line);

  DEBUG ("Number of molecules: %d", *numof_molecules);

  return status;
}

static RXI_STAT
get_frequencies (float *sfreq, float *efreq)
{
  DEBUG ("Get frequencies");

  char *line = malloc (RXI_STRING_MAX * sizeof (*line));
  CHECK (line && "Allocation error");
  if (!line)
    return RXI_ERR_ALLOC;

  RXI_STAT status = rxi_history_load ("freq.history");
  CHECK ((status == RXI_OK) && "Can't load history file");

  status = RXI_OK;
  while ((line = rxi_readline ("  >> ")) != NULL)
    {
      int n = sscanf (line, "%f %f", sfreq, efreq);
      if ((n < 3) && (*sfreq < *efreq))
        {
          rxi_history_save (line, "freq.history");
          break;
        }
      else
        {
          print_dialog_error ();
        }
    }

  free (line);
  return status;
}

static RXI_STAT
get_kin_temp (double *kin_temp)
{
  DEBUG ("Get kinetic temperature");

  char *line = malloc (RXI_STRING_MAX * sizeof (*line));
  CHECK (line && "Allocation error");
  if (!line)
    return RXI_ERR_ALLOC;

  RXI_STAT status = rxi_history_load ("kin_temp.history");
  CHECK ((status == RXI_OK) && "Can't load history file");

  status = RXI_OK;
  bool is_written = false;
  while (!is_written && ((line = rxi_readline ("  >> ")) != NULL))
    {
      *kin_temp = strtod (line, NULL);
      if ((*kin_temp > 0) && (*kin_temp < 10000))
        {
          is_written = true;
          rxi_history_save (line, "kin_temp.history");
        }
      else
        {
          print_dialog_error ();
        }
    }

  free (line);
  return status;
}

static RXI_STAT
get_bg_temp (double *bg_temp)
{
  DEBUG ("Get background temperature");

  char *line = malloc (RXI_STRING_MAX * sizeof (*line));
  CHECK (line && "Allocation error");
  if (!line)
    return RXI_ERR_ALLOC;

  RXI_STAT status = rxi_history_load ("bg_temp.history");
  CHECK ((status == RXI_OK) && "Can't load history file");

  status = RXI_OK;
  bool is_written = false;
  while (!is_written && ((line = rxi_readline ("  >> ")) != NULL))
    {
      *bg_temp = strtod (line, NULL);
      if ((*bg_temp > 0) && (*bg_temp < 10000))
        {
          is_written = true;
          rxi_history_save (line, "bg_temp.history");
        }
      else
        {
          print_dialog_error ();
        }
    }

  free (line);
  return status;
}

static RXI_STAT
get_coldens (double *coldens)
{
  DEBUG ("Get collision densities");

  char *line = malloc (RXI_STRING_MAX * sizeof (*line));
  CHECK (line && "Allocation error");
  if (!line)
    return RXI_ERR_ALLOC;

  RXI_STAT status = rxi_history_load ("coldens.history");
  CHECK ((status == RXI_OK) && "Can't load history file");

  status = RXI_OK;
  bool is_written = false;
  while (!is_written && ((line = rxi_readline ("  >> ")) != NULL))
    {
      *coldens = strtod (line, NULL);
      if ((*coldens > 0) && (*coldens < 1e25))
        {
          is_written = true;
          rxi_history_save (line, "coldens.history");
        }
      else
        {
          print_dialog_error ();
        }
    }

  free (line);
  return status;
}

static RXI_STAT
get_line_width (double *line_width)
{
  DEBUG ("Get collision densities");

  char *line = malloc (RXI_STRING_MAX * sizeof (*line));
  CHECK (line && "Allocation error");
  if (!line)
    return RXI_ERR_ALLOC;

  RXI_STAT status = rxi_history_load ("line_width.history");
  CHECK ((status == RXI_OK) && "Can't load history file");

  status = RXI_OK;
  bool is_written = false;
  while (!is_written && ((line = rxi_readline ("  >> ")) != NULL))
    {
      *line_width = strtod (line, NULL);
      if ((*line_width > 1e-3) && (*line_width < 1e3))
        {
          is_written = true;
          rxi_history_save (line, "line_width.history");
        }
      else
        {
          print_dialog_error ();
        }
    }

  free (line);
  return status;
}

static RXI_STAT
get_geometry (GEOMETRY *geom)
{
  DEBUG ("Get geometry");

  char *line = malloc (RXI_STRING_MAX * sizeof (*line));
  CHECK (line && "Allocation error");
  if (!line)
    return RXI_ERR_ALLOC;

  RXI_STAT status = rxi_history_load ("geometry.history");
  CHECK ((status == RXI_OK) && "Can't load history file");

  status = RXI_OK;
  bool is_written = false;
  while (!is_written && ((line = rxi_readline ("  >> ")) != NULL))
    {
      if (!strcmp (line, "slab"))
        *geom = SLAB;
      else if (!strcmp (line, "sphere"))
        *geom = SPHERE;
      else if (!strcmp (line, "lvg"))
        *geom = LVG;
      else
        *geom = OTHER;

      if (*geom != OTHER)
        {
          is_written = true;
          rxi_history_save (line, "geometry.history");
        }
      else
        {
          print_dialog_error ();
        }
    }

  free (line);
  return status;
}

static void
parse_collision_partners (char *line, COLL_PART *coll_part,
                          double *coll_part_dens, int8_t *n_coll_part)
{
  int8_t i = 0;
  for (char *tok = strtok (line, ";"); tok; tok = strtok (NULL, ";"))
    {
      DEBUG ("Parsing %s", tok);
      char *pair = malloc (RXI_STRING_MAX * sizeof (*pair));
      CHECK (pair);
      memcpy (pair, tok, RXI_STRING_MAX);

      int n = sscanf (tok, "%s %lf", pair, &coll_part_dens[i]);
      COLL_PART cp = nametonum (pair);
      if ((cp != NO_PARTNER) && (n < 3))
        {
          coll_part[i] = cp;
          ++i;
        }
    }
  *n_coll_part = i;
}

static bool
check_coll_partner (COLL_PART cp, struct rxi_db_molecule_info **mol_info,
                    int numof_molecules)
{
  bool result = false;
  for (int8_t j = 0; j < numof_molecules; ++j)
    {
      result = false;
      for (int8_t i = 0; i < mol_info[j]->numof_coll_part; ++i)
        {
          DEBUG ("Comparing collision partners %u and %u: ", cp,
                 mol_info[j]->coll_part[i]);
          if (cp == mol_info[j]->coll_part[i])
            {
              result = true;
              break;
            }
        }

      if (result == false)
        break;
    }
    return result;
}

static RXI_STAT
get_collision_partners (struct rxi_input_data *inp_data)
{
  DEBUG ("Get collision partners");

  char *line = malloc (RXI_STRING_MAX * sizeof (*line));
  CHECK (line && "Allocation error");
  if (!line)
    return RXI_ERR_ALLOC;

  RXI_STAT status = rxi_history_load ("coll_part.history");
  CHECK ((status == RXI_OK) && "Can't load history file");

  status = RXI_OK;
  struct rxi_db_molecule_info *mol_info[inp_data->numof_molecules];
  for (int i = 0; i < inp_data->numof_molecules; ++i)
    {
      status = rxi_db_molecule_info_malloc (&mol_info[i]);
      CHECK ((status == RXI_OK) && "Allocation error");
      if (status != RXI_OK)
        {
          free (line);
          return status;
        }

      status = rxi_db_read_molecule_info (inp_data->name_list[i], mol_info[i]);
      CHECK ((status == RXI_OK) && "Info file error");
    }

  bool is_written = false;
  while (!is_written && ((line = rxi_readline ("  >> ")) != NULL))
    {
      parse_collision_partners (line, inp_data->coll_part,
          inp_data->coll_part_dens, &inp_data->n_coll_partners);

      for (int8_t i = 0; i < inp_data->n_coll_partners; ++i)
        {
          DEBUG ("Checking collision partners");
          if (!check_coll_partner (inp_data->coll_part[i], mol_info,
                                   inp_data->numof_molecules))
            {
              print_dialog_error ();
              is_written = false;
              break;
            }
          else
            {
              is_written = true;
            }
        }
    }

  rxi_history_save (line, "coll_part.history");

  for (int i = 0; i < inp_data->numof_molecules; ++i)
    rxi_db_molecule_info_free (mol_info[i]);

  free (line);
  return status;
}

RXI_STAT
rxi_dialog_input (struct rxi_input_data *inp_data,
                  const struct rxi_options *opts)
{
  DEBUG ("Begin dialog with user");

  RXI_STAT status = RXI_OK;
  if (!opts->quite_start)
    {
      printf ("STARTING INFO\n");
    }

  printf ("  ## Enter molecule name\n");

  status = get_molecule_name (inp_data->names, inp_data->name_list,
                              &inp_data->numof_molecules);
  CHECK ((status == RXI_OK) && "Error getting molecule name");

  printf ("  ## Enter frequencies\n");

  status = get_frequencies (&inp_data->sfreq, &inp_data->efreq);
  CHECK ((status == RXI_OK) && "Error getting freqs");

  printf ("  ## Enter kinetic temperature\n");

  status = get_kin_temp (&inp_data->temp_kin);
  CHECK ((status == RXI_OK) && "Error getting kinetic temperature");

  printf ("  ## Enter background temperature\n");

  status = get_bg_temp (&inp_data->temp_bg);
  CHECK ((status == RXI_OK) && "Error getting background temperature");

  printf ("  ## Enter column density\n");

  status = get_coldens (&inp_data->col_dens);
  CHECK ((status == RXI_OK) && "Error getting column density");

  printf ("  ## Enter line width\n");

  status = get_line_width (&inp_data->line_width);
  CHECK ((status == RXI_OK) && "Error getting line width");

  printf ("  ## Enter geometry\n");

  status = get_geometry (&inp_data->geom);
  CHECK ((status == RXI_OK) && "Error getting geometry");

  printf ("  ## Enter collision partners and their densities");

  status = get_collision_partners (inp_data);
  CHECK ((status == RXI_OK) && "Error getting collision partners");

  return status;
}
