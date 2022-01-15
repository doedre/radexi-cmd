/**
 * @file rxi_common.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

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

RXI_STAT
rxi_db_molecule_enlev_malloc (struct rxi_db_molecule_enlev **mol_enl,
                              const size_t n_enlev)
{
  DEBUG ("Allocating memory for enlev");
  struct rxi_db_molecule_enlev *me = malloc (sizeof (*me));
  CHECK (me && "Allocation error");
  if (!me)
    goto malloc_error;

  char **qnum = malloc (n_enlev * sizeof (*qnum));
  CHECK (*qnum && "Allocation error");
  if (!*qnum)
    {
      free (me);
      goto malloc_error;
    }
  for (size_t i = 0; i < n_enlev; ++i)
    {
      qnum[i] = malloc (RXI_QNUM_MAX * sizeof (qnum[i]));
      CHECK (qnum[i] && "Allocation error");
      if (!qnum[i])
        {
          free (me);
          free (qnum);
          for (size_t j = 0; j < i; ++j)
            free (qnum[j]);
          goto malloc_error;
        }
    }

  int *level = malloc (n_enlev * sizeof (*level));
  CHECK (level && "Allocation error");
  if (!level)
    {
      free (me);
      goto malloc_error;
    }

  gsl_vector *energy = gsl_vector_calloc (n_enlev);
  CHECK (energy && "Allocation error");
  if (!energy)
    {
      free (me);
      free (level);
      goto malloc_error;
    }

  gsl_vector *weight = gsl_vector_calloc (n_enlev);
  CHECK (weight && "Allocation error");
  if (!weight)
    {
      free (me);
      free (level);
      gsl_vector_free (energy);
      goto malloc_error;
    }

  me->level = level;
  me->energy = energy;
  me->weight = weight;
  me->qnum = qnum;

  *mol_enl = me;

  return RXI_OK;

malloc_error:
  *mol_enl = NULL;
  return RXI_ERR_ALLOC;
}

void
rxi_db_molecule_enlev_free (struct rxi_db_molecule_enlev *mol_enl)
{
  DEBUG ("Free memory for enlev");
  free (mol_enl->level);
  gsl_vector_free (mol_enl->energy);
//  gsl_vector_free (mol_enl->weight);
  free (mol_enl->qnum);
  free (mol_enl);
}

RXI_STAT
rxi_db_molecule_radtr_malloc (struct rxi_db_molecule_radtr **mol_rat,
                              const size_t n_radtr)
{
  DEBUG ("Allocating memory for radtr");
  struct rxi_db_molecule_radtr *mr = malloc (sizeof (*mr));
  CHECK (mr && "Allocation error");
  if (!mr)
    goto malloc_error;

  int *up = malloc (n_radtr * sizeof (*up));
  CHECK (up && "Allocation error");
  if (!up)
    {
      free (mr);
      goto malloc_error;
    }

  int *low = malloc (n_radtr * sizeof (*low));
  CHECK (low && "Allocation error");
  if (!low)
    {
      free (mr);
      free (up);
      goto malloc_error;
    }

  gsl_vector *einst = gsl_vector_calloc (n_radtr);
  CHECK (einst && "Allocation error");
  if (!einst)
    {
      free (mr);
      free (up);
      free (low);
      goto malloc_error;
    }

  gsl_vector *freq = gsl_vector_calloc (n_radtr);
  CHECK (freq && "Allocation error");
  if (!freq)
    {
      free (mr);
      free (up);
      free (low);
      gsl_vector_free (einst);
      goto malloc_error;
    }

  gsl_vector *up_en = gsl_vector_calloc (n_radtr);
  CHECK (up_en && "Allocation error");
  if (!up_en)
    {
      free (mr);
      free (up);
      free (low);
      gsl_vector_free (einst);
      gsl_vector_free (freq);
      goto malloc_error;
    }

  mr->up = up;
  mr->low = low;
  mr->einst = einst;
  mr->freq = freq;
  mr->up_en = up_en;

  *mol_rat = mr;
  return RXI_OK;

malloc_error:
  *mol_rat = NULL;
  return RXI_ERR_ALLOC;
}

void
rxi_db_molecule_radtr_free (struct rxi_db_molecule_radtr *mol_rat)
{
  DEBUG ("Free memory for radtr");

  free (mol_rat->up);
  free (mol_rat->low);
  gsl_vector_free (mol_rat->einst);
  gsl_vector_free (mol_rat->freq);
//  gsl_vector_free (mol_rat->up_en);
  free (mol_rat);
}

RXI_STAT
rxi_db_molecule_coll_part_malloc (struct rxi_db_molecule_coll_part **mol_cp,
    const size_t n_cp_trans, const size_t n_temps)
{
  DEBUG ("Allocating memory for colision partner");
  struct rxi_db_molecule_coll_part *mp = malloc (sizeof (*mp));
  CHECK (mp && "Allocation error");
  if (!mp)
    goto malloc_error;

  int *up = malloc (n_cp_trans * sizeof (*up));
  CHECK (up && "Allocation error");
  if (!up)
    {
      free (mp);
      goto malloc_error;
    }

  int *low = malloc (n_cp_trans * sizeof (*low));
  CHECK (low && "Allocation error");
  if (!low)
    {
      free (mp);
      free (up);
      goto malloc_error;
    }

  gsl_matrix *coll_rates = gsl_matrix_calloc (n_cp_trans, n_temps);
  CHECK (coll_rates && "Allocation error");
  if (!coll_rates)
    {
      free (mp);
      free (up);
      free (low);
      goto malloc_error;
    }

  mp->up = up;
  mp->low = low;
  mp->coll_rates = coll_rates;

  *mol_cp = mp;
  return RXI_OK;

malloc_error:
  *mol_cp = NULL;
  return RXI_ERR_ALLOC;
}

void
rxi_db_molecule_coll_part_free (struct rxi_db_molecule_coll_part *mol_cp)
{
  DEBUG ("Free memory for collision partner");

  free (mol_cp->up);
  free (mol_cp->low);
  gsl_matrix_free (mol_cp->coll_rates);
  free (mol_cp);
}
