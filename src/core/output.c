/**
 * @file output.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "rxi_common.h"
#include "core/output.h"
#include "utils/cli_tools.h"
#include "utils/debug.h"

int
rxi_check_path (const char *path)
{
  struct stat sb;
  stat(path, &sb);

  if (S_ISREG (sb.st_mode))
      return 1;
  else if (S_ISDIR (sb.st_mode))
      return 2;

  FILE *check;
  check = fopen (path, "w");
  if (check)
    {
      fclose (check);
      return 0;
    }
  else
    {
      return -1;
    }
}

void
rxi_out_print (struct rxi_calc_data **data,
               const struct rxi_calc_results *results)
{
  DEBUG ("Print resuts to command line");
  char *geometry = geomtoname (data[0]->input.geom);

  printf ("* Radexi version             : " RXI_VERSION "\n");
  printf ("* Geometry                   : %s\n", geometry);
  free (geometry);
  for (int i = 0; i < data[0]->input.numof_molecules; ++i)
    printf ("* Molecule                   : %s\n", data[i]->input.name);
  printf ("* Kinetic temperature    [K] : %.3f\n", data[0]->input.temp_kin);
  printf ("* Background temperature [K] : %.3f\n", data[0]->input.temp_bg);
  printf ("* Column density      [cm-2] : %.3e\n", data[0]->input.col_dens);
  printf ("* Line width          [km/s] : %.3f\n", data[0]->input.line_width);
  for (int8_t i = 0; i < data[0]->input.n_coll_partners; ++i)
    {
      char *cp_name = numtoname (data[0]->input.coll_part[i]);
      printf ("* Density of %9s[cm-3] : %.3e\n", cp_name,
              data[0]->input.coll_part_dens[i]);
      free (cp_name);
    }

  printf ("*    LINE    MOLECULE      E_UP          FREQ         WAVEL        T_EX       \
  TAU         T_R         POP         POP        FLUX       FLUX\n");
  printf ("*                           [K]         [GHz]          [nm]         [K]       \
              [K]          UP         LOW    [K km/s]   [erg cm-2 s-1]\n");

  char line_format[120];
  strcpy (line_format, "%3d -> %3d %10s %8.1f  %12.4f  %12.4f  %10.3f  %10.3e \
 %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %2d\n");

  int size = 0;
  for (int8_t i = 0; i < data[0]->input.numof_molecules; ++i)
      size += data[i]->numof_radtr;

  for (int i = 0; i < size; ++i)
    {
      const int u = results[i].up - 1;
      const int l = results[i].low - 1;
      const float freq = results[i].spfreq;
      int m = 0;
      for (int j = 0; j < size; ++j)
        {
          if (!strcmp (results[i].name, data[j]->input.name))
            {
              m = j;
              break;
            }
        }

      if ((freq < data[0]->input.sfreq) || (freq > data[0]->input.efreq))
        continue;

      int blend_count = 0;
      float line_width_freq = data[0]->input.line_width * freq * 1e5 / RXI_SOL;
      for (int j = 0; j < size; ++j)
        {
          if (i == j)
            continue;

          float overlap = fabs (freq - results[j].spfreq);
          if (overlap <= line_width_freq)
            ++blend_count;

          if ((blend_count > 0) && (overlap > line_width_freq))
            break;
        }

      const double xt = gsl_pow_3 (results[i].xnu);
      printf (line_format, u, l, results[i].name, freq * 1e9 * RXI_HP / RXI_KB, freq,
          RXI_SOL / freq / 1e5, results[i].excit_temp, results[i].tau,
          results[i].antenna_temp, gsl_vector_get (data[m]->pop, u),
          gsl_vector_get (data[m]->pop, l),
          1.0645 * data[0]->input.line_width * results[i].antenna_temp,
          1.0645 * 8 * M_PI * RXI_KB * data[0]->input.line_width
            * results[i].antenna_temp * xt,
          blend_count);
    }
}

void
rxi_out_result_sort (struct rxi_calc_results *results, size_t results_size)
{
  DEBUG ("Sort results by spectral frequencies");
  for (size_t i = 0; i < results_size; ++i)
    {
      size_t e = i;
      for (size_t j = i; j < results_size; ++j)
        {
          if (results[j].spfreq < results[e].spfreq)
            e = j;
        }
      struct rxi_calc_results element = results[e];
      results[e] = results[i];
      results[i] = element;
    }
}

RXI_STAT
rxi_out_result (struct rxi_calc_data *data[RXI_MOLECULE_MAX],
                const struct rxi_options *opts)
{
  int size = 0;
  for (int8_t i = 0; i < data[0]->input.numof_molecules; ++i)
      size += data[i]->numof_radtr;

  struct rxi_calc_results output[size];
  int k = 0;
  for (int8_t i = 0; i < data[0]->input.numof_molecules; ++i)
    {
      for (size_t j = 0; j < data[i]->numof_radtr; ++j)
        {
          output[k].up = data[i]->up[j];
          output[k].low = data[i]->low[j];
          strcpy (output[k].name, data[i]->input.name);
          output[k].spfreq = gsl_matrix_get (data[i]->freq, output[k].up - 1,
                                             output[k].low - 1);
          output[k].xnu = gsl_vector_get (data[i]->term, output[k].up - 1)
                          - gsl_vector_get (data[i]->term, output[k].low - 1);
          output[k].tau = gsl_matrix_get (data[i]->tau, output[k].up - 1,
                                          output[k].low - 1);
          output[k].excit_temp =  gsl_matrix_get (data[i]->excit_temp,
              output[k].up - 1, output[k].low - 1);
          output[k].antenna_temp =  gsl_matrix_get (data[i]->antenna_temp,
              output[k].up - 1, output[k].low - 1);
          ++k;
        }
    }

  rxi_out_result_sort (output, size);

  if (opts->cmd_output)
    rxi_out_print (data, output);

  if (opts->no_result_file)
    return RXI_OK;

  char *path;
  path = malloc (RXI_PATH_MAX);
  CHECK (path);
  if (!path)
    return RXI_ERR_ALLOC;
  strcpy (path, opts->result_path);

  struct stat sb;
  stat (opts->result_path, &sb);
  if (S_ISDIR (sb.st_mode))
    {
      strcat (path, "/");
      strcat (path, "rxi_out.txt");
    }

  stat (path, &sb);
  if (S_ISREG (sb.st_mode) && !opts->force_fs)
    {
      printf ("  ## Specified file `%s` already exists. Rewrite it?\n", path);
      if (!rxi_readline_accept ())
        {
          free (path);
          return RXI_OK;
        }
    }

  DEBUG ("Opening `%s` file for output result", path);
  FILE *result_file;
  result_file = fopen (path, "w");
  if (!result_file)
    {
      free (path);
      return RXI_ERR_FILE;
    }

  char *geometry = geomtoname (data[0]->input.geom);

  fprintf (result_file, "* Radexi version             : " RXI_VERSION "\n");
  fprintf (result_file, "* Geometry                   : %s\n", geometry);
  free (geometry);
  for (int i = 0; i < data[0]->input.numof_molecules; ++i)
    {
      fprintf (result_file, "* Molecule                   : %s\n",
               data[i]->input.name);
    }
  fprintf (result_file, "* Kinetic temperature    [K] : %.3f\n",
           data[0]->input.temp_kin);
  fprintf (result_file, "* Background temperature [K] : %.3f\n",
           data[0]->input.temp_bg);
  fprintf (result_file, "* Column density      [cm-2] : %.3e\n",
           data[0]->input.col_dens);
  fprintf (result_file, "* Line width          [km/s] : %.3f\n",
           data[0]->input.line_width);
  for (int8_t i = 0; i < data[0]->input.n_coll_partners; i++)
    {
      char *cp_name = numtoname (data[0]->input.coll_part[i]);
      fprintf (result_file, "* Density of %9s[cm-3] : %.3e\n", cp_name,
               data[0]->input.coll_part_dens[i]);
      free (cp_name);
    }

  fprintf (result_file, "*    LINE    MOLECULE      E_UP          FREQ         WAVEL        T_EX       \
  TAU         T_R         POP         POP        FLUX       FLUX\n");
  fprintf (result_file, "*                           [K]         [GHz]          [nm]         [K]       \
              [K]          UP         LOW    [K km/s]   [erg cm-2 s-1]\n");

  char line_format[120];
  strcpy (line_format, "%3d -> %3d %10s %8.1f  %12.4f  %12.4f  %10.3f  %10.3e \
 %10.3e  %10.3e  %10.3e  %10.3e  %10.3e  %2d\n");

  for (int i = 0; i < size; ++i)
    {
      const int u = output[i].up - 1;
      const int l = output[i].low - 1;
      const float freq = output[i].spfreq;
      DEBUG ("%f", freq);
      int m = 0;
      for (int j = 0; j < data[0]->input.numof_molecules; ++j)
        {
          if (!strcmp (output[i].name, data[j]->input.name))
            {
              m = j;
              break;
            }
        }

      if ((freq < data[0]->input.sfreq) || (freq > data[0]->input.efreq))
        continue;

      int blend_count = 0;
      float line_width_freq = data[0]->input.line_width * freq * 1e5 / RXI_SOL;
      for (int j = 0; j < size; ++j)
        {
          if (i == j)
            continue;

          float overlap = fabs (freq - output[j].spfreq);
          if (overlap <= line_width_freq)
            ++blend_count;

          if ((blend_count > 0) && (overlap > line_width_freq))
            break;
        }

      const double xt = gsl_pow_3 (output[i].xnu);
      fprintf (result_file, line_format, u, l, output[i].name, freq * 1e9 * RXI_HP / RXI_KB, freq,
          RXI_SOL / freq / 1e5, output[i].excit_temp, output[i].tau,
          output[i].antenna_temp, gsl_vector_get (data[m]->pop, u),
          gsl_vector_get (data[m]->pop, l),
          1.0645 * data[0]->input.line_width * output[i].antenna_temp,
          1.0645 * 8 * M_PI * RXI_KB * data[0]->input.line_width
            * output[i].antenna_temp * xt,
          blend_count);
    }

  DEBUG ("Finish");
  fclose (result_file);
  free (path);
  return RXI_OK;
}
