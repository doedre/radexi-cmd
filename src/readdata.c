/*	readdata.c
 *
 * Used to control reading of database for the specified molecule. 
 *
 * ----------------------------------------------------------------------
 *
 * This file is part of the radexi-cmd program - C language implementation of 
 * RADEX software package. It is used to calculate molecular exitation and 
 * radiative transfer in a homogenuous medium. The purpose of this full 
 * refactoring is to speed up the process of calculation and implement some new
 * features like databases. It will make the program more stable and usable, 
 * which is crucial for newcomers.
 *
 * Documentation for radexi-cmd program is posted at 
 * https://github.com/doedre/radexi-cmd
 *
 * ----------------------------------------------------------------------
 *
 * From the RADEX software package:
 *    Documentation for the RADEX program is posted at
 *	  https://personal.sron.nl/~vdtak/radex/index.shtml 
 *
 *	  Although this program has been thoroughly tested, the authors do not 
 *	  claim that it is free of errors and gives correct results in all 
 *	  situations.
 *	
 *	  Publications using this program should make a reference to the paper:
 *	  A&A 468, 627 (2007).
 *
 * ---------------------------------------------------------------------*/

#include "radexi.h"

#include <string.h>
#include <math.h>

/* Reads .info file.
 * Returns:
 * - 0 on success
 * - -1 on .info file opening error
 * - 1 on wrong numbers in parameters  */
int 
read_info_file (struct rxi_input *inp)
{
  char *mi_path;
  mi_path = (char *) malloc (2*RXI_MOLECULE_MAX_SIZE + 12);
  sprintf (mi_path, "data/%s/%s.info", inp->name, inp->name);
  FILE *mi = fopen (mi_path, "r");
  free (mi_path);

  if (!mi)
    return 1;

  char *line;
  size_t n = 200;
  line = (char *) malloc (n);
  enum ColPart col_partner = NO_MOLECULE;
  while (fgets (line, n, mi))
    {
      char *p_ddots = strrchr (line, ':');
      char *parameter = p_ddots ? p_ddots + 1 : line;
      if (strstr (line, "Molecular weight:"))
        inp->weight = strtof (parameter, NULL);
      else if (strstr (line, "Number of energy levels:"))
        inp->numof_enlev = strtof (parameter, NULL);
      else if (strstr (line, "Number of radiative transitions:"))
        inp->numof_radtr = strtof (parameter, NULL);
      else if (strstr (line, "Number of collision partners:"))
        inp->max_numof_colpart = strtof (parameter, NULL);
      else if (strstr (line, "Collision partner:"))
        col_partner = strtol (parameter, NULL, 10);
      else if (strstr (line, "Number of collisional transitions:"))
        {
          for (int i = 0; i < RXI_MAX_COLL_PARTNERS; i++)
            {
              if (inp->cps[i].name != col_partner)
                continue;
              inp->cps[i].numof_coltr = strtof (parameter, NULL);
            }
        } 
      else if (strstr (line, "Collisional temperatures:"))
        {
          for (int i = 0; i < RXI_MAX_COLL_PARTNERS; i++)
            {
              if (inp->cps[i].name != col_partner)
                continue;
              char *end;
              int j = 0;
              for (float f = strtof (parameter, &end); 
                   parameter != end;
                   f = strtof (parameter, &end))
                {
                  inp->cps[i].temps[j++] = f;
                  parameter = end;
                }
              inp->cps[i].numof_temps = j;
            }
        }
    }
  free (line);

  fclose (mi);
  return 0;
}

struct rxi_data*
rxi_data_calloc (const struct rxi_input *inp)
{
  struct rxi_data *result;
  result = malloc (sizeof (struct rxi_data));

  /* Write all static info  */
  strncpy (result->mi.name, inp->name, RXI_MOLECULE_MAX_SIZE);
  result->mc.Tkin           = inp->Tkin;
  result->mc.Tbg            = inp->Tbg;
  result->mc.coldens        = inp->coldens;
  result->mc.geom           = inp->g;
  result->mc.line_width     = inp->fwhm;
  result->mc.numof_colpart  = inp->numof_colpart;
  result->mi.numof_enlev    = inp->numof_enlev;
  result->mi.numof_radtr    = inp->numof_radtr;
  result->mi.weight         = inp->weight;

  /* Allocate memory for collision partners first */
  result->mc.cp = malloc (inp->numof_colpart * sizeof (struct colpart));
  for (unsigned int i = 0; i < inp->numof_colpart; i++)
    {
      result->mc.cp[i].ucollev = malloc (inp->cps[i].numof_coltr *  \
                                         sizeof (unsigned int));
      result->mc.cp[i].lcollev = malloc (inp->cps[i].numof_coltr *  \
                                         sizeof (unsigned int));
      result->mc.cp[i].K = gsl_matrix_calloc (inp->numof_enlev, 
                                              inp->numof_enlev);
      
      result->mc.cp[i].name = inp->cps[i].name;
      result->mc.cp[i].dens = inp->cps[i].dens;
      result->mc.cp[i].numof_coltr = inp->cps[i].numof_coltr;
      result->mc.cp[i].numof_temps = inp->cps[i].numof_temps;
      for (unsigned int t = 0; t < inp->cps[i].numof_temps; t++)
        result->mc.cp[i].temps[t] = inp->cps[i].temps[t];
    }

  /* Then for molecular information in mol_info struct  */
  result->mi.el = malloc (inp->numof_enlev * sizeof (struct enlev));
  result->mi.rt = malloc (inp->numof_radtr * sizeof (struct radtr));
  result->mi.C = gsl_matrix_calloc (inp->numof_enlev, inp->numof_enlev);
  result->mi.Ctot = gsl_vector_calloc (inp->numof_enlev);

  /* For background field */
  result->bg.intens = malloc (inp->numof_radtr * sizeof (double));

  /* For results  */
  result->res.rates = gsl_matrix_calloc (inp->numof_enlev, inp->numof_enlev);
  result->res.Tex = malloc (inp->numof_radtr * sizeof (double));
  result->res.tau = malloc (inp->numof_radtr * sizeof (double));
  result->res.pop = malloc (inp->numof_enlev * sizeof (double));

  return result;
}

void
rxi_data_free (struct rxi_data *rxi)
{
  gsl_vector_free (rxi->mi.Ctot);
  gsl_matrix_free (rxi->mi.C);
  free (rxi->mi.rt);
  free (rxi->mi.el);
  for (unsigned int i = 0; i < rxi->mc.numof_colpart; i++)
    {
      free (rxi->mc.cp[i].ucollev);
      free (rxi->mc.cp[i].lcollev);
      gsl_matrix_free (rxi->mc.cp[i].K);
    }
  free (rxi->bg.intens);
  free (rxi->mc.cp);

  gsl_matrix_free (rxi->res.rates);
  free (rxi->res.Tex);
  free (rxi->res.tau);
  free (rxi->res.pop);

  free (rxi);
}

static int 
reading_enlev (struct rxi_data *rxi)
{
  char *el_path;
  el_path = (char *) malloc (RXI_MOLECULE_MAX_SIZE + 16);
  sprintf (el_path, "data/%s/enlev.csv", rxi->mi.name);
  FILE *enlev_file = fopen (el_path, "r");
  free (el_path);

  if (!enlev_file)
    return 1;

  char *line;
  size_t n = 250;
  line = (char *) malloc (n);
  unsigned int i = 0;
  while (fgets (line, n, enlev_file))
    {
      int t = 0;
      for (char *token = strtok (line, ","); token; token = strtok (NULL, ","))
        {
          if (t == 1)
            {
              rxi->energy_level[i].term = strtod (token, NULL);
            }
          else if (t == 2)
            {
              rxi->energy_level[i].statw = strtof (token, NULL);
            }
          else if (t > 2)
            {
              strcat (rxi->energy_level[i].qnums, token);
            }
          t++;
        }
      i++;
    }
  fclose (enlev_file);
  /* Compare with database's data on number of energy levels  */
  if (i != rxi->n_el)
    return -1;

  return 0;
}

static int 
reading_radtr (struct rxi_data *rxi)
{
  char *rt_path;
  rt_path = (char *) malloc (RXI_MOLECULE_MAX_SIZE + 16);
  sprintf (rt_path, "data/%s/radtr.csv", rxi->mi.name);
  FILE *radtr_file = fopen (rt_path, "r");
  free (rt_path);

  if (!radtr_file)
    return 1;

  char *line;
  size_t n = 250;
  line = (char *) malloc (n);
  unsigned int rtnum = 0;
  while (fgets (line, n, radtr_file))
    {
      int t = 0;
      unsigned int u = 0;
      unsigned int l = 0;
      for (char *token = strtok (line, ","); token; token = strtok (NULL, ","))
        {
          /* May deprecate ulev & llev later, if they will have no other need */
          if (t == 1)
            {
              u = atoi (token) - 1;
              rxi->rad_transfer[rtnum].ulev = u;
            }
          else if (t == 2)
            {
              l = atoi (token) - 1;
              rxi->rad_transfer[rtnum].llev = l;
            }
          else if (t == 3)
            {
              rxi->rad_transfer[rtnum].a_einst = strtod (token, NULL);
            }
          else if (t == 4)
            {
              rxi->rad_transfer[rtnum].spfreq = strtof (token, NULL);
            }
          else if (t == 5)
            {
              rxi->rad_transfer[rtnum].enup = strtof (token, NULL);
            }
          else if (t > 5)
            {
              strcat (rxi->rad_transfer[rtnum].tail, token);
            }
          t++;
        }
      rxi->rad_transfer[rtnum].xnu = rxi->energy_level[u].term - \
                                     rxi->energy_level[l].term;
      rtnum++;
    }
  fclose (radtr_file);
  /* Check final number of radiative transitions with the database's one  */
  if (rtnum != rxi->n_rt)
    return -1;

  return 0;
}

static int 
reading_colpart (struct rxi_data *rxi)
{
  for (unsigned int icp = 0; icp < rxi->mc.numof_colpart; icp++) 
    {
      char *cp_path = (char *) malloc (RXI_MOLECULE_MAX_SIZE + 20);
      char *cp_name = (char *) malloc (RXI_MOLECULE_MAX_SIZE);

      conv_int_to_name (rxi->coll_partner[icp].name, cp_name);
      sprintf (cp_path, "data/%s/%s.csv", rxi->mi.name, cp_name);
      FILE *cp_file = fopen (cp_path, "r");
      free (cp_path);
      free (cp_name);

      if (!cp_file)
        return 1;

      char *line;
      size_t n = 250;
      line = (char *) malloc (n);
      unsigned int icoef = 0;
      while (fgets (line, n, cp_file))
        {
          /* Need to make linear approximation for collision rates */
          int t = 0;
          double lcoef, ucoef = 0;
          double ltemp, utemp = 0;
          unsigned int l, u = 0;
          for ( char *token = strtok (line, ","); 
                token;
                token = strtok (NULL, ","))
            {
              if (t == 1)
                {
                  u = atoi (token) - 1;
                }
              else if (t == 2)
                {
                  l = atoi (token) - 1;
                }
              else if (t > 2)
                {
                  if (rxi->coll_partner[icp].temps[t-3] <= rxi->mc.Tkin)
                    {
                      lcoef = strtod (token, NULL);
                      ltemp = rxi->coll_partner[icp].temps[t-3];
                    }
                  else 
                    {
                      ucoef = strtod (token, NULL);
                      utemp = rxi->coll_partner[icp].temps[t-3];
                      break;
                    }
                }
              t++;
            }
          /* Setting collision rates matrice's elements  */
          if (ltemp > utemp)
            gsl_matrix_set (rxi->coll_partner[icp].K, u, l, lcoef);
          else if (t == 3)
            gsl_matrix_set (rxi->coll_partner[icp].K, u, l, ucoef);
          else 
            {
              double coef = lcoef + (ucoef - lcoef) / (utemp - ltemp) * \
                                    (rxi->mc.Tkin - ltemp);
              gsl_matrix_set (rxi->coll_partner[icp].K, u, l, coef);
            }
          icoef++;
        }
      fclose (cp_file);
      /* Compare if results align with the database one's */
      if (icoef != rxi->coll_partner[icp].numof_coltr)
        return -1;
    }
  return 0;
}

static int 
prepare_for_calculation (struct rxi_data *rxi)
{
  printf ("----making collision rate matrix\n");
  for (unsigned int cp = 0; cp < rxi->mc.numof_colpart; cp++)
    {
      rxi->mc.total_density += rxi->coll_partner[cp].dens;
/*
      for (unsigned int i = 1; i <= rxi->n_el; i++)
        for (unsigned int j = 1; j <= rxi->n_el; j++)
          {
            double crate = rxi->coll_partner[cp].dens * \
                           gsl_matrix_get (rxi->coll_partner[cp].K, i, j);

            gsl_matrix_set (rxi->mi.C, i, j, crate);
          }
          */
      gsl_matrix_scale (rxi->coll_partner[cp].K, rxi->coll_partner[cp].dens);
      gsl_matrix_add (rxi->mi.C, rxi->coll_partner[cp].K);
    }

  printf ("----recalculating it's coefficients by statistical balance\n");
  for (unsigned int i = 0; i < rxi->n_el; i++)
    {
      for (unsigned int j = 0; j < rxi->n_el; j++)
        {
          double ediff = rxi->energy_level[i].term -  \
                         rxi->energy_level[j].term;

          if ((ediff > 0) && (fk * ediff / rxi->mc.Tkin < 160))
            {
              double crate = rxi->energy_level[i].statw /            \
                             rxi->energy_level[j].statw *            \
                             exp (- fk * ediff / rxi->mc.Tkin) *     \
                             gsl_matrix_get (rxi->mi.C, i, j);
              gsl_matrix_set (rxi->mi.C, i, j, crate);
            }
          else 
            gsl_matrix_set (rxi->mi.C, i, j, 0);
        }
    }
  
  printf ("----conventing collisional rate matrix(not matching with radex but should)\n");
  for (unsigned int i = 0; i < rxi->n_el; i++)
    {
      gsl_vector *conv = gsl_vector_alloc (rxi->n_el);
      gsl_matrix_get_row (conv, rxi->mi.C, i);
      gsl_vector_add (rxi->mi.Ctot, conv);
      gsl_vector_free (conv);
    }

  return 0;
}

int
read_data (struct rxi_data *rxi)
{
  printf ("--reading energy levels...\n");
  reading_enlev (rxi);
  printf ("--reading radiative transitions...\n");
  reading_radtr (rxi);
  printf ("--reading collision partners...\n");
  reading_colpart (rxi);
  printf ("--preparing for calculations...\n");
  prepare_for_calculation (rxi);
  return 0;
}
