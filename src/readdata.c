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
static int 
reading_info (struct radexi_data *rxi)
{
  char *mi_path;
  mi_path = (char *) malloc (2*RXI_MOLECULE_MAX_SIZE + 12);
  sprintf (mi_path, "data/%s/%s.info", rxi->mi.name, rxi->mi.name);
  FILE *mi = fopen (mi_path, "r");
  free (mi_path);

  if (!mi)
    return 1;

  char *line;
  size_t n = 200;
  line = (char *) malloc (n);
  enum ColPartner col_partner = NO_MOLECULE;
  while (fgets (line, n, mi))
    {
      char *p_ddots = strrchr (line, ':');
      char *parameter = p_ddots ? p_ddots + 1 : line;
      if (strstr (line, "Molecular weight:"))
        rxi->mi.weight = strtof (parameter, NULL);
      else if (strstr (line, "Number of energy levels:"))
        rxi->mi.numof_enlev = strtof (parameter, NULL);
      else if (strstr (line, "Number of radiative transitions:"))
        rxi->mi.numof_radtr = strtof (parameter, NULL);
      else if (strstr (line, "Number of collision partners:"))
        rxi->mi.numof_colpart = strtof (parameter, NULL);
      else if (strstr (line, "Collision partner:"))
        col_partner = strtol (parameter, NULL, 10);
      else if (strstr (line, "Number of collisional transitions:"))
        {
          for (int i = 0; i < RXI_MAX_COLL_PARTNERS; i++)
            {
              if (rxi->mc_par.cps[i].name != col_partner)
                continue;
              rxi->mc_par.cps[i].numof_coltr = strtof (parameter, NULL);
            }
        } 
      else if (strstr (line, "Collisional temperatures:"))
        {
          for (int i = 0; i < RXI_MAX_COLL_PARTNERS; i++)
            {
              if (rxi->mc_par.cps[i].name != col_partner)
                continue;
              char *end;
              int j = 0;
              for (float f = strtof (parameter, &end); parameter != end;
                                             f = strtof (parameter, &end))
                {
                  rxi->mc_par.cps[i].temps[j++] = f;
                  parameter = end;
                }
              rxi->mc_par.cps[i].numof_temps = j;
            }
        }
    }
  free (line);

  fclose (mi);
  return 0;
}

static int 
reading_enlev (struct radexi_data *rxi)
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
  unsigned int el = 0;
  while (fgets (line, n, enlev_file))
    {
      int t = 0;
      for (char *token = strtok (line, ","); token; token = strtok (NULL, ","))
        {
          if (t == 1)
            {
              rxi->mi.enlev[el].term = strtod (token, NULL);
            }
          else if (t == 2)
            {
              rxi->mi.enlev[el].statw = strtof (token, NULL);
            }
          else if (t > 2)
            {
              strcat (rxi->mi.enlev[el].qnums, token);
            }
          t++;
        }
      el++;
    }
  fclose (enlev_file);
  if (el != rxi->mi.numof_enlev)
    return -1;

  return 0;
}

static int 
reading_radtr (struct radexi_data *rxi)
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
  unsigned int rt = 0;
  while (fgets (line, n, radtr_file))
    {
      int t = 0;
      for (char *token = strtok (line, ","); token; token = strtok (NULL, ","))
        {
          /* May deprecate ulev & llev later, if they will have no other need */
          if (t == 1)
            {
              rxi->mi.radtr[rt].ulev = atoi (token);
            }
          else if (t == 2)
            {
              rxi->mi.radtr[rt].llev = atoi (token);
            }
          else if (t == 3)
            {
              rxi->mi.radtr[rt].a_einst = strtod (token, NULL);
            }
          else if (t == 4)
            {
              rxi->mi.radtr[rt].spfreq = strtof (token, NULL);
            }
          else if (t == 5)
            {
              rxi->mi.radtr[rt].enup = strtof (token, NULL);
            }
          else if (t > 5)
            {
              strcat (rxi->mi.radtr[rt].tail, token);
            }
          t++;
        }
      rxi->mi.radtr[rt].xnu = rxi->mi.enlev[rxi->mi.radtr[rt].ulev].term - 
                              rxi->mi.enlev[rxi->mi.radtr[rt].llev].term;
      rt++;
    }
  fclose (radtr_file);
  if (rt != rxi->mi.numof_radtr)
    return -1;

  return 0;
}

static int 
reading_colpart (struct radexi_data *rxi)
{
  for (unsigned int icp = 0; icp < rxi->mi.numof_colpart; icp++) 
    {
      char *cp_path = (char *) malloc (RXI_MOLECULE_MAX_SIZE + 20);
      char *cp_name = (char *) malloc (RXI_MOLECULE_MAX_SIZE);

      conv_int_to_name (rxi->mc_par.cps[icp].name, cp_name);
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
          int t = 0;
          float lcoef, ucoef = 0;
          float ltemp, utemp = 0;
          unsigned int lcollev, ucollev = 0;
          for ( char *token = strtok (line, ","); 
                token;
                token = strtok (NULL, ","))
            {
              if (t == 1)
                {
                  ucollev = atoi (token);
                }
              else if (t == 2)
                {
                  lcollev = atoi (token);
                }
              else if (t > 2)
                {
                  if (rxi->mc_par.cps[icp].temps[t-3] <= rxi->mc_par.Tkin)
                    {
                      lcoef = strtof (token, NULL);
                      ltemp = rxi->mc_par.cps[icp].temps[t-3];
                    }
                  else 
                    {
                      ucoef = strtof (token, NULL);
                      utemp = rxi->mc_par.cps[icp].temps[t-3];
                      break;
                    }
                }
              t++;
            }
          if (ltemp > utemp)
            rxi->mc_par.cps[icp].coef[ucollev-1][lcollev-1] = lcoef;
          else if (t == 3)
            rxi->mc_par.cps[icp].coef[ucollev-1][lcollev-1] = ucoef;
          else 
            rxi->mc_par.cps[icp].coef[ucollev-1][lcollev-1] = lcoef       + \
                                        (ucoef - lcoef) / (utemp - ltemp) * \
                                        (rxi->mc_par.Tkin - ltemp);
          icoef++;
        }
      fclose (cp_file);
      if (icoef != rxi->mc_par.cps[icp].numof_coltr)
        return -1;
    }
  return 0;
}

static int 
prepare_for_calculation (struct radexi_data *rxi)
{
  float matrix[RXI_MAX_ENLEV][RXI_MAX_ENLEV];
  for (unsigned int cp = 0; rxi->mc_par.cps[cp].name != 0; cp++)
    {
      rxi->mc_par.total_density += rxi->mc_par.cps[cp].dens;
      for (unsigned int i = 0; i < rxi->mi.numof_enlev; i++)
        for (unsigned int j = 0; j < rxi->mi.numof_enlev; j++)
          matrix[i][j] += rxi->mc_par.cps[cp].dens                        * \
                          rxi->mc_par.cps[cp].coef[i][j];
    }

  for (unsigned int i = 0; i < rxi->mi.numof_enlev; i++)
    {
      for (unsigned int j = 0; j < rxi->mi.numof_enlev; j++)
        {
          double ediff = rxi->mi.enlev[i].term - rxi->mi.enlev[j].term;
          if ((ediff > 0) && (hP*sol/kB * ediff / rxi->mc_par.Tkin < 160))
             rxi->mc_par.uprate[j][i] = rxi->mi.enlev[i].statw            / \
                                        rxi->mi.enlev[j].statw            * \
               exp (- sol * hP * ediff / kB / rxi->mc_par.Tkin)           * \
                                        matrix[i][j];
          else 
            rxi->mc_par.uprate[j][i] = 0;
          printf ("Uprate %d %d -> %e\n", j, i, rxi->mc_par.uprate[j][i]);
        }
    }
  for (unsigned int i = 0; i < rxi->mi.numof_enlev; i++)
  {
    for (unsigned int j = 0; j < rxi->mi.numof_enlev; j++)
    {
      rxi->mc_par.totrate[i] += rxi->mc_par.uprate[i][j];
      printf ("+= %e \n", rxi->mc_par.uprate[i][j]);
    }
    printf ("res %e \n", rxi->mc_par.totrate[i]);
  }
  return 0;
}

int
read_data (struct radexi_data *rxi)
{
  reading_info (rxi);
  reading_enlev (rxi);
  reading_radtr (rxi);
  reading_colpart (rxi);
  prepare_for_calculation (rxi);
  for (unsigned int i = 0; i < rxi->mi.numof_enlev; i++)
    printf ("%e\n", rxi->mc_par.totrate[i]);
  return 0;
}
