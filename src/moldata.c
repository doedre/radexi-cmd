/*	moldata.c
 *
 * Used to create, rework and delete molecular data files from lamda database.
 * Original files are too big and required to be fully read before processing.
 * Now frequencies for the molecules will be sorted in the way to speed up the
 * calculation time.
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
#include "linenoise.h"

#include <sys/stat.h>
#include <ctype.h>
#include <string.h>

static void
write_molecular_energy_levels (FILE *molfile, const char *mol_name)
{
  char *line;
  size_t n = 50;
  line = (char *) malloc (n);

  char *molecular_enlev_file_name;
  molecular_enlev_file_name = (char *) malloc (RXI_MOLECULE_MAX_SIZE*2+16);
  sprintf  (molecular_enlev_file_name, 
            "data/%s/%s-enlev.csv", mol_name, mol_name);

  FILE *mev = fopen (molecular_enlev_file_name, "w");
  free (molecular_enlev_file_name);

//  fprintf (mev, "level,energy,weight,qn\n");
  while (fgets (line, n, molfile))
    {
      if (!strncmp (line, "!", 1))
        break;

      bool isFirst = true;
      for (char *token = strtok (line, " "); token; token = strtok (NULL, " "))
        {
          if (!isFirst)
              fprintf (mev, ",");
          /* Removes control characters  */
          char *token_clean = token;
          int i;
          for (i = 0; token[i] != '\0'; i++)
            {
              if (iscntrl ((unsigned char) token[i]))
                token_clean[i++] = ' ';
            }
          fprintf (mev, "%s", token_clean);
          isFirst = false;
        }
      fprintf (mev, "\n");
    }

  free (line);
  fclose (mev);
}

/* This function controls adding new molecules to the database. */
static void
add_molecular_file (char *mol_file_name, const struct rx_options *rx_opts)
{
  FILE *molfile = fopen (mol_file_name, "r");
  if (!molfile)
    {
      perror (mol_file_name);
      exit (EXIT_FAILURE);
    }

  char *folder_name;
  /* The name shouldn't contain more than 15 chars  */
  folder_name = (char *) malloc (20);
  strcpy (folder_name, "data/");
  strcat (folder_name, rx_opts->molecule_name);

  /* Create a folder if it doesn't exist  */
  struct stat sb;
  if (stat (folder_name, &sb) == -1)
      mkdir (folder_name, 0700);
  else 
    {
      printf ("\x1B[0;37;40m\x1B[1;35;40m  ## \x1B[0;37;40m");
      printf ("Specified molecule name already exists, rewrite it? (Y/n)");
      printf ("\x1B[0;37;40m\x1B[3;37;40m\n");

      char *choice;
      choice = (char *) malloc (1);
      while ((choice = linenoise ("  >> ")) != NULL)
        {
          if (!strcmp (choice, "y") || !strcmp (choice, "Y"))
            {
              printf ("TODO: Here it should delete all files\n");
              /* use stat + remove() or may stay the same */
              break;
            }
          else 
            {
              printf ("\x1B[0;37;40m\x1B[1;35;40m  ## \x1B[0;37;40m");
              printf ("Nothing to do, quitting...");
              printf ("\x1B[0;37;40m\n");
              exit (EXIT_SUCCESS);
            } 
        }
      free (choice);
    }

  char *molecular_info_file_name;
  molecular_info_file_name = (char *) malloc (RXI_MOLECULE_MAX_SIZE*2+10);
  sprintf (molecular_info_file_name, "data/%s/%s.info", 
                                      rx_opts->molecule_name, 
                                      rx_opts->molecule_name);

  FILE *mi = fopen (molecular_info_file_name, "w");
  free (molecular_info_file_name);

  char *line;
  size_t n = 50;
  line = (char *) malloc (n);

   for (int i = 0; i < 7; i++)
    {
      if (fgets (line, n, molfile))
        {
          if (i == 1)
            {
              fprintf (mi, "Molecule name: %s", line);
            } 
          else if (i == 3)
            {
              fprintf (mi, "Molecular weight: %s", line);
            }
          else if (i == 5)
            {
              fprintf (mi, "Number of energy levels: %s", line);
            }
          else 
            {
              continue;
            }
        }
    }  
   
  write_molecular_energy_levels (molfile, rx_opts->molecule_name);
  fgets (line, 15, molfile);
  fprintf (mi, "Number of radiative transitions: %s\n", line);

  fclose (molfile);
  fclose (mi);
  free (line);
  free (folder_name);
}

void 
operate_molecular_files (char *mol_file_name, const struct rx_options *rx_opts)
{
  if (rx_opts->usage_mode == UM_MOLECULAR_FILE_ADD)
    {
      add_molecular_file (mol_file_name, rx_opts);
    }
}
