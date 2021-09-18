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
write_csv  (FILE *molfile, const char *mol_name, const char *outfile_suf)
{
  char *line;
  size_t n = 200;
  line = (char *) malloc (n);

  char *csv_name;
  csv_name = (char *) malloc (RXI_MOLECULE_MAX_SIZE*2+16);
  sprintf  (csv_name, "data/%s/%s.csv", mol_name, outfile_suf);

  FILE *csv = fopen (csv_name, "w");
  free (csv_name);

  while (fgets (line, n, molfile))
    {
      if (!strncmp (line, "!", 1))
        break;

      bool isFirst = true;
      for (char *token = strtok (line, " "); token; token = strtok (NULL, " "))
        {
          if (!isFirst)
              fprintf (csv, ",");
          /* Removes control characters  */
          char *token_clean = token;
          int i;
          for (i = 0; token[i] != '\0'; i++)
            {
              if (iscntrl ((unsigned char) token[i]))
                token_clean[i++] = ' ';
            }
          fprintf (csv, "%s", token_clean);
          isFirst = false;
        }
      fprintf (csv, "\n");
    }

  free (line);
  fclose (csv);
}

/* Used to extract collision partner's name from the string without the user */
static enum ColPartner
extract_col_partner (char * line)
{
  enum ColPartner result = H2;
  /* Possible names */
  unsigned long numof_variants = 12;
  char *variants[] = { 
                        " p-h2 ", " para-h2 ", "para h2 ",
                        " o-h2 ", " ortho-h2 ", "ortho h2 ", 
                        " h2 ", "-h2 ", " h2-",
                        " he ", "-he ", "he-" 
                        }; 
  for (unsigned long i = 0; i < strlen (line); i++)
    {
      if (isalpha (line[i]))
        line[i] = tolower (line[i]);
    }

  unsigned long pos = 0;
  for (; pos < numof_variants; pos++)
    {
      if (strstr (line, variants[pos]))
        break;
    }

  if (pos < 3)
    result = PARA_H2;
  else if (pos < 6)
    result = ORTHO_H2;
  else if (pos < 9)
    result = H2;
  else if (pos < 12)
    result = He;
  else 
    result = NO_MOLECULE;

  return result;
}

static void
conv_int_to_name (int cp, char *cp_name)
{
  if (cp == H2)
    strcpy (cp_name, "H2");
  else if (cp == PARA_H2)
    strcpy (cp_name, "pH2");
  else if (cp == ORTHO_H2)
    strcpy (cp_name, "oH2");
  else if (cp == ELECTRONS)
    strcpy (cp_name, "electrons");
  else if (cp == HI)
    strcpy (cp_name, "HI");
  else if (cp == He)
    strcpy (cp_name, "He");
  else if (cp == HII)
    strcpy (cp_name, "HII");
  else 
    cp_name[0] = '\0';
}

/* Writes collision partner's parameters into the .info file.
 * returns:
 * - one of the collision partners as defined in ColPartner enum
 * - -1 if no collision partner was found
 * - 0 if it's an end of file  */
static int
define_collision_partner (FILE *molfile, FILE *mol_info_file)
{
  char *line;
  size_t n = 200;
  line = (char *) malloc (n);

  fgets (line, n, molfile);

  if (!strncmp (line, "!", 1))
    return 0;

  enum ColPartner cp = extract_col_partner (line);
  if (cp == NO_MOLECULE)
    return -1;

  fprintf (mol_info_file, "Collision partner: %d\n", cp);
  fgets (line, n, molfile);
  fgets (line, n, molfile);
  fprintf (mol_info_file, "Number of collisional transitions: %s", line);
  fgets (line, n, molfile);
  fgets (line, n, molfile);
  fgets (line, n, molfile);
  fgets (line, n, molfile);
  fprintf (mol_info_file, "Collisional temperatures: %s", line);
  fgets (line, n, molfile);

  free (line);
  return cp;
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

  /* File with information about the molecule? which cant be loaded in .csv 
   * database properly  */
  FILE *mi = fopen (molecular_info_file_name, "w");
  free (molecular_info_file_name);

  char *line;
  size_t n = 200;
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
   
  /* Writing .csv file with energy levels */
  write_csv (molfile, rx_opts->molecule_name, "enlev");

  fgets (line, n, molfile);
  fprintf (mi, "Number of radiative transitions: %s", line);
  fgets (line, n, molfile);

  /* Writing .csv file with radiative transitions */
  write_csv (molfile, rx_opts->molecule_name, "radtr");

  fgets (line, n, molfile);
  fprintf (mi, "Number of collision partners: %s\n", line);
  fgets (line, n, molfile);

  /* Writing the names of collision partners (user defines them by himself as 
   * there is no standart in LAMDA database)  */
  int cp = 0;
  while ((cp = define_collision_partner(molfile, mi)) > 0)
    {
      char *cp_name;
      cp_name = (char *) malloc (10);
      conv_int_to_name (cp, cp_name);
      write_csv (molfile, rx_opts->molecule_name, cp_name);
      free (cp_name);
    }
  if (cp == -1)
    printf ("Error with collision partners\n");
  else 
    printf ("All is ok\n");

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
