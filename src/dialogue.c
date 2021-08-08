/*	dialogue.c
 *
 * Contains realisations for functions, needed in dialogue w/ the user. Primary
 * goal is to collect all the needed info for calculations w/o troubles.
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

void
start_dialogue (struct rx_options *opts) 
{
  print_userguide ();
  printf ("Enter molecule's name ");
  printf ("('\x1B[4;37;40mlist\x1B[0;37;40m\x1B[3;37;40m'");
  printf (" to see the variants)\n");

  print_userinput ();

  int bytes_read;
  int size = 15;
  char *string;
  string = (char *) malloc (size);
  bytes_read = getdelim (&string, &size, ' ', stdin);

  if (bytes_read != -1)
    printf ("\r >> %s \t", string);

}

void
textformat (enum Attr attr, enum Color fg, enum Color bg)
{
  /* See https://www.linuxjournal.com/article/8603 for reference. */
  char command[13];
  sprintf (command, "%c[%d;%d;%dm", 0x1B, attr, fg + 30, bg + 40);
  printf ("%s", command);
}

void
print_userguide () 
{
  printf ("\x1B[0;37;40m\x1B[1;35;40m  ## \x1B[0;37;40m\x1B[3;37;40m");
}
  
void
print_userinput () 
{
  printf ("\x1B[0;37;40m\x1B[1;32;40m  >> \x1B[0;37;40m");
}
