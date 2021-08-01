/*	radexi.h
 *
 * Header for the whole radexi-cmd program. Contains function declarations, 
 * structures, enumerators and booleans for the flags.
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

#ifndef RADEXI_H
#define RADEXI_H

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#define PROGRAM_NAME  "radexi"


/* Giving live and colors to our program  */
enum Attr
{
  RESET = 0,
  BRIGHT,
  DIM,
  ITALIC,
  UNDERLINE,
};

enum Color
{
  BLACK = 0,
  RED,
  GREEN,
  YELLOW,
  BLUE,
  MAGENTA,
  CYAN,
  WHITE,
};


/* Use this function before and after the 'printf()' function to control the 
 * output flow  */
void textformat (enum Attr attr, enum Color fg, enum Color bg);

/* Defines how molecular cloud's parameters will be collected. */
enum UsageMode
{
  /* Whether via dialogue w/ the user (if no file done) */
  UM_DIALOGUE = 1,

  /* Or via special input file (if it was referenced)  */
  UM_FILE
};


/* This structure defines options for all of the possible processes in the 
 * program. It is being filled after reading options by 'getopt' in 
 * 'set_rx_options' functions */
struct rx_options
{
  /* How molecular cloud's parameters will be collecter: via dialogue or file.
   */
  enum UsageMode usage_mode;

  /* Path to the input file in case of UM_FILE usage_mode. */
  char inp_file_path[80];

  /* Path to the file w/ results. */
  char out_file_path[80];

  /* Whether resulting table should be printed in the terminal or not.
   * -o flag */
  bool cmd_output;

  /* If user defined a path for the file w/ resilts. 
   * -r: or --result flag */
  bool user_defined_out_file_path;

};


/* Exit statuses. */
enum
{
  RADEXI_NORMAL = 1,

  /* Option for major problems, like wrong options usage, incorrect path to the
   * files or even calculation problems.  */
  RADEXI_FAILURE
};


/* Checks what needs to be printed if there are some errors during the launch
 * (most commonly they appear because of the user's input). Whether prints a 
 * small message with future advice, or --help field output */
void cmd_usage (int status);

/* Defines how program will work (see comments for 'struct rx_options') and
 * returns the index of the first non-option argument (path to the input file).
 */
int set_rx_options (struct rx_options *rx_opts, int argc, char **argv);

/* Sets all parameters to defaults before proceeding w/ 'set_rx_options()'  */
void set_default_rx_options (struct rx_options *rx_opts);

/* Initiates a dialogue w/ the user if no other input had been set  */
void start_dialogue (struct rx_options *rx_opts);

#endif  // RADEXI_H
