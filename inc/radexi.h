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

#include <gsl/gsl_matrix.h>

#define PROGRAM_NAME          "radexi"
#define RXI_MOLECULE_MAX_SIZE 15  /* Max lenght for the molecule name 
                                     defined by the user                    */
#define RXI_MAX_EXCESS_SIZE   10  /* Max characters to read from .csv line  */
#define RXI_MAX_ENLEV         300 /* Max number of energy levels            */
#define RXI_MAX_RADTR         1000  /* Max number of radiative transitions  */
#define RXI_MAX_COLL_PARTNERS 7   /* Max number of collision partners       */
#define RXI_MAX_COLL_TEMPS    30  /* Max number of collisional temperatures */
#define RXI_MAX_COLL_TRANS    3000  /* Max number of collision transitions  */


/* Physical constants */
extern const double sol;              /* speed of light       [cm s-1]      */
extern const double hP;               /* Planck's constant    [erg s]       */
extern const double kB;               /* Boltzman's constant  [erg K-1]     */

extern const double fk;               /* hP * sol / kB                      */

/* Defines how molecular cloud's parameters will be collected. */
enum UsageMode
{
  /* Whether via dialogue w/ the user (if no file given). */
  UM_DIALOGUE = 1,

  /* Or via special input file (if it was referenced).  */
  UM_FILE,

  /* Otherwise it can be used to manipulate molecular data files. */
  UM_MOLECULAR_FILE_ADD,
  UM_MOLECULAR_FILE_DELETE,
  UM_MOLECULAR_FILE_LIST
};


/* This structure defines options for all of the possible processes in the 
 * program. It is being filled after reading options by 'getopt_long()' in 
 * 'set_rx_options' functions. */
struct rxi_options
{
  /* How molecular cloud's parameters will be collected: via dialogue or file.
   */
  enum UsageMode usage_mode;

  /* Path to the input file in case of UM_FILE usage_mode. */
  char inp_file_path[80];

  /* Path to the file w/ results. */
  char out_file_path[80];

  /* Molecule name defined by the user when adding one in case of 
   * UM_MOLECULAR_FILE usage_mode.  */
  char molecule_name[RXI_MOLECULE_MAX_SIZE];

  /* Whether resulting table should be printed in the terminal or not.
   * -o flag: disabled by default */
  bool cmd_output;

  /* Print message w/ information on the start or not.
   * -q flag: disabled by default */
  bool quite_start;

  /* If 'true', then instead of writing densities in scientific format user 
   * should only write power of 10.
   * -L of --log-density flag: disabled by default */
  bool dens_log_scale;

  /* if 'true', then use [GHz] line width instead of [km -s].
   * -H og --hz-width flag: disabled by default  */
  bool hz_width;

  /* If user defined a path for the file w/ resilts. 
   * -r or --result flag */
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


/* Names of the possible collision partners for RADEX.  */
enum ColPart
{
  H2 = 1,
  PARA_H2,
  ORTHO_H2,
  ELECTRONS,
  HI,
  He,
  HII,
  NO_MOLECULE
};

/* Used to define cloud's geometry, as RADEX is doing it's calculations in 
 * prediction of some uniformity. */
enum Geometry
{
  SPHERE = 1,
  SLAB,
  LVG
};

/* All needed information about the collision partner. May define several of 
 * them to cover all possible collision partners for the molecule.  */
struct colpart
{
  enum ColPart  name;                         /* collision partner's name 
                                                defined as constant integer */
  double        dens;                         /* collision partner's number
                                                density [cm-3]              */
  unsigned int  numof_coltr;                  /* num of collisional 
                                                transitions                 */
  unsigned int  numof_temps;                  /* num of collisional 
                                                temperatures                */
  unsigned int  *ucollev;                     /* upper collisional level    */
  unsigned int  *lcollev;                     /* lower collisional level    */
  float         temps[RXI_MAX_COLL_TEMPS];    /* collisional temperatures [K] */

  gsl_matrix    *K;                           /* collisional rate coefficients
                                                 (from database) [cm3 s-1].
                                                 It's dimensions are 
                                                 'numof_radtr X numof_radtr', 
                                                 but need to allocate memory
                                                 first.                     */
};

struct rxi_input
{
  char          name[RXI_MOLECULE_MAX_SIZE];
  double        Tkin;
  double        Tbg;
  double        coldens;
  double        fwhm;
  enum Geometry g;
  unsigned int  numof_colpart;
  unsigned int  max_numof_colpart;
  double        weight;
  unsigned int  numof_enlev;
  unsigned int  numof_radtr;
  struct colpart cps[RXI_MAX_COLL_PARTNERS];
};

/* Molecular cloud parameters used as starting data for calculations.  */
struct mc_par
{
  double        Tkin;                         /* kinetic temperature        */
  double        Tbg;                          /* background temparature     */
  double        coldens;                      /* molecular column density   */
  double        line_width;                   /* line width for all lines   */
  double        total_density;                /* total partner's density    */
  enum Geometry geom;                         /* Cloud's geometry           */
  unsigned int  numof_colpart;                /* num of collision partners  */
  struct colpart *cp;                         /* collision partners         */
};

struct bg_field
{
  double        intens[RXI_MAX_RADTR];        /* background intensity 
                                                 [erg s-1 cm-2 Hz-1 sr-1]   */
};

/* Data from enlev.csv will be stored here  */
struct enlev
{
  double        term;                         /* energy level [cm-1]        */
  double        statw;                        /* statistical weight         */
  char          qnums[RXI_MAX_EXCESS_SIZE];   /* quantum numbers in one string
                                              (writing excess strings here) */
};

/* Data from radtr.csv will be stored here  */
struct radtr
{
  unsigned int  ulev;                         /* upper level number         */
  unsigned int  llev;                         /* lower level number         */
  double        a_einst;                      /* Einstein's coefficient [s-1]  */
  float         spfreq;                       /* spectral frequency [GHz]   */
  double        enup;                         /* energy of upper level [K]  */
  double        xnu;                          /* line frequency in [cm-1]   */
  char          tail[RXI_MAX_EXCESS_SIZE];    /* excess info                */
};

/* Base information about the molecule. Used for data control: whether 
 * specified databases contain stated amount of parameters. */
struct mol_info
{
  char          name[RXI_MOLECULE_MAX_SIZE];  /* molecule name              */
  double        weight;                       /* molecular weight [a.m.u.]  */
  unsigned int  numof_enlev;                  /* num of energy levels       */
  unsigned int  numof_radtr;                  /* num of radiative transitions */ 
  struct enlev  *el;                          /* energy levels info         */
  struct radtr  *rt;                           /* radiative transitions info */

  gsl_matrix    *C;                           /* collision rates per second
                                                 for all collision partners 
                                                 [s-1]. Defined as 
                                                      C_ij += K_ij*ndens
                                                 To allocate memory needs 
                                                 numof_enlev dimensions.    */
  gsl_vector    *Ctot;                        /* total collision rate [s-1].
                                                 Defined as
                                                      Ctot_i = C_ij
                                                To allocate memory needs 
                                                numof_enlev dimensions.     */
};

/* Main structure which defines the input information for future calculations.
 * Be cautious with gsl_matrix and gsl_vector elements, as they need to be 
 * allocated before usage.                                                  */
struct rxi_data
{
  struct mc_par   mc;       /* Storing molecular clouds parameters here, which
                               were mostly entered by the user. Additional info
                               is written in collision partner's section, where
                               collisional coefficients and etc are defined.*/
  struct mol_info mi;       /* Any information in here is written from database
                               files. Pre-calculation is also done in here  */
  struct bg_field bg;       /* Structure for background field options.      */
};

/* Shorter/readable paths to frequently used variables.  */
# define n_el         mi.numof_enlev
# define n_rt         mi.numof_radtr
# define energy_level mi.el
# define rad_transfer mi.rt
# define coll_partner mc.cp

/* Results of the calculations will be stored here. Trying to avoid carrying
 * too much unnecessary information in here.                                */
struct radexi_results
{
};

struct rxi_data *rxi_data_calloc (const struct rxi_input *inp);
void rxi_data_free (struct rxi_data *rxi);

/* Checks what needs to be printed if there are some errors during the launch
 * (most commonly they appear because of the user's input). Whether prints a 
 * small message with future advice, or --help field output. */
void cmd_usage (int status);

/* Defines how program will work (see comments for 'struct rxi_options') and
 * returns the index of the first non-option argument (path to the input file).
 */
int set_rxi_options (struct rxi_options *rx_opts, int argc, char **argv);

/* Sets all parameters to defaults before proceeding w/ 'set_rxi_options()'  */
void set_default_rxi_options (struct rxi_options *rx_opts);

/* Initiates a dialogue w/ the user if no other input had been set. Fills 
 * sfreq with starting frequency, efreq with ending frequency, and also 
 * fills rxi.mc part with the basic information about the molecular cloud. */
void start_dialogue (float *sfreq, 
                     float *efreq,   
                     struct rxi_input *inp,
                     const struct rxi_options *rx_opts);

/* Reads .info file.
 * Returns:
 * - 0 on success */
int read_data (struct rxi_data *rxi);

int read_info_file (struct rxi_input *inp);

/* Used to convert enum with collision partners to their names.  */
void conv_int_to_name (int cp, char *cp_name);

int calculate_bg_field (struct rxi_data *rxi);

void main_calculations (struct rxi_data *rxi, struct radexi_results *rxi_res);

/* Reads usage_mode and decides what to do with the given molecular file  */
void operate_molecular_files (char *mol_file, const struct rxi_options *rx_opts);

#endif  // RADEXI_H
