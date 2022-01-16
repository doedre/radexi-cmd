/**
 * @file rxi_common.h
 * @brief Provides common data structures for calculations.
 */

#ifndef RXI_COMMON_H
#define RXI_COMMON_H

#include <stdint.h>
#include <stdbool.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_const_cgsm.h>

#define RXI_PATH_MAX 1024
#define RXI_STRING_MAX 512
#define RXI_QNUM_MAX 30

#define RXI_MOLECULE_MAX 15
#define RXI_COLL_TEMPS_MAX 50
#define RXI_COLL_PARTNERS_MAX 7

#define RXI_ELEMENTS_MAX RXI_COLL_TEMPS_MAX + 3

#define RXI_FK                                                                \
        GSL_CONST_CGSM_PLANCKS_CONSTANT_H * GSL_CONST_CGSM_SPEED_OF_LIGHT     \
        / GSL_CONST_CGSM_BOLTZMANN

#define RXI_HP GSL_CONST_CGSM_PLANCKS_CONSTANT_H
#define RXI_SOL GSL_CONST_CGSM_SPEED_OF_LIGHT

/// @brief Status codes for functions that may fail.
///
/// Status codes should be returned by all of the functions, that may fail
/// during memory allocation, filesystem errors or just local mistakes in
/// molecular databases. They are divided by two types: warnings and errors.
/// Errors may lead to unexpected behavior. Warnings are expected not to
/// terminate the program with segmentation faults.
typedef enum RXI_STAT
{
  RXI_OK = 0,             //!< Everything is ok.
  RXI_ERR_ALLOC,          //!< Error on memory allocation.
  RXI_ERR_OPTS,           //!< Error in command line options.
  RXI_ERR_FILE,           //!< File opening error.
  RXI_ERR_CONV,           //!< Error with type conversion.
  RXI_WARN_LIMITS = 10,   //!< 
  RXI_WARN_LAMDA,         //!< LAMDA's information mismatch.
  RXI_WARN_NOFILE,
  RXI_FILE_END
}
RXI_STAT;

/// @brief Used to specify how the program will be used.
///
/// Only used in @ref `struct rxi_options` to set the program in specified
/// state.
enum USAGE_MODE
{
  UM_NONE = 0,                //!< Undefined usage mode.
  UM_DIALOGUE,                //!< Standart dialog with the user.
  UM_FILE,                    //!< Read input parameters from the file.
  UM_MOLECULAR_FILE_ADD,      //!< Add molecular data file from LAMDA.
  UM_MOLECULAR_FILE_DELETE,   //!< Delete local molecular data file.
  UM_MOLECULAR_FILE_LIST,     //!< List local molecular data files.
  UM_HELP,                    //!< Print help information.
  UM_VERSION                  //!< Print version information.
};


/* This structure defines options for all of the possible processes in the 
 * program. It is being filled after reading options by 'getopt_long()' in 
 * 'set_rx_options' functions. */
struct rxi_options
{
  /* How molecular cloud's parameters will be collected: via dialogue or file.
   */
  enum USAGE_MODE usage_mode;

  enum RXI_STAT status;

  /* Molecule name defined by the user when adding one in case of 
   * UM_MOLECULAR_FILE usage_mode.  */
  char molecule_name[RXI_MOLECULE_MAX];

  /* Forces file rewrite if specified files for results of calculation 
   * already exist.
   * -f flag: disabled by default */
  bool force_fs;

  bool no_freq_limits;

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

  /* Specified path for user's output file  */
  char result_path[RXI_PATH_MAX];
};

/// @brief Get `$(HOME)/.local/share/radexi/` path.
///
/// Allocates memory for the returned string, so it should be freed after usage
/// to avoid memory leak.
/// @return On success returns `$(HOME)/.local/share/radexi/` path string. On
/// error returns `NULL`.
const char *rxi_database_path ();

/// @brief Get `$(HOME)/.config/radexi/` path.
///
/// Allocates memory for the returned string, so it should be freed after usage
/// to avoid memory leak.
/// @return On success returns `$(HOME)/.config/radexi/` path string. On error
/// returns `NULL`.
const char *rxi_config_path ();

/// @brief Names of the possible collision partners from LAMDA.
typedef enum COLL_PART
{
  H2 = 1,
  PARA_H2,
  ORTHO_H2,
  ELECTRONS,
  HI,
  He,
  HII,
  NO_PARTNER = 0
}
COLL_PART;

typedef enum GEOMETRY
{
  SPHERE = 1,
  SLAB,
  LVG,
  OTHER = 0
}
GEOMETRY;

/// @brief TODO
struct rxi_input_data
{
  char name[RXI_MOLECULE_MAX];
  float sfreq;
  float efreq;
  double temp_kin;
  double temp_bg;
  double col_dens;
  double line_width;
  GEOMETRY geom;
  int8_t n_coll_partners;

  COLL_PART coll_part[RXI_COLL_PARTNERS_MAX];
  double coll_part_dens[RXI_COLL_PARTNERS_MAX];
};

/// @brief Used to store molecular information from `*.info` file or LAMDA.
struct rxi_db_molecule_info
{
  char    *name;
  float   weight;
  int     numof_enlev;
  int     numof_radtr;
  int     numof_coll_part;

  COLL_PART *coll_part;
  int     *numof_coll_trans;
  int8_t  *numof_coll_temps; 
  gsl_matrix *coll_temps;
};

/// @brief TODO
RXI_STAT rxi_db_molecule_info_malloc (struct rxi_db_molecule_info **mol_info);

/// @brief TODO
void rxi_db_molecule_info_free (struct rxi_db_molecule_info *mol_info);

/// @brief Holds energy level information from database.
struct rxi_db_molecule_enlev
{
  int     *level;
  double  *term;
  double  *weight;
  char    **qnum;
};

/// @brief Memory allocation for `struct rxi_db_molecule_enlev`.
RXI_STAT rxi_db_molecule_enlev_malloc (struct rxi_db_molecule_enlev **mol_enl,
                                       const size_t n_enlev);

/// @brief TODO
void rxi_db_molecule_enlev_free (struct rxi_db_molecule_enlev *mol_enl);

/// @brief Holds radiative transfer information from database.
struct rxi_db_molecule_radtr
{
  int     *up;
  int     *low;
  double  *einst;
  double  *freq;
  double  *up_en;
};

/// @brief TODO
RXI_STAT rxi_db_molecule_radtr_malloc (struct rxi_db_molecule_radtr **mol_rat,
                                       const size_t n_radtr);

/// @brief TODO
void rxi_db_molecule_radtr_free (struct rxi_db_molecule_radtr *mol_rat);

/// @brief Holds collisional coefficients for specific molecule.
struct rxi_db_molecule_coll_part
{
  int   *up;
  int   *low;
  gsl_matrix *coll_rates;
};

/// @brief TODO
RXI_STAT rxi_db_molecule_coll_part_malloc (
    struct rxi_db_molecule_coll_part **mol_cp, const size_t n_cp_trans,
    const size_t n_temps);

/// @brief TODO
void rxi_db_molecule_coll_part_free (struct rxi_db_molecule_coll_part *mol_cp);


struct rxi_db_molecule
{
  struct rxi_db_molecule_info   *info;
  struct rxi_db_molecule_enlev  *enlev;
  struct rxi_db_molecule_radtr  *radtr;
  struct rxi_db_molecule_coll_part **coll_part;
};

/// @brief TODO
struct rxi_calc_data
{
  double *temp_kin;
  double *temp_bg;
  double *col_dens;
  double *line_width;
  GEOMETRY *geom;
  int *up;
  int *low;

  gsl_vector *term;
  gsl_vector *weight;
  gsl_matrix *einst;
  gsl_matrix *energy;
  gsl_matrix *coll_rates;
  gsl_vector *tot_rates;
  gsl_matrix *bgfield;

  gsl_matrix *rates;
  gsl_vector *pop;
  gsl_matrix *tau;
  gsl_matrix *exit_temp;
};

/// @brief TODO
RXI_STAT rxi_calc_data_malloc (struct rxi_calc_data **calc_data,
                               const size_t n_enlev);

/// @brief TODO
void rxi_calc_data_free (struct rxi_calc_data *calc_data);

#endif  // RXI_COMMON_H
