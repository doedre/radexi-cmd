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

//! Program version.
#define RXI_VERSION "0.2"

//! Maximum path size.
#define RXI_PATH_MAX 1024
//! Maximum string size to read from files.
#define RXI_STRING_MAX 512
//! Maximum string size for quantum numbers.
#define RXI_QNUM_MAX 30
//! Maximum molecule name size.
#define RXI_MOLECULE_MAX 50
//! Maximum number of collisional temperatures.
#define RXI_COLL_TEMPS_MAX 50
//! Maximum number of collisional parameters.
#define RXI_COLL_PARTNERS_MAX 7
//!
#define RXI_ELEMENTS_MAX 53

//!
#define RXI_FK                                                                \
        GSL_CONST_CGSM_PLANCKS_CONSTANT_H * GSL_CONST_CGSM_SPEED_OF_LIGHT     \
        / GSL_CONST_CGSM_BOLTZMANN

//! Planck constant
#define RXI_HP GSL_CONST_CGSM_PLANCKS_CONSTANT_H
//! Speed of light
#define RXI_SOL GSL_CONST_CGSM_SPEED_OF_LIGHT
//! Boltzmann constant
#define RXI_KB GSL_CONST_CGSM_BOLTZMANN

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
  UM_FIND_GOOD_FIT,           //!< Finding the best fit for entered intensities.
  UM_MOLECULAR_FILE_ADD,      //!< Add molecular data file from LAMDA.
  UM_MOLECULAR_FILE_DELETE,   //!< Delete local molecular data file.
  UM_MOLECULAR_FILE_LIST,     //!< List local molecular data files.
  UM_HELP,                    //!< Print help information.
  UM_VERSION                  //!< Print version information.
};

/// @brief Options to set program's global state.
///
/// This program acts like state machine and these options (defined through
/// flags on start) are used to set the global state. In my opinion this
/// behaviour is not good from the architecture perspective, but it may help
/// someone include this program as package in their own programs. Look for
/// `main.c` file for examples.
struct rxi_options
{
  //! Used to set program's main state. Look for `enum USAGE_MODE`.
  enum USAGE_MODE usage_mode;

  //! Write status during `rxi_set_options()` to exit on fail.
  enum RXI_STAT status;

  //! Molecule name for molecular file usage mods is written here.
  char molecule_name[RXI_MOLECULE_MAX];

  //! Forces every action without user interruption. `-f` option.
  bool force_fs;

  //! Disable frequency limits for calculations. `-l` option.
  bool no_freq_limits;

  //! Print result in `stdout`. `-o` option.
  bool cmd_output;

  //! Do not create file for results. `-x` option.
  bool no_result_file;

  //! Don't print starting message. `-q` option.
  bool quite_start;

  //! Enter powers of 10 instead of whole density. `-L` option.
  bool dens_log_scale;

  bool hz_width;

  //! Path to the file with results. `-r` or `--result` option.
  bool user_defined_out_file_path;

  //! Path for result file (`-r` option) is written here.
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

/// @brief Possible geometries for radiation fields.
typedef enum GEOMETRY
{
  SPHERE = 1,
  SLAB,
  LVG,
  OTHER = 0
}
GEOMETRY;

/// @brief Starting information may be written here.
///
/// All the needed information from user can be written in this structure for
/// future use. No memory should be allocated before.
struct rxi_input_data
{
  char    name[RXI_MOLECULE_MAX]; //!< Molecule name from local database.
  char    names[RXI_MOLECULE_MAX];//!< Molecule name from local database.
  char    name_list[10][15];      //!< Molecule name from local database.
  int8_t  numof_molecules;        //!< Number of molecules.
  float   sfreq;                  //!< Starting frequency for output [GHz].
  float   efreq;                  //!< Ending frequency for output [GHz].
  double  temp_kin;               //!< Kinetic temperature [K].
  double  temp_bg;                //!< Background temperature [K].
  double  col_dens;               //!< Column density [cm-2].
  double  line_width;             //!< FWHM width for all lines [km s-1].
  GEOMETRY geom;                  //!< Radiation field geometry.
  int8_t  n_coll_partners;        //!< Number of specified collision partners.

  COLL_PART coll_part[RXI_COLL_PARTNERS_MAX];   //!< Collision partner names.
  double coll_part_dens[RXI_COLL_PARTNERS_MAX]; //!< Partner densities [cm-3].
};

/// @brief Used to read molecular information from `*.info` file or LAMDA.
///
/// This structure shouldn't be filled by the user. It is used to store
/// information from database for future allocations. Should allocate memory
/// by `rxi_db_molecule_info_malloc()` before usage (TODO: allocation is not
/// needed).
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

/// @brief Allocate memory for `struct rxi_db_molecule_info`.
/// @param **mol_info -- pointer to a pointer to a structure for allocation.
/// @return `RXI_OK` on success; `RXI_ERR_ALLOC` on allocation error.
RXI_STAT rxi_db_molecule_info_malloc (struct rxi_db_molecule_info **mol_info);

/// @brief Free memory for `struct rxi_db_molecule_info`.
/// @param *mol_info -- pointer to a structure which needs to be freed.
void rxi_db_molecule_info_free (struct rxi_db_molecule_info *mol_info);

/// @brief Holds energy level information from database.
///
/// This structure shouldnt be filled by the user. It is used to store
/// information about energy levels from database. Should allocate memory
/// by `rxi_db_molecule_enlev_malloc()` before usage.
struct rxi_db_molecule_enlev
{
  int     *level;
  double  *term;
  double  *weight;
};

/// @brief Memory allocation for `struct rxi_db_molecule_enlev`.
/// @param **mol_enl -- pointer to a pointer to a structure for allocation;
/// @param n_enlev -- number of energy levels for current molecule (get it
/// from database by `rxi_db_read_molecule_info()` function).
/// @return `RXI_OK` on success; `RXI_ERR_ALLOC` on allocation error.
RXI_STAT rxi_db_molecule_enlev_malloc (struct rxi_db_molecule_enlev **mol_enl,
                                       const size_t n_enlev);

/// @brief Free memory for `struct rxi_db_molecule_enlev`.
/// @param *mol_enl -- pointer to a structure which needs to be freed.
void rxi_db_molecule_enlev_free (struct rxi_db_molecule_enlev *mol_enl);

/// @brief Holds radiative transfer information from database.
///
/// This structure shouldn't be filled by the user. It is used to store
/// information about radiative transfers from database. Should allocate memory
/// by `rxi_db_molecule_radtr_malloc()` before usage.
struct rxi_db_molecule_radtr
{
  int     *up;
  int     *low;
  double  *einst;
  double  *freq;
  double  *up_en;

  // These parameters are entered during finding the best fit.
  double  *intensity;
  double  *sigma;
  double  *fwhm;
};

/// @brief Memory allocation for `struct rxi_db_molecule_radtr`.
/// @param **mol_rat -- pointer to a pointer to a structure for allocation;
/// @param n_radtr -- number of radiative transfers for current molecule (get
/// it from database by `rxi_db_read_molecule_info()` function).
/// @return `RXI_OK` on success; `RXI_ERR_ALLOC` on allocation error.
RXI_STAT rxi_db_molecule_radtr_malloc (struct rxi_db_molecule_radtr **mol_rat,
                                       const size_t n_radtr);

/// @brief Free memory for `struct rxi_db_molecule_radtr`.
/// @param *mol_rat -- pointer to a structure which needs to be freed.
void rxi_db_molecule_radtr_free (struct rxi_db_molecule_radtr *mol_rat);

/// @brief Holds collisional coefficients for specific molecule.
///
/// This structure shouldn't be filled by the user. It is used to store
/// information about collisional partner from database. Should allocate memory
/// by `rxi_db_molecule_coll_part_malloc()` before usage.
struct rxi_db_molecule_coll_part
{
  int   *up;
  int   *low;
  gsl_matrix *coll_rates;
};

/// @brief Memory allocation for `struct rxi_db_molecule_coll_part`.
/// @param **mol_cp -- pointer to a pointer to a structure for allocation;
/// @param n_cp_trans -- number of collisional transitions for current
/// molecular partner (get it from database by `rxi_db_read_molecule_info()`
/// function).
/// @param n_temps -- number of colisional temperatures for current
/// molecular partner (get it from database by `rxi_db_read_molecule_info()`
/// function).
/// @return `RXI_OK` on success; `RXI_ERR_ALLOC` on allocation error.
RXI_STAT rxi_db_molecule_coll_part_malloc (
    struct rxi_db_molecule_coll_part **mol_cp, const size_t n_cp_trans,
    const size_t n_temps);

/// @brief Free memory for `struct rxi_db_molecule_coll_part`.
/// @param *mol_cp -- pointer to a structure which needs to be freed.
void rxi_db_molecule_coll_part_free (struct rxi_db_molecule_coll_part *mol_cp);

/// @brief Holds all information about the molecule from database.
///
/// This structure shouldn't be filled by the user. All components should be
/// allocated before including in this structure.
struct rxi_db_molecule
{
  struct rxi_db_molecule_info   *info;
  struct rxi_db_molecule_enlev  *enlev;
  struct rxi_db_molecule_radtr  *radtr;
  struct rxi_db_molecule_coll_part **coll_part;
};

/// @brief Holds all information for calculation and output.
struct rxi_calc_data
{
  struct rxi_input_data input;
  size_t numof_enlev;
  size_t numof_radtr;
  double chisq;
  int *up;
  int *low;

  gsl_vector *term;
  gsl_vector *weight;
  gsl_matrix *einst;
  gsl_matrix *freq;
  gsl_matrix *coll_rates;
  gsl_vector *tot_rates;
  gsl_matrix *bgfield;

  gsl_matrix *rates_archive;
  gsl_matrix *rates;
  gsl_vector *pop;
  gsl_matrix *tau;
  gsl_matrix *excit_temp;
  gsl_matrix *antenna_temp;
  gsl_matrix *radiation_temp;
};

/// @brief Memory allocation for `struct rxi_calc_data`.
/// @param **calc_data -- pointer to a pointer to a structure for allocation;
/// @param n_enlev -- number of energy levels for current molecule (get it from
/// database by `rxi_db_read_molecule_info()` function).
/// @param n_radtr -- number of radiative transfers for current molecule (get
/// it from database by `rxi_db_read_molecule_info()` function).
/// @return `RXI_OK` on success; `RXI_ERR_ALLOC` on allocation error.
RXI_STAT rxi_calc_data_malloc (struct rxi_calc_data **calc_data,
                               const size_t n_enlev, const size_t n_radtr);

/// @brief Free memory for `struct rxi_db_molecule_coll_part`.
/// @param *mol_cp -- pointer to a structure which needs to be freed.
void rxi_calc_data_free (struct rxi_calc_data *calc_data);

/// @brief For output results sorting.
struct rxi_calc_results
{
  int up;
  int low;
  char name[RXI_MOLECULE_MAX];
  double xnu;
  double spfreq;
  double tau;
  double population;
  double excit_temp;
  double antenna_temp;
  double upop;
  double lpop;
};

/// @brief Converts `enum GEOMETRY` to string.
char *geomtoname (GEOMETRY geom);

/// @brief Converts collisional partner number to string.
///
/// Allocates memory for the returned string, so it should be freed after
/// usage. Convert `COLL_PART` to predefined string names.
/// @param cp -- collision partner of `enum COLL_PART` type.
/// @return Allocated string with collision partner's name.
char *numtoname (COLL_PART cp);

/// @brief Converts string to collisional partner number.
///
/// Converts collision partner's name to `enum COLL_PART`.
/// @param *name -- string with name.
/// @return One of the collisional partners or `NO_MOLECULE` on error.
COLL_PART nametonum (const char *name);

/// @brief Get number for specified collision partner from database info file.
///
/// Searches database (through `struct rxi_db_molecule_info`) for specified
/// collision partner number in database (not `enum COLL_PART`).
/// @param *mol_info -- structure with information about the molecule from
/// local database;
/// @param cp -- collision partner which number needs to be defined.
/// @return Count number for specified collision partner.
int8_t cptonum (const struct rxi_db_molecule_info *mol_info, COLL_PART cp);

/// @brief Remove spaces from passed string.
void remove_spaces (char *str);

#endif  // RXI_COMMON_H
