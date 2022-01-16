/**
 * @file core/calculation.h
 * @brief Defines high level functions for main calculations.
 */

#include "rxi_common.h"

/// @brief Fills `struct rxi_calc_data` to start calculations.
RXI_STAT rxi_calc_data_init (struct rxi_calc_data *calc_data,
                             const struct rxi_input_data *inp_data,
                             const struct rxi_db_molecule_info *mol_info);

/// @brief TODO
RXI_STAT rxi_calc_data_fill (const struct rxi_input_data *inp_data,
                             const struct rxi_db_molecule_info *mol_info,
                             const struct rxi_db_molecule_enlev *mol_enlev,
                             const struct rxi_db_molecule_radtr *mol_radtr, 
                             struct rxi_db_molecule_coll_part **mol_cp,
                             struct rxi_calc_data *calc_data);

RXI_STAT rxi_calc_find_rates (struct rxi_calc_data *data, const int n_enlev,
                              const int n_radtr);

double rxi_calc_crate (const double istat, const double jstat,
    const double ediff, const double kin_temp, const double crate);

double rxi_calc_escape_prob (const double tau, const GEOMETRY geom);

double rxi_calc_optical_depth (const double coldens, const double line_width,
    const double energy, const double einst, const double ustat,
    const double lstat, const double upop, const double lpop);
