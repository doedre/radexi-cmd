/**
 * @file core/calculation.h
 * @brief Defines high level functions for main calculations.
 */

#include "rxi_common.h"

/// @brief Fills `struct rxi_calc_data` to start calculations.
RXI_STAT rxi_calc_data_init (const struct rxi_input_data *inp_data);

/// @brief TODO
RXI_STAT rxi_calc_data_fill (const struct rxi_input_data *inp_data,
                             const struct rxi_db_molecule_info *mol_info,
                             const struct rxi_db_molecule_enlev *mol_enlev,
                             const struct rxi_db_molecule_radtr *mol_radtr, 
                             struct rxi_db_molecule_coll_part **mol_cp,
                             struct rxi_calc_data *calc_data);
