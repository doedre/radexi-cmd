/**
 * @file utils/csv.h
 * @brief Simple csv manager for program needs.
 */

#include <stdio.h>

#include "rxi_common.h"

/// @brief Parse line and write it to csv (coma separated).
RXI_STAT rxi_csv_write_line (FILE *csv, const char *line);

/// @brief Read parsed line from the file to buffer.
RXI_STAT rxi_csv_read_line (FILE *csv, void **buff);
