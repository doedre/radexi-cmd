/**
 * @file utils/csv.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>

#include "utils/csv.h"

#include "rxi_common.h"
#include "utils/debug.h"

RXI_STAT
rxi_csv_write_line (FILE *csv, const char *line)
{
  DEBUG ("Writing '%s' to .csv file", line);

  char *nline = malloc (RXI_STRING_MAX * sizeof (*nline));
  CHECK (nline && "Allocation failed");
  if (!nline)
    return RXI_ERR_ALLOC;

  strcpy (nline, line);
  bool is_first = true;
  for (char *token = strtok (nline, " "); token; token = strtok (NULL, " "))
    {
      if (!is_first)
        fprintf (csv, ",");

      // Remove control characters
      char *token_clean = token;
      for (int i = 0; token[i] != '\0'; ++i)
        {
          if (iscntrl ((uint8_t)token[i]))
            token_clean[i++] = ' ';
        }
      fprintf (csv, "%s", token_clean);
      is_first = false;
    }
  fprintf (csv, "\n");

  return RXI_OK;
}
