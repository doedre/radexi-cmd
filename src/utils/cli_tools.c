/**
 * @file utils/cli_tools.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <ctype.h>

#include "linenoise/linenoise.h"

#include "utils/cli_tools.h"
#include "utils/debug.h"

char *rxi_readline (const char *prompt)
{
  DEBUG ("Call readline");

  return linenoise (prompt);
}

bool rxi_readline_accept ()
{
  DEBUG ("Call readline");

  bool result = false;
  char *line;
  if ((line = linenoise ("  [Y/n]>> ")) != NULL)
    {
      uint8_t symbol = tolower (line[0]);
      DEBUG ("Got %s -> %c", line, symbol);
      if (symbol == 'y')
        result = true;

      free (line);
    }

  return result;
}

RXI_STAT rxi_history_save (const char *line, const char *filename)
{
  DEBUG ("Save command '%s' to `%s'", line, filename);

  int ln_status = linenoiseHistoryAdd (line);
  if (ln_status != 1)
    {
      DEBUG ("Command '%s' cannot be added to history", line);
      return RXI_OK;
    }

  ln_status = linenoiseHistorySave (filename); 
  CHECK ((ln_status == 0) && "linenoise error");
  if (ln_status != 0)
    return RXI_ERR_FILE;

  return RXI_OK;
}

RXI_STAT rxi_history_load (const char *filename)
{
  DEBUG ("Load command history from `%s'", filename);

  linenoiseHistoryReset ();
  int ln_status = linenoiseHistoryLoad (filename);
  CHECK ((ln_status == 0) && "linenoise error");
  if (ln_status != 0)
    return RXI_ERR_FILE;

  return RXI_OK;
}
