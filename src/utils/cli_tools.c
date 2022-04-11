/**
 * @file utils/cli_tools.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
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
      DEBUG ("Got %s", line);
      uint8_t symbol = tolower (line[0]);
      DEBUG ("Got %s -> %c", line, symbol);
      if (symbol == 'y')
        result = true;
    }

  free (line);
  return result;
}

RXI_STAT rxi_history_save (const char *line, const char *filename)
{
  const char *history_path = rxi_config_path ();
  CHECK (history_path && "Allocation error");
  if (!history_path)
    return RXI_ERR_ALLOC;

  char *history_file = malloc (RXI_PATH_MAX * sizeof (*history_file));
  CHECK (history_file && "Allocation error");
  if (!history_file)
  {
      free ((void*)history_path);
      return RXI_ERR_ALLOC;
  }

  strcpy (history_file, history_path);
  strcat (history_file, filename);

  DEBUG ("Save command '%s' to `%s'", line, history_file);

  int ln_status = linenoiseHistoryAdd (line);
  if (ln_status != 1)
    {
      DEBUG ("Command '%s' cannot be added to history", line);
      free ((void*)history_path);
      free (history_file);
      return RXI_OK;
    }

  ln_status = linenoiseHistorySave (history_file); 
  CHECK ((ln_status == 0) && "linenoise error");

  free ((void*)history_path);
  free (history_file);
  if (ln_status != 0)
    return RXI_ERR_FILE;

  return RXI_OK;
}

RXI_STAT rxi_history_load (const char *filename)
{
  const char *history_path = rxi_config_path ();
  CHECK (history_path && "Allocation error");
  if (!history_path)
    return RXI_ERR_ALLOC;

  char *history_file = malloc (RXI_PATH_MAX * sizeof (*history_file));
  CHECK (history_file && "Allocation error");
  if (!history_file)
  {
      free ((void*)history_path);
      return RXI_ERR_ALLOC;
  }

  strcpy (history_file, history_path);
  strcat (history_file, filename);

  DEBUG ("Load command history from `%s'", history_file);

  linenoiseHistoryReset ();
  int ln_status = linenoiseHistoryLoad (history_file);
  CHECK ((ln_status == 0) && "linenoise error");

  free ((void*)history_path);
  free (history_file);

  if (ln_status != 0)
    return RXI_ERR_FILE;

  return RXI_OK;
}
