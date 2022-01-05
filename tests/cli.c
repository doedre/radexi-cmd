#include <stdio.h>

#include "utils/cli_tools.h"
#include "linenoise/linenoise.h"

int main (void)
{
  char *line;

//  linenoiseHistorySetMaxLen(10);
  rxi_history_load ("./.history/test1.txt");
  while ((line = rxi_readline ("1> ")) != NULL)
    {
      if (line[0] == '\0')
        break;

      printf ("%s\n", line);
      rxi_history_save (line, "./.history/test1.txt");
    }

  rxi_history_load ("./.history/test2.txt");
  while ((line = rxi_readline ("2> ")) != NULL)
    {
      if (line[0] == '\0')
        break;

      printf ("%s\n", line);
      rxi_history_save (line, "./.history/test2.txt");
    }


  while (!rxi_readline_accept ())
    {
      printf ("Exit?\n");
    }

  return 0;
}
