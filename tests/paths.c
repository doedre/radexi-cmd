#include <stdlib.h>
#include <stdio.h>

#include "rxi_common.h"

int main(void)
{
  for (int i = 0; i < 1000000; i++)
    {
      const char *db_path = rxi_database_path ();
      const char *config_path = rxi_config_path ();

      printf ("Database: %s\nConfig: %s\n", db_path, config_path);

      free ((void*)db_path);
      free ((void*)config_path);
      printf ("Database: %s\nConfig: %s\n", db_path, config_path);
    }
  return 0;
}
