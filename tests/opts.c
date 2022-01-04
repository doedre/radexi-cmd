#include <stdlib.h>
#include <stdio.h>

#include "rxi_common.h"
#include "utils/options.h"

int
main (int argc, char **argv)
{
  struct rxi_options opts;

  int index = rxi_set_options (&opts, argc, argv);
  printf ("Optind: %d\n", index);

	exit (EXIT_SUCCESS);
}
