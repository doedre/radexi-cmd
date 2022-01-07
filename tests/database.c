#include <stdio.h>

#include "utils/database.h"

#include "rxi_common.h"

int main (void)
{
  RXI_STAT status = rxi_add_molecule ("oa", "tests/oatom.dat");
  printf ("Status %u", status);
  return 0;
}
