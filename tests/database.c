#include <stdio.h>

#include "utils/database.h"

#include "rxi_common.h"
#include "utils/debug.h"

int main (void)
{
  RXI_STAT status = rxi_add_molecule ("oa", "tests/oatom.dat");
  ASSERT (status == RXI_OK);
  status = rxi_add_molecule ("oa_delete", "tests/oatom.dat");
  ASSERT (status == RXI_OK);
  status = rxi_delete_molecule ("oa_delete");
  ASSERT (status == RXI_OK);

  return 0;
}
