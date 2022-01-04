#include <stdio.h>

#include "utils/debug.h"

int main ()
{
  DEBUG ("%s, %d", "string", 4);
  CHECK (0);
  ASSERT (0);
  return 0;
}
