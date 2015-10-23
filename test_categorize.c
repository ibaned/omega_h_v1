#include <stdio.h>

#include "loop.h"
#include "global.h"

int main()
{
  unsigned const a[6] = {7,7,2,2,13,13};
  unsigned* cats;
  unsigned* cat_offsets;
  categorize_by_part(a, 6, &cats, &cat_offsets);
  for (unsigned i = 0; i < 6; ++i)
    printf("%u %u\n", cats[i], cat_offsets[i]);
  loop_free(cats);
  loop_free(cat_offsets);
}
