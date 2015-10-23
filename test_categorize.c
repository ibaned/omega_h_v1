#include <stdio.h>

#include "loop.h"
#include "global.h"

int main()
{
  unsigned const a[8] = {7,13,2,7,2,2,13,13};
  unsigned n = 8;
  unsigned* cats;
  unsigned* cat_indices;
  unsigned ncats;
  unsigned* cat_parts;
  unsigned* cat_counts;
  categorize_by_part(a, n, &cats, &cat_indices,
      &ncats, &cat_parts, &cat_counts);
  printf("ncats %u\n", ncats);
  for (unsigned i = 0; i < ncats; ++i)
    printf("%u %u\n", cat_parts[i], cat_counts[i]);
  loop_free(cat_parts);
  loop_free(cat_counts);
  printf("\n");
  for (unsigned i = 0; i < n; ++i)
    printf("%u %u\n", cats[i], cat_indices[i]);
  loop_free(cats);
  loop_free(cat_indices);
}
