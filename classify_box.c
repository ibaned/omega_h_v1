#include "classify_box.h"
#include <math.h>
#include <stdlib.h>

#define EPSILON 1e-10

unsigned* classify_box(
    unsigned elem_dim,
    unsigned nverts,
    double const* coords)
{
  unsigned* out = malloc(sizeof(unsigned) * nverts);
  for (unsigned i = 0; i < nverts; ++i) {
    double const* x = coords + i * 3;
    unsigned classif_dim = 3;
    for (unsigned j = 0; j < 3; ++j)
      if ((fabs(x[j] - 0) < EPSILON) ||
          (fabs(x[j] - 1) < EPSILON))
        --classif_dim;
    out[i] = classif_dim;
  }
  return out;
}
