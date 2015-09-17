#include "inertia.h"
#include <assert.h>
#include <math.h>

int main()
{
  double A[3][3] = {
    {0.8147,    0.9134,    0.2785},
    {0.9058,    0.6324,    0.5469},
    {0.1270,    0.0975,    0.9575}
  };
  double kl[3] = {-0.1879, 1.7527, 0.8399};
  double l[3];
  unsigned n = eigenvals_3x3(A, l);
  assert(n == 3);
  for (unsigned i = 0; i < 3; ++i) {
    unsigned j = 0;
    for (j = 0; j < 3; ++j)
      if (fabs(l[j] - kl[i]) < 1e-3)
        break;
    assert(j != 3);
  }
}
