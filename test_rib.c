#include <assert.h>
#include <stdlib.h>

#include "cloud.h"
#include "comm.h"
#include "loop.h"
#include "vtk.h"

static struct cloud* make_shuffled_cloud(unsigned n)
{
  struct cloud* c = new_cloud(n);
  double* coords = LOOP_MALLOC(double, 3 * n);
  for (unsigned i = 0; i < n; ++i) {
    coords[i * 3 + 0] = (double) i + (1.0 / comm_size()) * comm_rank();
    coords[i * 3 + 1] = 0;
    coords[i * 3 + 2] = 0;
  }
  cloud_add_tag(c, TAG_F64, "coordinates", 3, coords);
  return c;
}

int main(int argc, char** argv)
{
  comm_init();
  assert(argc == 2);
  unsigned n = (unsigned) atoi(argv[1]);
  struct cloud* c = make_shuffled_cloud(n);
  write_parallel_vtu_cloud(c, "before.pvtu");
  free_cloud(c);
  comm_fini();
}
