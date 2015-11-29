#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "cloud.h"
#include "comm.h"
#include "doubles.h"
#include "ints.h"
#include "loop.h"
#include "migrate_cloud.h"
#include "parallel_inertial_bisect.h"
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
  double* coords = doubles_copy(
      cloud_find_tag(c, "coordinates")->d.f64, n * 3);
  unsigned* orig_ranks = uints_filled(n, comm_rank());
  unsigned* ones = uints_filled(n, 1);
  unsigned* orig_ids = uints_exscan(ones, n);
  loop_free(ones);
  recursive_inertial_bisect(&n, &coords, 0, &orig_ranks, &orig_ids);
  loop_free(coords);
  migrate_cloud(&c, n, orig_ranks, orig_ids);
  loop_free(orig_ranks);
  loop_free(orig_ids);
  write_parallel_vtu_cloud(c, "after.pvtu");
  free_cloud(c);
  comm_fini();
}
