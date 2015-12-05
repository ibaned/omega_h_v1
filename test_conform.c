#include <assert.h>

#include "bcast.h"
#include "doubles.h"
#include "comm.h"
#include "global.h"
#include "loop.h"
#include "mesh.h"
#include "migrate_mesh.h"
#include "parallel_mesh.h"
#include "vtk.h"

static struct mesh* make_2_tri_parallel(void)
{
  struct mesh* m = 0;
  assert(comm_size() == 2);
  if (comm_rank() == 0) {
    m = new_box_mesh(2);
    mesh_number_simply(m);
  }
  m = bcast_mesh_metadata(m);
  if (comm_rank() == 0) {
    unsigned n = 1;
    unsigned recvd_elem_ranks[1] = {0};
    unsigned recvd_elem_ids[1] = {0};
    migrate_mesh(&m, n, recvd_elem_ranks, recvd_elem_ids);
  } else {
    unsigned n = 1;
    unsigned recvd_elem_ranks[1] = {0};
    unsigned recvd_elem_ids[1] = {1};
    migrate_mesh(&m, n, recvd_elem_ranks, recvd_elem_ids);
  }
  return m;
}

int main()
{
  comm_init();
  struct mesh* m = make_2_tri_parallel();
  double* data = doubles_filled(3, (double) (comm_rank() + 1));
  mesh_add_tag(m, 0, TAG_F64, "field", 1, data);
  write_parallel_vtu(m, "one.pvtu");
  mesh_accumulate_tag(m, 0, "field");
  write_parallel_vtu(m, "two.pvtu");
  mesh_conform_tag(m, 0, "field");
  write_parallel_vtu(m, "three.pvtu");
  free_mesh(m);
  comm_fini();
}
