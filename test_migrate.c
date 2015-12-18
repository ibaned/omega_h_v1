#include <assert.h>

#include "bcast.h"
#include "comm.h"
#include "ints.h"
#include "loop.h"
#include "mesh.h"
#include "migrate_mesh.h"
#include "parallel_mesh.h"
#include "subset.h"
#include "vtk.h"

int main()
{
  comm_init();
  assert(comm_size() == 2);
  struct mesh* m = 0;
  if (comm_rank() == 0) {
    m = new_box_mesh(2);
    mesh_count(m, 1); //trigger edge creation
  }
  m = bcast_mesh_metadata(m);
  mesh_number_simply(m, 0);
  mesh_number_simply(m, 1);
  write_parallel_vtu(m, "before.pvtu");
  if (comm_rank() == 0) {
    unsigned n = 1;
    unsigned recvd_elem_ranks[1] = {0};
    unsigned recvd_elem_ids[1] = {0};
    migrate_mesh(&m, n, recvd_elem_ranks, recvd_elem_ids);
  } else {
    unsigned n = 2;
    unsigned recvd_elem_ranks[2] = {0,0};
    unsigned recvd_elem_ids[2] = {1,0};
    migrate_mesh(&m, n, recvd_elem_ranks, recvd_elem_ids);
  }
  unsigned* offsets = uints_linear(mesh_count(m, 1) + 1, 1);
  struct mesh* sm = subset_mesh(m, 1, offsets);
  loop_free(offsets);
  free_mesh(m);
  write_parallel_vtu(sm, "after.pvtu");
  free_mesh(sm);
  comm_fini();
}
