#include <assert.h>

#include "bcast.h"
#include "comm.h"
#include "global.h"
#include "loop.h"
#include "mesh.h"
#include "migrate_mesh.h"
#include "vtk.h"

int main()
{
  comm_init();
  assert(comm_size() == 2);
  struct mesh* m = 0;
  if (comm_rank() == 0) {
    m = new_box_mesh(2);
    mesh_number_simply(m);
  }
  m = bcast_mesh_metadata(m);
  write_parallel_vtu(m, "before.pvtu");
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
  write_parallel_vtu(m, "after.pvtu");
  free_mesh(m);
  comm_fini();
}