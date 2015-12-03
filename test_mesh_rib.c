#include "comm.h"
#include "global.h"
#include "mesh.h"
#include "parallel_inertial_bisect.h"
#include "parallel_mesh.h"
#include "refine_by_size.h"
#include "vtk.h"

int main()
{
  comm_init();
  struct mesh* m = 0;
  if (comm_rank() == 0) {
    m = new_box_mesh(3);
  //for (unsigned i = 0; i < 5; ++i)
  //  uniformly_refine(&m);
  } else
    m = new_empty_mesh(3);
  mesh_number_simply(m);
  write_parallel_vtu(m, "before.pvtu");
  balance_mesh_inertial(&m);
  mesh_global_renumber(m, 0);
  write_parallel_vtu(m, "after.pvtu");
  free_mesh(m);
  comm_fini();
}
