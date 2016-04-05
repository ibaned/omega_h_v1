#include "parallel_to_serial.hpp"

#include "arrays.hpp"
#include "comm.hpp"
#include "loop.hpp"
#include "mesh.hpp"
#include "parallel_mesh.hpp"
#include "exchanger.hpp"
#include "migrate_mesh.hpp"

namespace omega_h {

LOOP_KERNEL(copy_globals,
    unsigned long const* globals,
    unsigned* dest_locals)
  dest_locals[i] = static_cast<unsigned>(globals[i]);
}

void parallel_to_serial(struct mesh* m)
{
  unsigned dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, dim);
  unsigned long total_elems = comm_add_ulong(nelems);
  unsigned long const* globals = mesh_ask_globals(m, dim);
  unsigned* dest_locals = LOOP_MALLOC(unsigned, nelems);
  LOOP_EXEC(copy_globals, nelems, globals, dest_locals);
  unsigned* dest_ranks = filled_array<unsigned>(nelems, 0);
  struct exchanger* elem_push = new_exchanger(nelems, dest_ranks);
  loop_free(dest_ranks);
  unsigned ndests = 0;
  if (comm_rank() == 0)
    ndests = static_cast<unsigned>(total_elems);
  set_exchanger_dests(elem_push, ndests, dest_locals);
  loop_free(dest_locals);
  migrate_mesh(m, elem_push);
  free_exchanger(elem_push);
}

}
