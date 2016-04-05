#include "parallel_to_serial.hpp"

#include "arrays.hpp"
#include "comm.hpp"
#include "loop.hpp"
#include "ints.hpp"
#include "mesh.hpp"
#include "parallel_mesh.hpp"
#include "exchanger.hpp"
#include "migrate_mesh.hpp"
#include "reorder_mesh.hpp"

namespace omega_h {

LOOP_KERNEL(copy_globals,
    unsigned long const* globals,
    unsigned* out)
  out[i] = static_cast<unsigned>(globals[i]);
}

static unsigned* old_to_new_by_global(struct mesh* m, unsigned dim)
{
  unsigned long const* globals = mesh_ask_globals(m, dim);
  unsigned* out = LOOP_MALLOC(unsigned, mesh_count(m, dim));
  LOOP_EXEC(copy_globals, mesh_count(m, dim), globals, out);
  return out;
}

/* migration brings the elements in the right order but
   not the other entities.
   we'll just reorder them. */
static void reorder_lower_dims(struct mesh* m)
{
  unsigned dim = mesh_dim(m);
  unsigned const* old_to_new_ents[4] = {0,0,0,0};
  unsigned* to_free[4] = {0,0,0,0};
  for (unsigned d = 0; d <= dim; ++d) {
    if (!mesh_has_dim(m, d))
      continue;
    old_to_new_ents[d] = to_free[d] = old_to_new_by_global(m, d);
  }
  reorder_mesh(m, old_to_new_ents);
  for (unsigned d = 0; d <= dim; ++d)
    loop_free(to_free[d]);
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
  unsigned* sent_of_srcs_offsets = uints_linear(nelems + 1, 1);
  set_exchanger_srcs(elem_push, nelems, sent_of_srcs_offsets);
  loop_free(sent_of_srcs_offsets);
  migrate_mesh(m, elem_push);
  free_exchanger(elem_push);
  reorder_lower_dims(m);
}

}
