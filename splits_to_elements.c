#include "splits_to_elements.h"

#include <assert.h>

#include "ints.h"
#include "loop.h"
#include "mesh.h"
#include "tables.h"

void project_splits_to_elements(
    unsigned elem_dim,
    unsigned src_dim,
    unsigned nelems,
    unsigned const* srcs_of_elems,
    unsigned const* gen_offset_of_srcs,
    unsigned const* gen_vert_of_srcs,
    unsigned** gen_offset_of_elems_out,
    unsigned** gen_direction_of_elems_out,
    unsigned** gen_vert_of_elems_out)
{
  assert(elem_dim >= src_dim);
  unsigned srcs_per_elem = the_down_degrees[elem_dim][src_dim];
  unsigned* elem_will_split = uints_filled(nelems, 0);
  unsigned* gen_direction_of_elems = LOOP_MALLOC(unsigned, nelems);
  unsigned* gen_vert_of_elems = LOOP_MALLOC(unsigned, nelems);
  for (unsigned i = 0; i < nelems; ++i) {
    unsigned const* srcs_of_elem = srcs_of_elems + i * srcs_per_elem;
    for (unsigned j = 0; j < srcs_per_elem; ++j) {
      unsigned src = srcs_of_elem[j];
      if (gen_offset_of_srcs[src] == gen_offset_of_srcs[src + 1])
        continue;
      elem_will_split[i] = 1;
      gen_direction_of_elems[i] = j;
      gen_vert_of_elems[i] = gen_vert_of_srcs[src];
      break;
    }
  }
  *gen_offset_of_elems_out = uints_exscan(elem_will_split, nelems);
  *gen_direction_of_elems_out = gen_direction_of_elems;
  *gen_vert_of_elems_out = gen_vert_of_elems;
  loop_free(elem_will_split);
}

void mesh_splits_to_elements(
  struct mesh* m,
  unsigned src_dim,
  unsigned const* gen_offset_of_srcs,
  unsigned const* gen_vert_of_srcs,
  unsigned** p_gen_offset_of_elems,
  unsigned** p_gen_direction_of_elems,
  unsigned** p_gen_vert_of_elems)
{
  unsigned elem_dim = mesh_dim(m);
  unsigned nelems = mesh_count(m, elem_dim);
  unsigned const* srcs_of_elems = mesh_ask_down(m, elem_dim, src_dim);
  project_splits_to_elements(elem_dim, src_dim, nelems,
      srcs_of_elems, gen_offset_of_srcs, gen_vert_of_srcs,
      p_gen_offset_of_elems, p_gen_direction_of_elems, p_gen_vert_of_elems);
}
