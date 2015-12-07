#include "ghost_mesh.h"

#include "arrays.h"
#include "exchanger.h"
#include "loop.h"
#include "mesh.h"
#include "parallel_mesh.h"

static void get_elem_use_owners_of_verts(
    struct mesh* m,
    unsigned** p_use_own_ranks,
    unsigned** p_use_own_ids,
    unsigned** p_uses_of_verts_offsets)
{
  unsigned dim = mesh_dim(m);
  unsigned const* elems_of_verts_offsets =
    mesh_ask_up(m, 0, dim)->offsets;
  unsigned const* elems_of_verts =
    mesh_ask_up(m, 0, dim)->adj;
  unsigned const* elem_own_ranks = mesh_ask_own_ranks(m, dim);
  unsigned const* elem_own_ids = mesh_ask_own_ids(m, dim);
  unsigned nverts = mesh_count(m, 0);
  unsigned nuses_in = elems_of_verts_offsets[nverts];
  unsigned* use_ranks_in = LOOP_MALLOC(unsigned, nuses_in);
  unsigned* use_ids_in = LOOP_MALLOC(unsigned, nuses_in);
  for (unsigned i = 0; i < nverts; ++i) {
    unsigned f = elems_of_verts_offsets[i];
    unsigned e = elems_of_verts_offsets[i + 1];
    for (unsigned j = f; j < e; ++j) {
      unsigned elem = elems_of_verts[j];
      use_ranks_in[j] = elem_own_ranks[elem];
      use_ids_in[j] = elem_own_ids[elem];
    }
  }
  unsigned const* vert_own_ranks = mesh_ask_own_ranks(m, 0);
  unsigned const* vert_own_ids = mesh_ask_own_ids(m, 0);
  unsigned* dest_ranks = uints_expand(nverts,
      vert_own_ranks, 1, elems_of_verts_offsets);
  unsigned* dest_ids = uints_expand(nverts,
      vert_own_ids, 1, elems_of_verts_offsets);
  struct exchanger* ex = new_exchanger(nuses_in, dest_ranks);
  loop_free(dest_ranks);
  set_exchanger_dests(ex, nverts, dest_ids);
  loop_free(dest_ids);
  *p_use_own_ranks = exchange_uints(ex, 1, use_ranks_in,
      EX_FOR, EX_ITEM);
  loop_free(use_ranks_in);
  *p_use_own_ids = exchange_uints(ex, 1, use_ids_in,
      EX_FOR, EX_ITEM);
  loop_free(use_ids_in);
  *p_uses_of_verts_offsets = uints_copy(
      ex->items_of_roots_offsets[EX_REV], nverts + 1);
  free_exchanger(ex);
}

void ghost_mesh(struct mesh** p_m, unsigned nlayers)
{
  struct mesh* m = *p_m;
  unsigned dim = mesh_dim(m);
  /* we assume the input mesh has no ghosting, so we
     can "simply" number the element to determine
     their ownership */
  mesh_number_simply(m, dim);
  (void) nlayers;
  /* get use owners of elements used by vertices */
  /* close_partition to determine 1st layer elements */
  /* 1st layer element used vertex owners, pulled from
     input layer used vertex owners */
  /* close_partition to determine 1st layer vertices */
}
